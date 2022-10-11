/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                    W. Liu                                                  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file algo_aitken.h
 * @author W. Liu
 * @date 05 October 2022
 * @brief Aitken coupling algorithm
 */

#ifndef MUI_COUPLING_ALGORITHM_AITKEN_H_
#define MUI_COUPLING_ALGORITHM_AITKEN_H_

#include "../util.h"
#include "../config.h"
#include <iterator>

namespace mui {

template<typename CONFIG=default_config> class algo_aitken {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	using point_type = typename CONFIG::point_type;

	algo_aitken( REAL undRelxFac = 1.0,
 	  REAL undRelxFacMax = 1.0,
 	  std::vector<std::pair<point_type, REAL>> ptsVluInit = 
		std::vector<std::pair<point_type, REAL>>(),
 	  REAL resL2NormNM1 = 0.0):
	  initUndRelxFac_(undRelxFac),
	  undRelxFacMax_(undRelxFacMax){

		undRelxFac_.insert(undRelxFac_.begin(),
			std::make_pair(std::make_pair(
				std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest()), initUndRelxFac_
			)
		);

		if (!ptsVluInit.empty()) {
			ptsTimeVlu_.insert(ptsTimeVlu_.begin(),
				std::make_pair(std::make_pair(
					std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest()),ptsVluInit
				)
			);
		}

		if(resL2NormNM1 != 0.0){
			residualL2Norm_.insert(residualL2Norm_.begin(),
				std::make_pair(std::make_pair(
					std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest()), 
					(std::make_pair(static_cast<INT>(0), resL2NormNM1)
					)
				)
			);
		}
	}

	//- relaxation based on single time value
	template<typename OTYPE>
	OTYPE relaxation(std::pair<time_type,time_type> t, point_type focus, OTYPE filteredValue) const {

		OTYPE filteredOldValue = 0.0;

		if (ptsTimeVlu_.empty()) {

			assert(ptsTimeRes_.empty());

			filteredOldValue = 0.0;

			std::vector<std::pair<point_type, REAL>> ptsVluTemp{
				std::make_pair(focus,calculate_relaxed_value(t,filteredValue,filteredOldValue))
				};

			ptsTimeVlu_.insert(ptsTimeVlu_.begin(),std::make_pair(t,ptsVluTemp));

			std::vector<std::pair<point_type, REAL>> ptsResTemp{
				std::make_pair(focus,calculate_point_residual(t,filteredValue,filteredOldValue))
				};

			ptsTimeRes_.insert(ptsTimeRes_.begin(),std::make_pair(t,ptsResTemp));

			return calculate_relaxed_value(t,filteredValue,filteredOldValue);

		} else {

				auto presentIter = std::find_if(ptsTimeVlu_.begin(), ptsTimeVlu_.end(),
					[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
						return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

				auto previousIter = std::find_if(ptsTimeVlu_.begin(), ptsTimeVlu_.end(),
					[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
						return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});

				auto presentResIter = std::find_if(ptsTimeRes_.begin(), ptsTimeRes_.end(),
					[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
						return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

				auto previousResIter = std::find_if(ptsTimeRes_.begin(), ptsTimeRes_.end(),
					[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
						return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});

				if ((presentIter == std::end(ptsTimeVlu_)) &&
					(previousIter == std::end(ptsTimeVlu_)) ) {

					assert((presentResIter == std::end(ptsTimeRes_)) || ptsTimeRes_.empty());

					std::cerr << "Non-monotonic time marching does not (yet) supported for the Aitken coupling method! " << std::endl;

				} else if ((presentIter != std::end(ptsTimeVlu_)) &&
					(previousIter == std::end(ptsTimeVlu_)) ) {

						assert((((presentIter->first.first - presentResIter->first.first) < std::numeric_limits<REAL>::epsilon()) && 
								((presentIter->first.second - presentResIter->first.second) < std::numeric_limits<REAL>::epsilon())) || ptsTimeRes_.empty());

						if (!ptsTimeRes_.empty()) {
							assert(!residualL2Norm_.empty());
						}

						auto ptsRelxValIter = std::find_if(presentIter->second.begin(),
							presentIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						auto ptsRelxResIter = std::find_if(presentResIter->second.begin(),
							presentResIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						if ( ptsRelxValIter == std::end(presentIter->second) ) {

							assert((ptsRelxResIter == std::end(presentResIter->second)) || ptsTimeRes_.empty());

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < presentIter->second.size() ; i++ ) {
								REAL dr2 = normsq( focus - presentIter->second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = presentIter->second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = presentIter->second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							auto relaxedValueTemp = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

							presentIter->second.insert(presentIter->second.begin(),std::make_pair(focus, relaxedValueTemp));

							if (ptsTimeRes_.empty()) {

								std::vector<std::pair<point_type, REAL>> ptsResTemp{
									std::make_pair(focus, (filteredValue-relaxedValueTemp))
									};

								ptsTimeRes_.insert(ptsTimeRes_.begin(),std::make_pair(t,ptsResTemp));

							} else {
								presentResIter->second.insert(presentResIter->second.begin(),std::make_pair(focus, (filteredValue-relaxedValueTemp)));
							}

							return relaxedValueTemp;

						} else {

							assert((normsq(ptsRelxValIter->first - ptsRelxResIter->first) < std::numeric_limits<REAL>::epsilon()) || ptsTimeRes_.empty());

							return ptsRelxValIter->second;

						}

				} else if ((presentIter == std::end(ptsTimeVlu_)) &&
					(previousIter != std::end(ptsTimeVlu_)) ) {

						assert((presentResIter == std::end(ptsTimeRes_)) || ptsTimeRes_.empty());

						auto ptsRelxValIter = std::find_if(previousIter->second.begin(),
							previousIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						if ( ptsRelxValIter == std::end(previousIter->second) ) {

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < previousIter->second.size() ; i++ ) {
								REAL dr2 = normsq( focus - previousIter->second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = previousIter->second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = previousIter->second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							filteredOldValue = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

						} else {

							filteredOldValue = ptsRelxValIter->second;

						}

						std::vector<std::pair<point_type, REAL>> ptsVluTemp{
							std::make_pair(focus,calculate_relaxed_value(t,filteredValue,filteredOldValue))
						};

						ptsTimeVlu_.insert(ptsTimeVlu_.begin(),std::make_pair(t, ptsVluTemp));

						std::vector<std::pair<point_type, REAL>> ptsResTemp{
							std::make_pair(focus,calculate_point_residual(t,filteredValue,filteredOldValue))
							};

						ptsTimeRes_.insert(ptsTimeRes_.begin(),std::make_pair(t,ptsResTemp));

						presentIter = std::find_if(ptsTimeVlu_.begin(), ptsTimeVlu_.end(),
							[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
								return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
										((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

						previousIter = std::find_if(ptsTimeVlu_.begin(), ptsTimeVlu_.end(),
							[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
								return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
									(b.first.first < t.first);});

						presentResIter = std::find_if(ptsTimeRes_.begin(), ptsTimeRes_.end(),
							[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
								return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
										((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

						previousResIter = std::find_if(ptsTimeRes_.begin(), ptsTimeRes_.end(),
							[t](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
							return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
									(b.first.first < t.first);});

						auto ptsResidualL2NormIter = std::find_if(residualL2Norm_.begin(),
							residualL2Norm_.end(), [previousIter](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
						return ((previousIter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
								((previousIter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

						if(ptsResidualL2NormIter == std::end(residualL2Norm_)) {

							if (previousResIter!=std::end(ptsTimeRes_)) {

								assert(((previousIter->first.first - previousResIter->first.first) < std::numeric_limits<REAL>::epsilon()) && 
										((previousIter->first.second - previousResIter->first.second) < std::numeric_limits<REAL>::epsilon()) );

								REAL localResidualMagSqSumTemp = 0.0;

								for (auto & elementPair : previousResIter->second) {
									localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
								}

								if((localResidualMagSqSumTemp != 0) || (!residualL2Norm_.empty())){
									residualL2Norm_.insert(residualL2Norm_.begin(),
										std::make_pair(
											previousResIter->first, (
												std::make_pair(
													static_cast<INT>(
														previousResIter->second.size()
													), std::sqrt(localResidualMagSqSumTemp)
												)
											)
										)
									);

									ptsResidualL2NormIter = std::find_if(residualL2Norm_.begin(),
										residualL2Norm_.end(), [previousIter](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
									return ((previousIter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
											((previousIter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});
								}
							}

						} else {

							if(ptsResidualL2NormIter->second.first != 0) {

								auto ptsTimeResIter = std::find_if(ptsTimeRes_.begin(),
								ptsTimeRes_.end(), [previousResIter](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
								return ((previousResIter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
										((previousResIter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

								assert(ptsTimeResIter != std::end(ptsTimeRes_) );

								if(ptsTimeResIter->second.size() != ptsResidualL2NormIter->second.first) {

									REAL localResidualMagSqSumTemp = 0.0;

									for (auto & elementPair : ptsTimeResIter->second) {
										localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
									}
									ptsResidualL2NormIter->second.second = std::sqrt(localResidualMagSqSumTemp);
								}
							}
						}

						return calculate_relaxed_value(t,filteredValue,filteredOldValue);

				} else {

						assert((((presentIter->first.first - presentResIter->first.first) < std::numeric_limits<REAL>::epsilon()) && 
								((presentIter->first.second - presentResIter->first.second) < std::numeric_limits<REAL>::epsilon())) || ptsTimeRes_.empty());

						auto ptsRelxValIter = std::find_if(previousIter->second.begin(),
							previousIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						if ( ptsRelxValIter == std::end(previousIter->second) ) {

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < previousIter->second.size() ; i++ ) {
								REAL dr2 = normsq( focus - previousIter->second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = previousIter->second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = previousIter->second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							filteredOldValue = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

						} else {

							filteredOldValue = ptsRelxValIter->second;

						}

						auto ptsPresentRelxValIter = std::find_if(presentIter->second.begin(),
							presentIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						auto ptsRelxResIter = std::find_if(presentResIter->second.begin(),
							presentResIter->second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						});

						if ( ptsPresentRelxValIter == std::end(presentIter->second) ) {

							assert((ptsRelxResIter == std::end(presentResIter->second)) || ptsTimeRes_.empty());

							presentIter->second.insert(presentIter->second.begin(),
								std::make_pair(focus,calculate_relaxed_value(
									t,filteredValue,filteredOldValue)
								)
							);

							if (ptsTimeRes_.empty()) {

								std::vector<std::pair<point_type, REAL>> ptsResTemp{
									std::make_pair(focus,calculate_point_residual(
										t,filteredValue,filteredOldValue))
									};

								ptsTimeRes_.insert(ptsTimeRes_.begin(),std::make_pair(t,ptsResTemp));

							} else {

								presentResIter->second.insert(presentResIter->second.begin(),
									std::make_pair(focus,calculate_point_residual(
										t,filteredValue,filteredOldValue)
									)
								);
							}
						} else {

							assert((normsq(ptsPresentRelxValIter->first - ptsRelxResIter->first) < std::numeric_limits<REAL>::epsilon()) || ptsTimeRes_.empty());

							ptsPresentRelxValIter->second = calculate_relaxed_value(t,filteredValue,filteredOldValue);

							if (ptsTimeRes_.empty()) {

								std::vector<std::pair<point_type, REAL>> ptsResTemp{
									std::make_pair(focus,calculate_point_residual(
										t,filteredValue,filteredOldValue))
									};

								ptsTimeRes_.insert(ptsTimeRes_.begin(),std::make_pair(t,ptsResTemp));

							} else {
								ptsRelxResIter->second = calculate_point_residual(t,filteredValue,filteredOldValue);
							}
						}

						return calculate_relaxed_value(t,filteredValue,filteredOldValue);

				}

		}
	}

	REAL get_under_relaxation_factor (time_type t_single) {

		std::pair<time_type,time_type> t = std::make_pair(t_single,std::numeric_limits<time_type>::lowest());

		update_undRelxFac(t);

		auto undRelxPresentIter = std::find_if(undRelxFac_.begin(),
			undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
			return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
			((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

		assert(undRelxPresentIter != std::end(undRelxFac_) );

		return undRelxPresentIter->second;

	}

	REAL get_under_relaxation_factor (std::pair<time_type,time_type> t) {

		update_undRelxFac(t);

		auto undRelxPresentIter = std::find_if(undRelxFac_.begin(),
			undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
			return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
			((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

		assert(undRelxPresentIter != std::end(undRelxFac_) );

		return undRelxPresentIter->second;

	}

	REAL get_residual_L2_Norm (time_type t_single) {

		std::pair<time_type,time_type> t = std::make_pair(t_single,std::numeric_limits<time_type>::lowest());

		update_undRelxFac(t);

		auto resL2NormNM1Iter = std::find_if(residualL2Norm_.begin(),
			residualL2Norm_.end(), [t](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
		return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
			(b.first.first < t.first);});

		if(resL2NormNM1Iter != std::end(residualL2Norm_) ) {
			return resL2NormNM1Iter->second.second;
		} else {
			return 0.0;
		}
	}

	REAL get_residual_L2_Norm (std::pair<time_type,time_type> t) {

		update_undRelxFac(t);

		auto resL2NormNM1Iter = std::find_if(residualL2Norm_.begin(),
			residualL2Norm_.end(), [t](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
		return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
			(b.first.first < t.first);});

		if(resL2NormNM1Iter != std::end(residualL2Norm_) ) {
			return resL2NormNM1Iter->second.second;
		} else {
			return 0.0;
		}
	}

private:
	template<typename OTYPE>
	OTYPE calculate_relaxed_value(std::pair<time_type,time_type> t, OTYPE filteredValue, OTYPE filteredOldValue) {
 
		update_undRelxFac(t);

		auto undRelxPresentIter = std::find_if(undRelxFac_.begin(),
			undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
			return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
			((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

		assert(undRelxPresentIter != std::end(undRelxFac_) );

		return (undRelxPresentIter->second * filteredValue) + 
			((1 - undRelxPresentIter->second) * filteredOldValue);

	}

	void update_undRelxFac(std::pair<time_type,time_type> t) {

		auto undRelxPresentIter = std::find_if(undRelxFac_.begin(),
			undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
		return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
			((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

		auto undRelxPrevIter = std::find_if(undRelxFac_.begin(),
			undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
		return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
			(b.first.first < t.first);});

		auto resL2NormNM1Iter = std::find_if(residualL2Norm_.begin(),
			residualL2Norm_.end(), [t](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
		return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
			(b.first.first < t.first);});

		auto resL2NormNM2Iter = std::find_if(residualL2Norm_.begin(),
			residualL2Norm_.end(), [resL2NormNM1Iter](std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>> b) {
		return (((b.first.first - resL2NormNM1Iter->first.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < resL2NormNM1Iter->first.second)) || 
				(b.first.first < resL2NormNM1Iter->first.first);});

		if (undRelxPresentIter==std::end(undRelxFac_)) {

			if(resL2NormNM2Iter == std::end(residualL2Norm_) ) {
				undRelxFac_.insert(undRelxFac_.begin(),std::make_pair(t, initUndRelxFac_));

				undRelxPresentIter = std::find_if(undRelxFac_.begin(),
					undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
				return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
						((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

				undRelxPrevIter = std::find_if(undRelxFac_.begin(),
					undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
				return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});

			} else {

				assert(resL2NormNM1Iter != std::end(residualL2Norm_) );

				if(resL2NormNM2Iter->second.first != 0 ) {
					auto ptsTimeResNM2Iter = std::find_if(ptsTimeRes_.begin(),
					ptsTimeRes_.end(), [resL2NormNM2Iter](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
					return ((resL2NormNM2Iter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((resL2NormNM2Iter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

					assert(ptsTimeResNM2Iter != std::end(ptsTimeRes_) );

					if(ptsTimeResNM2Iter->second.size() != resL2NormNM2Iter->second.first) {

						REAL localResidualMagSqSumTemp = 0.0;

						for (auto & elementPair : ptsTimeResNM2Iter->second) {
							localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
						}

						resL2NormNM2Iter->second.second = std::sqrt(localResidualMagSqSumTemp);

						update_undRelxFac(resL2NormNM2Iter->first);
					}

					auto ptsTimeResNM1Iter = std::find_if(ptsTimeRes_.begin(),
					ptsTimeRes_.end(), [resL2NormNM1Iter](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
					return ((resL2NormNM1Iter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((resL2NormNM1Iter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

					assert(ptsTimeResNM1Iter != std::end(ptsTimeRes_) );

					if(ptsTimeResNM1Iter->second.size() != resL2NormNM1Iter->second.first) {

						REAL localResidualMagSqSumTemp = 0.0;

						for (auto & elementPair : ptsTimeResNM1Iter->second) {
							localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
						}

						resL2NormNM1Iter->second.second = std::sqrt(localResidualMagSqSumTemp);

						update_undRelxFac(resL2NormNM1Iter->first);
					}
				}

				double nominator   = resL2NormNM2Iter->second.second;  
				double denominator = resL2NormNM1Iter->second.second - resL2NormNM2Iter->second.second;

				if (denominator != 0.0 ) {

					if (undRelxPrevIter==std::end(undRelxFac_)) {
						undRelxFac_.insert(undRelxFac_.begin(),
							std::make_pair(
								t, calculate_aitken_constraint_PNControl_zeroControl(initUndRelxFac_)
							)
						);

						undRelxPresentIter = std::find_if(undRelxFac_.begin(),
							undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
						return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());
						});

						undRelxPrevIter = std::find_if(undRelxFac_.begin(),
							undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
						return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});
					} else {
						undRelxFac_.insert(undRelxFac_.begin(),
							std::make_pair(
								t, calculate_aitken_constraint_PNControl_zeroControl(
									-undRelxPrevIter->second * (nominator/denominator)
								)
							)
						);

						undRelxPresentIter = std::find_if(undRelxFac_.begin(),
							undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
						return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

						undRelxPrevIter = std::find_if(undRelxFac_.begin(),
							undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
						return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});
					}

				} else {
					undRelxFac_.insert(undRelxFac_.begin(),std::make_pair(t, initUndRelxFac_));

					undRelxPresentIter = std::find_if(undRelxFac_.begin(),
						undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
					return ((t.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((t.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

					undRelxPrevIter = std::find_if(undRelxFac_.begin(),
						undRelxFac_.end(), [t](std::pair<std::pair<time_type,time_type>, REAL> b) {
					return (((b.first.first - t.first) < std::numeric_limits<REAL>::epsilon()) && (b.first.second < t.second)) || 
							(b.first.first < t.first);});
				}
			}
		} else { 

			if(resL2NormNM2Iter == std::end(residualL2Norm_) ) {
				if(undRelxPresentIter->second != initUndRelxFac_) {
					std::cout << "change under Relx Factor to its initial value." << std::endl;
					undRelxPresentIter->second = initUndRelxFac_;
				}
			} else {
				assert(resL2NormNM1Iter != std::end(residualL2Norm_) );

				if(resL2NormNM2Iter->second.first != 0 ) {
					auto ptsTimeResNM2Iter = std::find_if(ptsTimeRes_.begin(),
					ptsTimeRes_.end(), [resL2NormNM2Iter](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
					return ((resL2NormNM2Iter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((resL2NormNM2Iter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon());});

					assert(ptsTimeResNM2Iter != std::end(ptsTimeRes_) );

					if(ptsTimeResNM2Iter->second.size() != resL2NormNM2Iter->second.first) {

						REAL localResidualMagSqSumTemp = 0.0;

						for (auto & elementPair : ptsTimeResNM2Iter->second) {
							localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
						}

						resL2NormNM2Iter->second.second = std::sqrt(localResidualMagSqSumTemp);

						update_undRelxFac(resL2NormNM2Iter->first);
					}

					auto ptsTimeResNM1Iter = std::find_if(ptsTimeRes_.begin(),
					ptsTimeRes_.end(), [resL2NormNM1Iter](std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>> b) {
					return (((resL2NormNM1Iter->first.first - b.first.first) < std::numeric_limits<REAL>::epsilon()) && 
							((resL2NormNM1Iter->first.second - b.first.second) < std::numeric_limits<REAL>::epsilon()));});

					assert(ptsTimeResNM1Iter != std::end(ptsTimeRes_) );

					if(ptsTimeResNM1Iter->second.size() != resL2NormNM1Iter->second.first) {

						REAL localResidualMagSqSumTemp = 0.0;

						for (auto & elementPair : ptsTimeResNM1Iter->second) {
							localResidualMagSqSumTemp += std::pow(elementPair.second, 2);
						}

						resL2NormNM1Iter->second.second = std::sqrt(localResidualMagSqSumTemp);

						update_undRelxFac(resL2NormNM1Iter->first);
					}
				}

				double nominator   = resL2NormNM2Iter->second.second;  
				double denominator = resL2NormNM1Iter->second.second - resL2NormNM2Iter->second.second;

				if (denominator != 0.0 ) {

					if (undRelxPrevIter==std::end(undRelxFac_)) {
						if(undRelxPresentIter->second != calculate_aitken_constraint_PNControl_zeroControl(initUndRelxFac_)
						) {
							std::cout << "Update under Relx Factor." << std::endl;
							undRelxPresentIter->second = calculate_aitken_constraint_PNControl_zeroControl(initUndRelxFac_);
						}
					} else {
						if(undRelxPresentIter->second != calculate_aitken_constraint_PNControl_zeroControl(
							-undRelxPrevIter->second * (nominator/denominator))
						) {
							std::cout << "Update under Relx Factor." << std::endl;
							undRelxPresentIter->second = calculate_aitken_constraint_PNControl_zeroControl(
								-undRelxPrevIter->second * (nominator/denominator)
							);
						}
					}
				} else {
					if(undRelxPresentIter->second != initUndRelxFac_) {
							std::cout << "Update under Relx Factor." << std::endl;
							undRelxPresentIter->second = initUndRelxFac_;
					}
				}
			}
		}
	}

	REAL calculate_aitken_constraint(REAL undRelxfactor) {

		return (sign(undRelxfactor) * std::min(std::abs(undRelxfactor),undRelxFacMax_));

	}

	REAL calculate_aitken_constraint_PNControl(REAL undRelxfactor) {

		return (std::min(std::abs(undRelxfactor),undRelxFacMax_));
	}

	REAL calculate_aitken_constraint_PNControl_zeroControl(REAL undRelxfactor) {

		return ((std::min(std::abs(undRelxfactor),undRelxFacMax_)) < initUndRelxFac_) ? initUndRelxFac_ : (std::min(std::abs(undRelxfactor),undRelxFacMax_));
	}

	template<typename OTYPE>
	OTYPE calculate_point_residual(std::pair<time_type,time_type> t, OTYPE filteredValue, OTYPE filteredOldValue) {

		return (filteredValue - calculate_relaxed_value(t, filteredValue, filteredOldValue));

	}

	INT sign(REAL value)
	{
		return (value < 0) ? -1 : ((value > 0) ? 1 : 0);
	}

protected:
	const REAL initUndRelxFac_;
	const REAL undRelxFacMax_;
	mutable std::vector<std::pair<std::pair<time_type,time_type>, REAL>> undRelxFac_;
	
	mutable std::vector<std::pair<std::pair<time_type,time_type>,std::pair<INT, REAL>>> residualL2Norm_;

	mutable std::vector<std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>>> ptsTimeVlu_;
	mutable std::vector<std::pair<std::pair<time_type,time_type>,std::vector<std::pair<point_type, REAL>>>> ptsTimeRes_;

};

}

#endif /* MUI_COUPLING_ALGORITHM_AITKEN_H_ */
