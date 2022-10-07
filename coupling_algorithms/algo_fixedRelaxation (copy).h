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
 * @file algo_fixedRelaxation.h
 * @author W. Liu
 * @date 05 October 2022
 * @brief Fixed Relaxation coupling algorithm
 */

#ifndef MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_
#define MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_

#include "../util.h"
#include "../config.h"
#include <iterator>

namespace mui {

template<typename CONFIG=default_config> class algo_fixedRelaxation {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	using point_type = typename CONFIG::point_type;

	algo_fixedRelaxation( REAL undRelxFac = 1.0, 
 	  std::vector<std::pair<point_type, REAL>> ptsVluInit = 
		std::vector<std::pair<point_type, REAL>>() ):
 	  initUndRelxFac_(undRelxFac) {

		undRelxFac_ = initUndRelxFac_;

		if (!ptsVluInit.empty()) {
			ptsTimeVlu_.push_back(
				std::make_pair(
				std::numeric_limits<time_type>::lowest(),ptsVluInit
				)
			);
		}
	}

	//- relaxation based on single time value
	template<typename OTYPE>
	OTYPE relaxation(time_type t, point_type focus, OTYPE filteredValue) const {

		OTYPE filteredOldValue = 0.0;

			if (ptsTimeVlu_.empty()) {
				
				filteredOldValue = 0.0

				std::vector<std::pair<point_type, REAL>> ptsVluTemp{
					std::make_pair(focus,get_relaxed_value(filteredValue,filteredOldValue))
					};

				ptsTimeVlu_.push_back(std::make_pair(t,ptsVluTemp));

				return get_relaxed_value(filteredValue,filteredOldValue);

			} else { // ptsTimeVlu_ not empty

				auto presentIter = std::find_if(ptsTimeVlu_.rbegin(), ptsTimeVlu_.rend(), 
					[t](std::pair<time_type,std::vector<std::pair<point_type, REAL>>> b) {
						return (t - b.first) < std::numeric_limits<REAL>::epsilon();
					}

				auto previousIter = std::find_if(ptsTimeVlu_.rbegin(), ptsTimeVlu_.rend(), 
					[t](std::pair<time_type,std::vector<std::pair<point_type, REAL>>> b) {
						return b.first < t;
					}

				if ((presentIter == std::rend(ptsTimeVlu_)) && 
					(previousIter == std::rend(ptsTimeVlu_)) ) {

					std::cerr << "Non-monotonic time marching does not (yet) supported \
						for the Fixed Relaxation coupling method! " << std::end;

				} else if ((presentIter != std::rend(ptsTimeVlu_)) && 
					(previousIter == std::rend(ptsTimeVlu_)) ) {

						auto ptsRelxValIter = std::find_if(presentIter.second.begin(), 
							presentIter.second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						};

						if ( ptsRelxValIter == std::end(presentIter.second) ) {

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < presentIter.second.size() ; i++ ) {
								REAL dr2 = normsq( focus - presentIter.second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = presentIter.second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = presentIter.second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							auto relaxedValueTemp = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

							presentIter.second.push_back(std::make_pair(focus, relaxedValueTemp));

							return relaxedValueTemp;		

						} else {

							return ptsRelxValIter.second;

						}

				} else if ((presentIter == std::rend(ptsTimeVlu_)) && 
					(previousIter != std::rend(ptsTimeVlu_)) ) {

						auto ptsRelxValIter = std::find_if(previousIter.second.begin(), 
							previousIter.second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						};

						if ( ptsRelxValIter == std::end(previousIter.second) ) {

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < previousIter.second.size() ; i++ ) {
								REAL dr2 = normsq( focus - previousIter.second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = previousIter.second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = previousIter.second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							filteredOldValue = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

						} else {

							filteredOldValue = ptsRelxValIter.second;

						}

						std::vector<std::pair<point_type, REAL>> ptsVluTemp{
							std::make_pair(focus,get_relaxed_value(filteredValue,filteredOldValue))
						};

						ptsTimeVlu_.push_back(std::make_pair(t, ptsVluTemp));

						return get_relaxed_value(filteredValue,filteredOldValue);

				} else {

						auto ptsRelxValIter = std::find_if(previousIter.second.begin(), 
							previousIter.second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						};

						if ( ptsRelxValIter == std::end(previousIter.second) ) {

							// Interpolate the relaxed value by N2_linear
							REAL r2min_1st = std::numeric_limits<REAL>::max();
							REAL r2min_2nd = std::numeric_limits<REAL>::max();
							OTYPE value_1st = 0, value_2nd = 0;
							for( size_t i = 0 ; i < previousIter.second.size() ; i++ ) {
								REAL dr2 = normsq( focus - previousIter.second[i].first );
								if ( dr2 < r2min_1st ) {
									r2min_2nd = r2min_1st;
									value_2nd = value_1st;
									r2min_1st = dr2;
									value_1st = previousIter.second[i].second ;
								} else if ( dr2 < r2min_2nd ) {
									r2min_2nd = dr2;
									value_2nd = previousIter.second[i].second ;
								}
							}

							REAL r1 = std::sqrt( r2min_1st );
							REAL r2 = std::sqrt( r2min_2nd );

							filteredOldValue = ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );

						} else {

							filteredOldValue = ptsRelxValIter.second;

						}

						auto ptsRelxValIter = std::find_if(presentIter.second.begin(), 
							presentIter.second.end(), [focus](std::pair<point_type, REAL> b) {
						return normsq(focus - b.first) < std::numeric_limits<REAL>::epsilon();
						};

						if ( ptsRelxValIter == std::end(presentIter.second) ) {

							presentIter.second.push_back(
								std::make_pair(focus,get_relaxed_value(
									filteredValue,filteredOldValue)
								)
							);

						} else {

							ptsRelxValIter.second = get_relaxed_value(filteredValue,filteredOldValue);

						}

						return get_relaxed_value(filteredValue,filteredOldValue);

				}

			}

private:
	template<typename OTYPE>
	OTYPE get_relaxed_value(OTYPE filteredValue, OTYPE filteredOldValue) {

		return (undRelxFac_ * filteredValue) + ((1 - undRelxFac_) * filteredOldValue);

	}


	}

protected:
	REAL initUndRelxFac_;
	REAL undRelxFac_;

	mutable std::vector<std::pair<time_type,std::vector<std::pair<point_type, REAL>>>> ptsTimeVlu_;

};

}

#endif /* MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_ */
