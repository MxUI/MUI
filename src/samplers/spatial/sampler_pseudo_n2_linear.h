/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
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
 * @file sampler_pseudo_n2_linear.h
 * @author Y. H. Tang
 * @date 10 October 2014
 * @brief Spatial sampler that provides a value at a point using a
 * pseudo-linear n^2 interpolation.
 */

#ifndef MUI_SAMPLER_PSEUDO_N2_LINEAR_H_
#define MUI_SAMPLER_PSEUDO_N2_LINEAR_H_

#include "../../config.h"
#include "../sampler.h"

namespace mui {

template<typename CONFIG=default_config, typename O_TP=typename CONFIG::REAL, typename I_TP=O_TP>
class sampler_pseudo_n2_linear {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_pseudo_n2_linear( REAL h_ ) : h(h_) {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		REAL r2min_1st = std::numeric_limits<REAL>::max();
		REAL r2min_2nd = std::numeric_limits<REAL>::max();
		OTYPE value_1st = 0, value_2nd = 0;
		for( size_t i = 0 ; i < data_points.size() ; i++ ) {
			REAL dr2 = normsq( focus - data_points[i].first );
			if ( dr2 < r2min_1st ) {
				r2min_2nd = r2min_1st;
				value_2nd = value_1st;
				r2min_1st = dr2;
				value_1st = data_points[i].second ;
			} else if ( dr2 < r2min_2nd ) {
				r2min_2nd = dr2;
				value_2nd = data_points[i].second ;
			}
		}
		REAL r1 = std::sqrt( r2min_1st );
		REAL r2 = std::sqrt( r2min_2nd );
		return ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );
	}

	inline geometry::any_shape<CONFIG> support( point_type focus, REAL domain_mag ) const {
		return geometry::sphere<CONFIG>( focus, h );
	}
protected:
	REAL h;
};

}

#endif /* MUI_SAMPLER_NN_H_ */
