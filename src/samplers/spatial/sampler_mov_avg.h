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
 * @file sampler_mov_avg.h
 * @author Y. H. Tang
 * @date 10 February 2014
 * @brief Spatial sampler that provides a value at a point using a moving
 * average interpolation.
 */

#ifndef MUI_SAMPLER_MOVING_AVG_H_
#define MUI_SAMPLER_MOVING_AVG_H_

#include "../../config.h"
#include "../sampler.h"
#include <cmath>

namespace mui {

template<typename CONFIG=default_config, typename O_TP=typename CONFIG::REAL, typename I_TP=O_TP>
class sampler_moving_average {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_moving_average( point_type bbox_ ) {
		bbox  = bbox_;
	}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
	    size_t n(0);
		OTYPE vsum(0);
		for( size_t i = 0 ; i < data_points.size() ; i++ ) {
            point_type dx(REAL(0.0));
            for (INT j = 0 ; j < CONFIG::D ; j++) {
                dx[j] = std::fabs(data_points[i].first[j] - focus[j]);
            }
			bool within = true;
			for( INT i = 0 ; within && i < CONFIG::D ; i++ ) within = within && ( dx[i] < bbox[i] );
			if ( within ) {
				vsum += data_points[i].second;
				n++;
			}
		}
		if (CONFIG::DEBUG) assert( n!=0 );
		return n ? ( vsum / OTYPE(n) ): OTYPE(0.);
	}

	inline geometry::any_shape<CONFIG> support( point_type focus, REAL domain_mag ) const {
	    return geometry::box<CONFIG>( focus - REAL(0.5) * bbox, focus + REAL(0.5) * bbox );
	}

protected:
	point_type bbox;
};

}

#endif /* MUI_SAMPLER_MOVING_AVG_H_ */
