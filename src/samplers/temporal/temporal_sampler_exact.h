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
 * @file temporal_sampler_exact.h
 * @author Y. H. Tang
 * @date 19 February 2014
 * @brief Temporal sampler that samples at exactly the time specified and performs no
 * interpolation.
 */

#ifndef MUI_TEMPORAL_SAMPLER_EXACT_H_
#define MUI_TEMPORAL_SAMPLER_EXACT_H_

#include <limits>
#include "../../general/util.h"
#include "../../config.h"

namespace mui {

template<typename CONFIG=default_config> class temporal_sampler_exact {
public:
	using REAL       	= typename CONFIG::REAL;
	using INT        	= typename CONFIG::INT;
	using time_type  	= typename CONFIG::time_type;
	using iterator_type = typename CONFIG::iterator_type;

	temporal_sampler_exact( time_type tol = time_type(std::numeric_limits<time_type>::epsilon()) ) {
	  int exponent;
      frexp10<time_type>( std::numeric_limits<time_type>::max(), exponent );
      time_type real_precision_time = static_cast<time_type>( exponent );
      iterator_type real_precision_it = static_cast<iterator_type>( exponent );
      tolerance_time = tol*real_precision_time;
      tolerance_it = tol*real_precision_it;
	}

	//- Filter based on time input
	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
	    for( auto i: points ) {
			time_type dt = std::abs(i.first.first - focus);
	    	if ( dt <= tolerance_time ) {
				return i.second;
			}
		}

		return TYPE();
	}

	//- Filter based on time and iterator input - both used
	template<typename TYPE>
	TYPE filter( std::pair<time_type,iterator_type> focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
		for( auto i: points ) {
			time_type dt = std::abs(i.first.first - focus.first);
			iterator_type di = std::abs(i.first.second - focus.second);
			if ( dt <= tolerance_time && di <= tolerance_it ) {
				return i.second;
			}
		}

		return TYPE();
	}

	template<typename TYPE>
	TYPE get_upper_bound( TYPE focus ) const {
		return focus;
	}

	template<typename TYPE>
	TYPE get_lower_bound( TYPE focus ) const {
		return focus;
	}

private:
	time_type tolerance_time;
	time_type tolerance_it;
};

}

#endif /* MUI_TEMPORAL_SAMPLER_EXACT_H_ */
