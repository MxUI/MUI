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
 * @file temporal_sampler_sum.h
 * @author Y. H. Tang
 * @date 15 April 2014
 * @brief Temporal sampler that sums in time ranging from
 * [ now - left, now + right ].
 */

#ifndef MUI_TEMPORAL_SAMPLER_SUM_H_
#define MUI_TEMPORAL_SAMPLER_SUM_H_

#include "../../general/util.h"
#include "../../config.h"

namespace mui {

template<typename CONFIG=default_config> class temporal_sampler_sum {
public:
	using REAL       	= typename CONFIG::REAL;
	using INT        	= typename CONFIG::INT;
	using time_type  	= typename CONFIG::time_type;
	using iterator_type = typename CONFIG::iterator_type;
	
	temporal_sampler_sum( time_type left = time_type(0), time_type right = time_type(0) ) {
		left_ = left;
		right_ = right;
	}

	//- Filter based on time input
	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
		TYPE sum = TYPE(0);

		for( auto i: points ) {
			if ( i.first.first <= focus + right_ && i.first.first >= focus - left_ ) {
				sum += i.second;
			}
		}

		return sum;
	}

	//- Filter based on time and iterator input - only time used
	template<typename TYPE>
	TYPE filter( std::pair<time_type,iterator_type> focus, const std::vector<std::pair<std::pair<time_type,iterator_type>, TYPE> > &points ) const {
		TYPE sum = TYPE(0);

		for( auto i: points ) {
			if ( i.first.first <= focus.first + right_ && i.first.first >= focus.first - left_ ) {
				sum += i.second;
			}
		}

		return sum;
	}

	time_type get_upper_bound( time_type focus ) const {
		return focus + right_;
	}

	time_type get_lower_bound( time_type focus ) const {
		return focus - left_;
	}

	time_type tolerance() const {
		return time_type(0);
	}	

protected:
	time_type left_;
	time_type right_;
};

}

#endif /* MUI_TEMPORAL_SAMPLER_SUM_H_ */
