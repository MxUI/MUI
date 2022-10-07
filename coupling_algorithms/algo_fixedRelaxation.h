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

namespace mui {

template<typename CONFIG=default_config> class algo_fixedRelaxation {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	using point_type = typename CONFIG::point_type;

	algo_fixedRelaxation( REAL undRelxFac = 1.0, 
 	  std::vector<std::pair<point_type, REAL>> ptsVluInit = std::vector<std::pair<point_type, REAL>>() ):
 	  initUndRelxFac_(undRelxFac),
 	  ptsVluPositive_(ptsVluInit) {

	    posNegCounter_ = 1;
	    ptsVluNegative_= std::vector<std::pair<point_type, REAL>>();

	}
 
	//- relaxation based on single time value
	template<typename OTYPE>
	OTYPE relaxation(point_type focus, const OTYPE &filteredValue) const {

		posNegCounter_ *= -1;

		OTYPE relaxedValue = filteredValue;

		return relaxedValue;
	}

protected:
	REAL initUndRelxFac_;

	std::vector<std::pair<point_type, REAL>> ptsVluPositive_;
	std::vector<std::pair<point_type, REAL>> ptsVluNegative_;

	mutable INT posNegCounter_;
};

}

#endif /* MUI_COUPLING_ALGORITHM_FIXED_RELAXATION_H_ */
