/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2021 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    S. M. Longshaw                                          *
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
 * @file config_c_wrapper.h
 * @author S. M. Longshaw
 * @date 12 July 2021
 * @brief File containing data structures defining all data types used
 *        by an interface using the C wrapper - based on original config.h.
 */

#ifndef CONFIG_C_WRAPPER_H
#define CONFIG_C_WRAPPER_H

namespace mui {

struct mui_c_wrapper_1D {
	using EXCEPTION = exception_segv;			//- Exception handling type

	static const bool DEBUG = false;			//- Enable extra debugging output
	static const int D = 1;						//- Dimensionality of the domain
	static const bool FIXEDPOINTS = false;		//- Enable optimisations for problems with fixed point structure
	static const bool QUIET = false;			//- If the library is quiet then it will only issue critical warning messages

	using REAL = double;						//- REAL data type
	using INT = int;							//- INT data type

	using time_type = REAL;						//- INT for iteration coupling, REAL for time-based coupling
	using iterator_type = INT;					//- Typically INT for sub-iteration count
	using point_type = point<REAL, D>;			//- "point" data type and dimensionality
	using data_types = type_list<uint32_t,		//- Data types that can be used in the interface
								 uint64_t,
								 int32_t,
								 int64_t,
								 double,
								 float,
								 std::string
								>;
};

struct mui_c_wrapper_2D {
	using EXCEPTION = exception_segv;			//- Exception handling type

	static const bool DEBUG = false;			//- Enable extra debugging output
	static const int D = 2;						//- Dimensionality of the domain
	static const bool FIXEDPOINTS = false;		//- Enable optimisations for problems with fixed point structure
	static const bool QUIET = false;			//- If the library is quiet then it will only issue critical warning messages

	using REAL = double;						//- REAL data type
	using INT = int;							//- INT data type

	using time_type = REAL;						//- INT for iteration coupling, REAL for time-based coupling
	using iterator_type = INT;					//- Typically INT for sub-iteration count
	using point_type = point<REAL, D>;			//- "point" data type and dimensionality
	using data_types = type_list<uint32_t,		//- Data types that can be used in the interface
								 uint64_t,
								 int32_t,
								 int64_t,
								 double,
								 float,
								 std::string
								>;
};

struct mui_c_wrapper_3D {
	using EXCEPTION = exception_segv;			//- Exception handling type

	static const bool DEBUG = false;			//- Enable extra debugging output
	static const int D = 3;						//- Dimensionality of the domain
	static const bool FIXEDPOINTS = false;		//- Enable optimisations for problems with fixed point structure
	static const bool QUIET = false;			//- If the library is quiet then it will only issue critical warning messages

	using REAL = double;						//- REAL data type
	using INT = int;							//- INT data type

	using time_type = REAL;						//- INT for iteration coupling, REAL for time-based coupling
	using iterator_type = INT;					//- Typically INT for sub-iteration count
	using point_type = point<REAL, D>;			//- "point" data type and dimensionality
	using data_types = type_list<uint32_t,		//- Data types that can be used in the interface
								 uint64_t,
								 int32_t,
								 int64_t,
								 double,
								 float,
								 std::string
								>;
};

}

#endif
