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
 * @file mui.h
 * @author Y. H. Tang
 * @date 19 February 2014
 * @brief The main header file for MUI. Usually the only file that needs to
 * be included in order to integrate into an application.
 */

#ifndef MUI_H_
#define MUI_H_

//Include spatial samplers
#include "samplers/spatial/sampler_exact.h"
#include "samplers/spatial/sampler_gauss.h"
#include "samplers/spatial/sampler_mov_avg.h"
#include "samplers/spatial/sampler_nn.h"
#include "samplers/spatial/sampler_pseudo_nn.h"
#include "samplers/spatial/sampler_pseudo_n2_linear.h"
#include "samplers/spatial/sampler_sum_quintic.h"
#include "samplers/spatial/sampler_sph_quintic.h"
#include "samplers/spatial/sampler_shepard_quintic.h"
#include "samplers/spatial/sampler_rbf.h"

//Include temporal samplers
#include "samplers/temporal/temporal_sampler_exact.h"
#include "samplers/temporal/temporal_sampler_gauss.h"
#include "samplers/temporal/temporal_sampler_mean.h"
#include "samplers/temporal/temporal_sampler_sum.h"

//Include coupling algorithms
#include "samplers/algorithm/algo_fixed_relaxation.h"
#include "samplers/algorithm/algo_aitken.h"

//Include other library headers
#include "samplers/sampler.h"
#include "communication/comm.h"
#include "communication/comm_mpi.h"
#include "communication/comm_mpi_smart.h"
#include "general/dim.h"
#include "communication/lib_mpi_split.h"
#include "communication/lib_mpi_multidomain.h"
#include "uniface.h"
#include "general/util.h"
#include <string>

namespace mui {

#define DECLARE_SAMPLER_1ARG(SAMPLER,SUFFIX,CONFIG)	\
	template<typename T> using SAMPLER ## SUFFIX = SAMPLER<CONFIG,T,T>;
#define DECLARE_SAMPLER_0ARG(SAMPLER,SUFFIX,CONFIG)	\
	using SAMPLER ## SUFFIX = SAMPLER<CONFIG>;

#define SPECIALIZE(SUFFIX,REALTYPE,INTTYPE,DIM) \
		typedef struct config_##SUFFIX {\
  	  	  	using EXCEPTION = exception_segv;\
  	  	  	static const bool DEBUG = false;\
			static const int D = DIM;\
			static const bool FIXEDPOINTS = false;\
			static const bool QUIET = false;\
			using REAL = REALTYPE;\
			using INT = INTTYPE;\
			using point_type = point<REAL,D>;\
			using time_type = REAL;\
			using iterator_type = INT;\
			using data_types = type_list<uint32_t,uint64_t,int32_t,int64_t,double,float,std::string>;\
		} mui_config_##SUFFIX;\
		using uniface##SUFFIX = uniface<config_##SUFFIX>;\
		using point##SUFFIX = point<config_##SUFFIX::REAL,config_##SUFFIX::D>;\
		DECLARE_SAMPLER_1ARG(sampler_sum_quintic,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_sph_quintic,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_shepard_quintic,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_nearest_neighbor,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest_neighbor,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_n2_linear,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_moving_average,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_exact,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_gauss,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_rbf,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_0ARG(temporal_sampler_exact,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_gauss,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_sum,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_mean,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(algo_fixed_relaxation,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(algo_aitken,SUFFIX,config_##SUFFIX);\
		namespace geometry {\
			using point##SUFFIX = point<config_##SUFFIX>;\
			using sphere##SUFFIX = sphere<config_##SUFFIX>;\
			using box##SUFFIX = box<config_##SUFFIX>;\
			using or_set##SUFFIX = or_set<config_##SUFFIX>;\
		}

SPECIALIZE(1d,double,int32_t,1);
SPECIALIZE(2d,double,int32_t,2);
SPECIALIZE(3d,double,int32_t,3);
SPECIALIZE(1dx,double,int64_t,1);
SPECIALIZE(2dx,double,int64_t,2);
SPECIALIZE(3dx,double,int64_t,3);
SPECIALIZE(1f,float,int32_t,1);
SPECIALIZE(2f,float,int32_t,2);
SPECIALIZE(3f,float,int32_t,3);
SPECIALIZE(1fx,float,int64_t,1);
SPECIALIZE(2fx,float,int64_t,2);
SPECIALIZE(3fx,float,int64_t,3);

#undef SPECIALIZE

// usage: SPECIALIZE( _your_custom_suffix, your_custom_config_structure )
//        uniface_your_custom_suffix interface; ...
#define SPECIALIZE(SUFFIX,CONFIG) \
		namespace mui {\
		using uniface##SUFFIX = uniface<CONFIG>;\
		using point##SUFFIX = point<CONFIG::REAL,CONFIG::D>;\
		DECLARE_SAMPLER_1ARG(sampler_nearest_neighbor,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest_neighbor,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_n2_linear,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_moving_average,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_exact,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_rbf,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_gauss,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_0ARG(temporal_sampler_exact,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_gauss,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_sum,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(temporal_sampler_mean,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(algo_fixed_relaxation,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(algo_aitken,SUFFIX,CONFIG);\
		}

}

#endif /* MUI_H_ */
