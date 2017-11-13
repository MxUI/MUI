/*
Multiscale Universal Interface Code Coupling Library

Copyright (C) 2017 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis

This software is jointly licensed under the Apache License, Version 2.0
and the GNU General Public License version 3, you may use it according
to either.

** Apache License, version 2.0 **

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

** GNU General Public License, version 3 **

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

** File Details **

Filename: mui.h
Created: Feb 19, 2014
Author: Y. H. Tang
Description: Gateway header for MUI
*/

#ifndef MUI_H_
#define MUI_H_

#include "chrono_sampler_exact.h"
#include "chrono_sampler_gauss.h"
#include "chrono_sampler_sum.h"
#include "chrono_sampler_mean.h"
#include "comm.h"
#include "comm_mpi.h"
#include "comm_mpi_smart.h"
#include "dim.h"
#include "lib_mpi_split.h"
#include "lib_mpi_multidomain.h"
#include "sampler_exact.h"
#include "sampler_gauss.h"
#include "sampler.h"
#include "sampler_mov_avg.h"
#include "sampler_nn.h"
#include "sampler_pseudo_nn.h"
#include "sampler_pseudo_n2_linear.h"
#include "sampler_sum_quintic.h"
#include "sampler_sph_quintic.h"
#include "sampler_shepard_quintic.h"
#include "uniface.h"
#include "util.h"

namespace mui {

#define DECLARE_SAMPLER_1ARG(SAMPLER,SUFFIX,CONFIG)	\
	template<typename T> using SAMPLER ## SUFFIX = SAMPLER<T,T,CONFIG>;
#define DECLARE_SAMPLER_0ARG(SAMPLER,SUFFIX,CONFIG)	\
	using SAMPLER ## SUFFIX = SAMPLER<CONFIG>;

#define SPECIALIZE(SUFFIX,REALTYPE,INTTYPE,DIM) \
		typedef struct config_##SUFFIX {\
			static const int D = DIM;\
			using REAL = REALTYPE;\
			using INT  = INTTYPE;\
			using point_type = point<REAL,D>;\
			using time_type  = REAL;\
			static const bool DEBUG = false;\
			using data_types = type_list<int,double,float>;\
			using EXCEPTION = exception_segv;\
		} mui_config_##SUFFIX;\
		using uniface##SUFFIX = uniface<config_##SUFFIX>;\
		using point##SUFFIX = point<config_##SUFFIX::REAL,config_##SUFFIX::D>;\
		DECLARE_SAMPLER_1ARG(sampler_nearest_neighbor,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest_neighbor,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest2_linear,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_moving_average,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_exact,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_gauss,SUFFIX,config_##SUFFIX)\
		DECLARE_SAMPLER_0ARG(chrono_sampler_exact,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_gauss,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_sum,SUFFIX,config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_mean,SUFFIX,config_##SUFFIX);\
		namespace geometry {\
			using point##SUFFIX = point<config_##SUFFIX>;\
			using sphere##SUFFIX = sphere<config_##SUFFIX>;\
			using box##SUFFIX = box<config_##SUFFIX>;\
		}


SPECIALIZE(1d,double,int,1);
SPECIALIZE(2d,double,int,2);
SPECIALIZE(3d,double,int,3);
SPECIALIZE(1dx,double,int64_t,1);
SPECIALIZE(2dx,double,int64_t,2);
SPECIALIZE(3dx,double,int64_t,3);
SPECIALIZE(1f,float,int,1);
SPECIALIZE(2f,float,int,2);
SPECIALIZE(3f,float,int,3);
SPECIALIZE(1fx,float,int64_t,1);
SPECIALIZE(2fx,float,int64_t,2);
SPECIALIZE(3fx,float,int64_t,3);

#undef SPECIALIZE

// sample configuration class
// modify to adapt to specific simulation scenario
/******************************************************************************
struct config_custom {
	static const int D = 3;
	using REAL = double;
	using INT  = int;
	using point_type = point<REAL,D>;
	using time_type  = REAL;
	static const bool DEBUG = false;
	using data_types = type_list<int,double,float>;
	using EXCEPTION = exception_segv;
};
******************************************************************************/

// usage: SPECIALIZE( _your_custom_suffix, your_custom_config_structure )
//        uniface_your_custom_suffix interface; ...
#define SPECIALIZE(SUFFIX,CONFIG) \
		namespace mui {\
		using uniface##SUFFIX = uniface<CONFIG>;\
		using point##SUFFIX = point<CONFIG::REAL,CONFIG::D>;\
		DECLARE_SAMPLER_1ARG(sampler_nearest_neighbor,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest_neighbor,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest2_linear,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_moving_average,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_exact,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_1ARG(sampler_gauss,SUFFIX,CONFIG)\
		DECLARE_SAMPLER_0ARG(chrono_sampler_exact,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_gauss,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_sum,SUFFIX,CONFIG);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_mean,SUFFIX,CONFIG);\
		}

}

#endif /* MUI_H_ */
