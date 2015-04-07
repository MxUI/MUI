/*
 * mui.h
 *
 *  Created on: Feb 19, 2014
 *      Author: ytang
 *
 *  gateway header for MUI
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
		struct mui_config_##SUFFIX {\
			static const int D = DIM;\
			using REAL = REALTYPE;\
			using INT  = INTTYPE;\
			using point_type = point<REAL,D>;\
			using time_type  = REAL;\
			static const bool DEBUG = false;\
			using data_types = type_list<int,double,float>;\
			using EXCEPTION = exception_segv;\
		};\
		using uniface##SUFFIX = uniface<mui_config_##SUFFIX>;\
		using point##SUFFIX = point<mui_config_##SUFFIX::REAL,mui_config_##SUFFIX::D>;\
		DECLARE_SAMPLER_1ARG(sampler_nearest_neighbor,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest_neighbor,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_pseudo_nearest2_linear,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_moving_average,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_exact,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_1ARG(sampler_gauss,SUFFIX,mui_config_##SUFFIX)\
		DECLARE_SAMPLER_0ARG(chrono_sampler_exact,SUFFIX,mui_config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_gauss,SUFFIX,mui_config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_sum,SUFFIX,mui_config_##SUFFIX);\
		DECLARE_SAMPLER_0ARG(chrono_sampler_mean,SUFFIX,mui_config_##SUFFIX);\


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

}

#endif /* MUI_H_ */
