/*
 * sampler_nn.h
 *
 *  Created on: Feb 10, 2014
 *      Author: ytang
 *
 *  nearest neighbor sampler
 */

#ifndef MUI_SAMPLER_NN_H_
#define MUI_SAMPLER_NN_H_

#include "config.h"
#include "sampler.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_nearest_neighbor {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_nearest_neighbor() {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		REAL r2min = std::numeric_limits<REAL>::max();
		OTYPE value = 0;
		for(INT i = 0 ; i < data_points.size() ; i++) {
			REAL dr2 = normsq( focus - data_points[i].first );
			if ( dr2 < r2min ) {
				r2min = dr2;
				value = data_points[i].second ;
			}
		}
		return value;
	}
	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		return geometry::point<CONFIG>( focus );
	}
};

}

#endif /* MUI_SAMPLER_NN_H_ */
