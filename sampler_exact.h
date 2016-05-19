/*
 * sampler_exact.h
 *
 *  Created on: Feb 10, 2014
 *      Author: ytang
 *
 *  exact point sampler
 */

#ifndef MUI_SAMPLER_EXACT_H_
#define MUI_SAMPLER_EXACT_H_

#include <limits>
#include "config.h"
#include "sampler.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_exact {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_exact() {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		for(INT i = 0 ; i < data_points.size() ; i++) {
			if ( norm( focus - data_points[i].first ) < std::numeric_limits<REAL>::epsilon() ) {
				return data_points[i].second;
			}
		}
		std::cerr << "sampler exact: hit nothing\n";
		return OTYPE(0.);
	}
	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		return geometry::sphere<CONFIG>( focus, 1 );
	}
};

}

#endif /* MUI_SAMPLER_EXACT_H_ */
