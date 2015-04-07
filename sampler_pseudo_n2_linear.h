/*
 * sampler_pseudo_nn.h
 *
 *  Created on: Oct 10, 2014
 *      Author: ytang
 *
 *  nearest neighbor sampler
 */

#ifndef MUI_SAMPLER_PSEUDO_N2_LINEAR_H_
#define MUI_SAMPLER_PSEUDO_N2_LINEAR_H_

#include "config.h"
#include "sampler.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_pseudo_nearest2_linear {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_pseudo_nearest2_linear( REAL h_ ) : h(h_) {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		REAL r2min_1st = std::numeric_limits<REAL>::max();
		REAL r2min_2nd = std::numeric_limits<REAL>::max();
		OTYPE value_1st = 0, value_2nd = 0;
		for(INT i = 0 ; i < data_points.size() ; i++) {
			REAL dr2 = normsq( focus - data_points[i].first );
			if ( dr2 < r2min_1st ) {
				r2min_2nd = r2min_1st;
				value_2nd = value_1st;
				r2min_1st = dr2;
				value_1st = data_points[i].second ;
			} else if ( dr2 < r2min_2nd ) {
				r2min_2nd = dr2;
				value_2nd = data_points[i].second ;
			}
		}
		REAL r1 = std::sqrt( r2min_1st );
		REAL r2 = std::sqrt( r2min_2nd );
		return ( value_1st * r2 + value_2nd * r1 ) / ( r1 + r2 );
	}
	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		return geometry::sphere<CONFIG>( focus, h );
	}
protected:
	REAL h;
};

}

#endif /* MUI_SAMPLER_NN_H_ */
