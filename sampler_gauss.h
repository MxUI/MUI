/*
 * sampler_gauss.h
 *
 *  Created on: Feb 10, 2014
 *      Author: ytang
 *
 *  Gaussian sampler
 */

#ifndef MUI_SAMPLER_GAUSS_H_
#define MUI_SAMPLER_GAUSS_H_

#include "util.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_gauss {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_gauss( REAL r_, REAL h_ ) : r(r_), h(h_), nh(std::pow(2*PI*h,-0.5*CONFIG::D)) {}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		REAL  wsum = 0;
		OTYPE vsum = 0;
		for(INT i = 0 ; i < data_points.size() ; i++) {
			auto d = normsq( focus - data_points[i].first );
			if ( d < r*r ) {
				REAL w = nh * std::exp( (-0.5/h) * d );
				vsum += data_points[i].second * w;
				wsum += w;
			}
		}
		if ( wsum ) return vsum / wsum;
		else return 0.;
	}

	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		return geometry::sphere<CONFIG>( focus, r );
	}

protected:
	REAL r;
	REAL h;
	REAL nh;
};

}

#endif /* MUI_SAMPLER_GAUSS_H_ */
