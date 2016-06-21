/*
 * sampler_mov_avg.h
 *
 *  Created on: Feb 10, 2014
 *      Author: ytang
 *
 *  moving average sampler
 */

#ifndef MUI_SAMPLER_MOVING_AVG_H_
#define MUI_SAMPLER_MOVING_AVG_H_

#include "config.h"
#include "sampler.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_moving_average {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_moving_average( point_type bbox_ ) {
		bbox  = bbox_;
	}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		INT n(0);
		OTYPE vsum(0);
		for(INT i = 0 ; i < data_points.size() ; i++) {
			point_type dx = apply( data_points[i].first - focus, abs );
			bool within = true;
			for(INT i = 0 ; within && i < CONFIG::D ; i++ ) within = within && ( dx[i] < bbox[i] );
			if ( within ) {
				vsum += data_points[i].second;
				n++;
			}
		}
		if (CONFIG::DEBUG) assert( n!=0 );
		return n ? ( vsum / OTYPE(n) ): OTYPE(0.);
	}

	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		return geometry::box<CONFIG>( focus - 0.5 * bbox, focus + 0.5 * bbox );
	}

protected:
	point_type bbox;
};

}

#endif /* MUI_SAMPLER_MOVING_AVG_H_ */
