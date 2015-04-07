/* $Id: sampler_null.h 1458 2014-11-05 14:21:38Z ytang $
 * Dummy sampler, intended as a file template for creating new samplers
 */

#ifndef MUI_SAMPLER_NULL_H_
#define MUI_SAMPLER_NULL_H_

#include "util.h"

namespace mui {

template<typename O_TP, typename I_TP=O_TP, typename CONFIG=default_config>
class sampler_null {
public:
	using OTYPE      = O_TP;
	using ITYPE      = I_TP;
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using point_type = typename CONFIG::point_type;

	sampler_null() {
		// to do: initialization
	}

	template<template<typename,typename> class CONTAINER>
	inline OTYPE filter( point_type focus, const CONTAINER<ITYPE,CONFIG> &data_points ) const {
		// to do: interpolation algorithm
	}

	inline geometry::any_shape<CONFIG> support( point_type focus ) const {
		// to do: scope of points needed expressed in terms of mui::geometry::...
	}
};

}

#endif /* MUI_SAMPLER_NULL_H_ */
