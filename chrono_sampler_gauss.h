/*
 * sampler_time_gauss.h
 *
 *  Created on: Feb 14, 2014
 *      Author: ytang
 *
 *  Gaussian sampler in time, symmetric on past & future
 */

#ifndef MUI_SAMPLER_TIME_GAUSS_H_
#define MUI_SAMPLER_TIME_GAUSS_H_

#include "util.h"
#include "config.h"

namespace mui {

template<typename CONFIG=default_config> class chrono_sampler_gauss {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	
	chrono_sampler_gauss( time_type newcutoff, REAL newsigma ) {
		sigma  = newsigma;
		cutoff = newcutoff;
	}

	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<time_type, TYPE> > &points ) const {
		REAL wsum = REAL(0);
		TYPE vsum = TYPE(0);
		for( auto i: points ) {
			time_type dt = std::abs(i.first - focus);
			if ( dt < cutoff ) {
				REAL w = pow( 2*PI*sigma, -0.5 ) * exp( -0.5 * dt * dt / sigma );
				vsum += i.second * w;
				wsum += w;
			}
		}
		return ( wsum > std::numeric_limits<REAL>::epsilon() ) ? ( vsum / wsum ) : TYPE(0);
	}
	time_type get_upper_bound( time_type focus ) const {
		return focus + cutoff;
	}
	time_type get_lower_bound( time_type focus ) const {
		return focus - cutoff;
	}

protected:
	time_type cutoff;
	REAL sigma;
};

}

#endif /* MUI_SAMPLER_GAUSS_H_ */
