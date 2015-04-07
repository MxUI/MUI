/*
 * chrono_sampler_exact.h
 *
 *  Created on: Feb 19, 2014
 *      Author: ytang
 *
 *  sampling exact the time frame as specified
 */

#ifndef MUI_CHRONO_SAMPLER_EXACT_H_
#define MUI_CHRONO_SAMPLER_EXACT_H_

#include "util.h"
#include "config.h"

namespace mui {

template<typename CONFIG=default_config> class chrono_sampler_exact {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	
	chrono_sampler_exact( time_type tol = time_type(0) ) {
		tolerance = tol;
	}

	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<time_type, TYPE> > &points ) const {
		for( auto i: points ) {
			if ( std::abs(i.first - focus) <= tolerance ) {
				return i.second;
			}
		}
		return TYPE(0);
	}
	time_type get_upper_bound( time_type focus ) const {
		return focus + tolerance;
	}
	time_type get_lower_bound( time_type focus ) const {
		return focus - tolerance;
	}

protected:
	time_type tolerance;
};

}

#endif /* MUI_CHRONO_SAMPLER_EXACT_H_ */
