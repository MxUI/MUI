/*
 * sampler_time_mean.h
 *
 *  Created on: Oct 12, 2014
 *      Author: ytang
 *
 *  average in time, range from [ now - left, now + right ]
 */

#ifndef MUI_SAMPLER_TIME_MEAN_H_
#define MUI_SAMPLER_TIME_MEAN_H_

#include "util.h"
#include "config.h"

namespace mui {

template<typename CONFIG=default_config> class chrono_sampler_mean {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	
	chrono_sampler_mean( time_type newleft = time_type(0), time_type newright = time_type(0) ) {
		left   = newleft;
		right  = newright;
	}

	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<time_type, TYPE> > &points ) const {
		TYPE sum = TYPE(0);
		for( auto i: points ) {
			if ( i.first <= focus + right && i.first >= focus - left ) {
				sum += i.second;
			}
		}
		if ( points.size() ) return sum / TYPE(points.size());
		else return TYPE(0);
	}
	time_type get_upper_bound( time_type focus ) const {
		return focus + right;
	}
	time_type get_lower_bound( time_type focus ) const {
		return focus - left;
	}

protected:
	time_type left, right;
};

}

#endif /* MUI_SAMPLER_TIME_SUM_H_ */
