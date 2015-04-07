/* $Id: chrono_sampler_null.h 1458 2014-11-05 14:21:38Z ytang $
 * Dummy sampler, intended as a file template for creating new samplers
 */

#ifndef MUI_SAMPLER_TIME_NULL_H_
#define MUI_SAMPLER_TIME_NULL_H_

#include "util.h"
#include "config.h"

namespace mui {

template<typename CONFIG=default_config> class chrono_sampler_null {
public:
	using REAL       = typename CONFIG::REAL;
	using INT        = typename CONFIG::INT;
	using time_type  = typename CONFIG::time_type;
	
	chrono_sampler_null() {
		// to do: initialization
	}

	template<typename TYPE>
	TYPE filter( time_type focus, const std::vector<std::pair<time_type, TYPE> > &points ) const {
		// to do: interpolation algorithm
	}
	time_type get_upper_bound( time_type focus ) const {
		// to do: return newest time needed with regard to focus
	}
	time_type get_lower_bound( time_type focus ) const {
		// to do: return oldest time needed with regard to focus
	}
};

}

#endif /* MUI_SAMPLER_TIME_NULL_H_ */
