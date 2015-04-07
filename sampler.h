/*
 * sampler.h
 *
 *  Created on: Feb 10, 2014
 *      Author: ytang
 *
 *  This is just a reference file for making custom samplers
 *  The new sampler does not have to derive from this class, it
 *  just need to implement all the interfaces with the signatures
 *  specified
 *
 *  Only very limited template programming is required
 */

#ifndef MUI_SAMPLER_H_
#define MUI_SAMPLER_H_

#include "config.h"
#include "geometry.h"
#include "virtual_container.h"

namespace mui {

#if 0
class sampler {
public:
	sampler() {
		r = 1.0;
	}

	double filter( point<double,3> focus, const virtual_container<double> &data_points ) const {
		double sum = 0.;
		int n = 0;
		for(int i = 0 ; i < data_points.size() ; i++) {
			if ( ( focus - data_points[i].first ).norm() < r ) {
				sum += data_points[i].second;
				n++;
			}
		}
		if (n) return sum / n;
		else return 0;
	}
	inline span<> support() const {
		return span<>() || geometry::sphere<>( point<double,3>(0), r );
	}

protected:
	double r;
};

#endif
}

#endif /* MUI_SAMPLER_H_ */
