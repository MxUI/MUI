/* -*- c++ -*-
 * config.h
 *
 *  Created on: Mar 10, 2013
 *      Author: skudo
 */

#ifndef CONFIG_H
#define CONFIG_H

#include "util.h"
#include "dim.h"
#include "exception.h"

namespace mui {

template<typename... TYPES> struct type_list {};

struct crunch {
	static const int D = 3;

	using REAL = double;
	using INT  = int64_t;
	using point_type = point<REAL,D>;
	using time_type  = REAL; // INT-typed time stamp might be an alternative
	using data_types = type_list<int,double,float>;

	static const bool DEBUG = false;
	using EXCEPTION = exception_segv;
};

// backward-compatibility
struct default_config : crunch {};

/*
 * user can define his own config like this:
 * struct my_config: mui::crunch {
 *   typedef unsigned int INT;
 *
 *   typedef std::array<REAL,D> point_type;
 *  };
 */
}

#endif
