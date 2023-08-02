/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis  *
*                                                                            *
* This software is jointly licensed under the Apache License, Version 2.0    *
* and the GNU General Public License version 3, you may use it according     *
* to either.                                                                 *
*                                                                            *
* ** Apache License, version 2.0 **                                          *
*                                                                            *
* Licensed under the Apache License, Version 2.0 (the "License");            *
* you may not use this file except in compliance with the License.           *
* You may obtain a copy of the License at                                    *
*                                                                            *
* http://www.apache.org/licenses/LICENSE-2.0                                 *
*                                                                            *
* Unless required by applicable law or agreed to in writing, software        *
* distributed under the License is distributed on an "AS IS" BASIS,          *
* WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
* See the License for the specific language governing permissions and        *
* limitations under the License.                                             *
*                                                                            *
* ** GNU General Public License, version 3 **                                *
*                                                                            *
* This program is free software: you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation, either version 3 of the License, or          *
* (at your option) any later version.                                        *
*                                                                            *
* This program is distributed in the hope that it will be useful,            *
* but WITHOUT ANY WARRANTY; without even the implied warranty of             *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              *
* GNU General Public License for more details.                               *
*                                                                            *
* You should have received a copy of the GNU General Public License          *
* along with this program.  If not, see <http://www.gnu.org/licenses/>.      *
*****************************************************************************/

/**
 * @file util.h
 * @author Y. H. Tang
 * @date 9 December 2013
 * @brief Provides a number of utility functions used through the rest of the
 * library.
 *
 */

#ifndef UTIL_H_
#define UTIL_H_

// C headers
#include <cassert>
#include <cctype>
#include <climits>
#include <cmath>
#include <csignal>
#include <cstring>
#include <cstdio>
#include <cstdint>
#include <chrono>

// C++ headers
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <functional>
#include <stdexcept>
#include <atomic>

// STL headers
#include <algorithm>
#include <deque>
#include <list>
#include <map>
#include <set>
#include <unordered_map>
#include <vector>
#include <complex>
#include <numeric>
#include <array>

#ifdef __OMP
#include <omp.h>
#endif

#include "../storage/stream.h"
#include "../geometry/point.h"

namespace mui {

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned long long llong;

const static double PI = 3.1415926535897932385;

static bool _quiet = false;

inline void set_quiet(bool q) {
    _quiet = q;
}

template<typename REAL> inline REAL clamp( REAL x, REAL l, REAL r )
{
	return (x < l) ? l : ( (x > r) ? r : x );
}

template<typename REAL> inline REAL sgn( REAL x )
{
	return (x>0)-(x<0);
}

template<uint N> inline double powr( const double x ) {
	return powr<2>( powr<N/2>( x ) ) * powr<N%2>(x);
}

template<> inline double powr<2>( const double x ) {
	return x * x;
}

template<> inline double powr<1>( const double x ) {
	return x;
}

template<class T> inline bool almost_equal(T x, T y) {
	return (x == y) ||
		   (std::fabs(x-y) < std::numeric_limits<T>::epsilon() * std::fabs(x+y)) ||
		   (std::fabs(x-y) < std::numeric_limits<T>::min());
}

template<typename T> inline T frexp10(T arg, int &exp) {
	if(almost_equal(arg, static_cast<T>(0))) exp = 0;
	else exp = 1 + static_cast<int>(std::floor(std::log10(std::fabs(arg))));
	return arg * std::pow(10, -(exp));
}

template<typename T> inline T frexp10(T arg, long &exp) {
	if(almost_equal(arg, static_cast<T>(0))) exp = 0;
	else exp = 1 + static_cast<long>(std::floor(std::log10(std::fabs(arg))));
	return arg * std::pow(10, -(exp));
}

template<class T> inline T threshold(T x) {
	return std::numeric_limits<T>::epsilon() * std::fabs(x);
}

#ifdef __GNUC__
template<> inline double powr<0>( __attribute__((unused)) const double x ) {
#else
template<> inline double powr<0>( const double x ) {
#endif
	return 1.0;
}

template<typename T, uint D>
ostream& operator<<( ostream& stream, const point<T,D>& p )
{
	for( uint i=0; i<D; ++i ) stream << p[i];
	return stream;
}

template<typename T, uint D>
istream& operator>>( istream& stream, point<T,D>& p )
{
	for( uint i=0; i<D; ++i ) stream >> p[i];
	return stream;
}

template<typename T1, typename T2, typename T3>
struct triple {
	T1 i;
	T2 j;
	T3 k;
	triple() {}
	triple( T1 _i, T2 _j, T3 _k ) : i(_i), j(_j), k(_k) {}
	triple( const triple &another ) : i(another.i), j(another.j), k(another.k) {}

	inline triple& operator = ( const triple &another ) {
		i = another.i;
		j = another.j;
		k = another.k;
		return *this;
	}

	inline bool operator != ( const triple &another ) const {
		return i != another.i || j != another.j || k != another.k;
	}
	inline bool operator == ( const triple &another ) const {
		return !(operator ==(another) );
	}
	inline bool operator < ( const triple &another ) const {
		return i < another.i || ( i == another.i && ( j < another.j || ( j == another.j && k < another.k ) ) );
	}
};

}

#endif /* UTIL_H_ */
