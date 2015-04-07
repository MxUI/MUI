/* -*- c++ -*-
 * util.h
 *
 *  Created on: Dec 9, 2013
 *      Author: ytang
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

#include "stream.h"

#include "point.h"

namespace mui {

typedef unsigned int uint;
typedef unsigned long long ull;
typedef unsigned long long llong;

const static double PI = 3.1415926535897931160;

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

#ifdef __GNUC__
template<> inline double powr<0>( __attribute__((unused)) const double x ) {
#else
template<> inline double powr<0>( const double x ) {
#endif
	return 1.0;
}

///****************************************************************************
//MATRIX
//*****************************************************************************/
//
//template<typename TYPE, uint M=3, uint N=3>
//struct matrix {
//	matrix() {
//		for(uint i=0; i<M; i++) for(uint j=0; j<N; j++) _e[i][j] = 0.;
//	}
//	matrix( const matrix& other ) {
//		*this = other;
//	}
//	matrix& operator = ( const matrix& other ) {
//		for(uint i=0; i<M; i++) for(uint j=0; j<N; j++) _e[i][j] = other._e[i][j];
//		return *this;
//	}
//
//	inline const TYPE& operator () (uint row, uint col) const {
//		assert(row<M&&col<N);
//		return _e[row][col];
//	}
//	inline TYPE& operator () (uint row, uint col) {
//		assert(row<M&&col<N);
//		return _e[row][col];
//	}
//	inline int m() const {
//		return M;
//	}
//	inline int n() const {
//		return N;
//	}
//protected:
//	TYPE _e[M][N];
//};
//
//template<typename TYPE, uint M=3, uint N=3>
//struct eye: public matrix<TYPE,M,N> {
//	eye() {
//		assert(M==N);
//		for(uint i=0; i<M; i++) this->_e[i][i] = 1.;
//	}
//};
//
//template<typename TYPE, uint N, uint M, uint K>
//matrix<TYPE,N,M> operator * ( const matrix<TYPE,N,K> m1, const matrix<TYPE,K,M> &m2 )
//{
//	matrix<TYPE,N,M> r;
//	for(uint i=0; i<N; i++)
//		for(uint j=0; j<M; j++)
//			for(uint k=0; k<K; k++)
//				r(i,j) += m1(i,k) * m2(k,j);
//	return r;
//}
//
///****************************************************************************
//VECTOR/POINT/COORDINATE
//*****************************************************************************/
//
//template<typename TYPE, uint D>
//struct point {
//	TYPE x[D];
//
//	point() {}
//	point( const point &other ) {
//		*this = other;
//	}
//	point(TYPE r) {
//		for(uint i = 0 ; i < D ; i++) x[i] = r;
//	}
//	point(TYPE *a) {
//		for(uint i = 0 ; i < D ; i++) x[i] = *a++;
//	}
//	/*point(array<REAL,D> &r) {
//	    for(uint i = 0 ; i < D ; i++) x[i] = r[i];
//	    }*/
//	inline uint d() const {
//		return D;
//	}
//	inline point& operator = ( const point &other ) {
//		for(uint i = 0 ; i < D ; i++) x[i] = other.x[i];
//		return *this;
//	}
//	inline point& operator = ( const matrix<TYPE,D,1>& m ) {
//		for(uint i = 0 ; i < D ; i++ ) x[i] = m(i,0);
//		return *this;
//	}
//	inline operator matrix<TYPE,D,1> () const {
//		matrix<TYPE,D,1> m;
//		for(uint i = 0 ; i < D ; i++ ) m(i,0) = x[i];
//		return m;
//	}
//
//	inline void operator += ( const point &v ) {
//		for(uint i = 0 ; i < D ; i++) x[i] += v.x[i];
//	}
//	inline void operator -= ( const point &v ) {
//		for(uint i = 0 ; i < D ; i++) x[i] -= v.x[i];
//	}
//	inline void operator *= ( const TYPE k ) {
//		for(uint i = 0 ; i < D ; i++) x[i] *= k;
//	}
//	inline void operator *= ( const point &v ) { // element-wise product
//		for(uint i = 0 ; i < D ; i++) x[i] *= v.x[i];
//	}
//	inline void operator /= ( const TYPE k ) {
//		for(uint i = 0 ; i < D ; i++) x[i] /= k;
//	}
//	inline void operator /= ( const point &v ) { // element-wise division
//		for(uint i = 0 ; i < D ; i++) x[i] /= v.x[i];
//	}
//	inline TYPE dot( const point &v ) const {
//		TYPE s = 0.;
//		for(uint i = 0 ; i < D ; i++) s += x[i] * v.x[i];
//		return s;
//	}
//	inline point cross( const point &v ) const {
//		assert(D==3);
//		point n;
//		n[0] = x[1] * v.x[2] - v.x[1] * x[2];
//		n[1] = x[2] * v.x[0] - v.x[2] * x[0];
//		n[2] = x[0] * v.x[1] - v.x[0] * x[1];
//		return n;
//	}
//
//	inline const TYPE& operator [] ( uint i ) const { // C style
//		assert(i<D);
//		return x[i];
//	}
//	inline TYPE& operator [] ( uint i ) { // C style
//		assert(i<D);
//		return x[i];
//	}
//	inline const TYPE& operator () ( uint i ) const { // MATLAB style
//		assert(i-1<D);
//		return x[i-1];
//	}
//	inline TYPE& operator () ( uint i ) { // MATLAB style
//		assert(i-1<D);
//		return x[i-1];
//	}
//
//	inline TYPE normsq() const {
//		return dot(*this);
//	}
//	inline TYPE norm() const {
//		return std::sqrt( normsq() );
//	}
//
//	template<typename FUNCTOR>
//	inline void apply( FUNCTOR func ) {
//		for(uint i = 0 ; i < D ; i++) x[i] = func( x[i] );
//	}
//
////    template<typename TYPE>
////    inline operator TYPE () {
////    	return dynamic_cast<TYPE>(*x);
////    }
//
//	friend inline point operator + ( point v1 , point v2 ) {
//		v1 += v2;
//		return v1;
//	}
//	friend inline point operator - ( point v1 , point v2 ) {
//		v1 -= v2;
//		return v1;
//	}
//	friend inline point operator - ( point v1 ) {
//		point s(0) ;
//		s -= v1;
//		return s;
//	}
//	friend inline point operator * ( TYPE k, point v ) {
//		v *= k;
//		return v;
//	}
//	friend inline point operator * ( point v, TYPE k ) {
//		v *= k;
//		return v;
//	}
//	friend inline point operator * ( point v1, point v2 ) { // element-wise product
//		v1 *= v2;
//		return v1;
//	}
//	friend inline TYPE dot( point v1, point v2 ) { // dot product
//		return v1.dot(v2);
//	}
//	friend inline point operator / ( point v, TYPE k ) {
//		v *= 1/k;
//		return v;
//	}
//	friend inline point operator / ( TYPE k, point v ) {
//		point u(k);
//		u /= v;
//		return u;
//	}
//	friend inline point operator / ( point v1, point v2 ) {
//		v1 /= v2;
//		return v1;
//	}
//	friend inline point cross( point v1, point v2) {
//		return v1.cross(v2);
//	}
//	friend inline point floor( point v ) {
//		for(uint i = 0 ; i < D ; i++) v.x[i] = std::floor( v.x[i] );
//		return v;
//	}
//	friend inline point abs( point v ) {
//		for(uint i = 0 ; i < D ; i++) v.x[i] = std::abs( v.x[i] );
//		return v;
//	}
//	friend inline point clamp( point v, point l, point r ) {
//		for(uint i = 0 ; i < D ; i++) v.x[i] = clamp( v.x[i], l[i], r[i] );
//		return v;
//	}
//	friend inline std::ostream& operator << ( std::ostream& out, point v ) {
//		out<<'(';
//		for(int i = 0 ; i < D ; i++) out<<v[i]<<(i!=D-1?",":"");
//		out<<')';
//		return out;
//	}
//};

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
