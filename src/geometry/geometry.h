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
 * @file geometry.h
 * @author S. Kudo
 * @date 12 March 2014
 * @brief Base classes for creating geometries, primarily used by spatial
 * interpolation methods and for defining smart MPI communication map.
 */

#ifndef MUI_GEOMETRY_H
#define MUI_GEOMETRY_H

#include <memory>

#include "../config.h"
#include "../storage/stream.h"
#include "../storage/stream_string.h"

namespace mui {
namespace geometry{

enum struct shape_type: std::int8_t {
	universe = 0,
	or_set,
	point,
	sphere,
	box
};

template<typename CONFIG> class shape;
template<typename CONFIG> class any_shape;
template<typename CONFIG> class box;
template<typename CONFIG> std::unique_ptr<shape<CONFIG> > shape_factory( shape_type );

template<typename CONFIG>
class shape {
protected: // make those ctors & copies protected to prevent slicing
	shape() = default;
	shape(shape&&) noexcept = default;
	shape(const shape&) = default;
	shape& operator=(shape&&) noexcept = default;
	shape& operator=(const shape&) = default;
public:
	virtual ~shape(){}
	virtual shape* clone() const = 0;
	virtual shape_type type() const noexcept = 0;
	virtual box<CONFIG> bbox() const = 0;

	virtual void serialize(ostream& stream) const = 0;
	virtual void deserialize(istream& stream) = 0;
};

template<typename CONFIG>
class any_shape {
public:
	any_shape() = default;
	any_shape( const any_shape& rhs )
		: content(rhs.content ? rhs.content->clone() : 0) {}
	any_shape( any_shape&& ) noexcept = default;
	any_shape( const shape<CONFIG>& rhs ) : content(rhs.clone()) {}
	any_shape& operator=( any_shape rhs ){
		rhs.swap(*this);
		return *this;
	}

	explicit operator bool() const noexcept { return static_cast<bool>(content); }
	bool empty() const noexcept { return !static_cast<bool>(content); }
	void swap( any_shape& rhs ) noexcept { content.swap(rhs.content); }

	shape<CONFIG>& get() { return *content; }
	const shape<CONFIG>& get() const { return *content; }
	
	shape_type type() const noexcept { return content ? content->type() : shape_type::universe;  }
	box<CONFIG> bbox() const;

	friend ostream& operator<<(ostream& stream, const any_shape& obj){
		stream << static_cast<std::int8_t>(obj.type());
		if( !obj.empty() ) obj.get().serialize(stream);
		return stream;
	}
	friend istream& operator>>(istream& stream, any_shape& obj){
		std::int8_t tp;
		stream >> tp;
		if( static_cast<shape_type>(tp) == shape_type::universe ) obj = any_shape();
		else {
			auto tmp = shape_factory<CONFIG>(static_cast<shape_type>(tp));
			tmp->deserialize(stream);
			obj.content.swap(tmp);
		}
		return stream;
	}
	
private:
	std::unique_ptr<shape<CONFIG> > content;
};
template<typename CONFIG> any_shape<CONFIG> get_universe_set() { return any_shape<CONFIG>(); }

template<typename CONFIG>
class point: public shape<CONFIG> {
	typedef typename CONFIG::point_type coordinate_type;
	typedef typename CONFIG::REAL REAL;
public:
	point() = default;
	point(const coordinate_type& c__) : center(c__) {}

	coordinate_type& get_center() { return center; }
	coordinate_type get_center() const { return center; }

	shape<CONFIG>* clone() const { return static_cast<shape<CONFIG>*>(new point(*this)); }
	shape_type type() const noexcept { return shape_type::point; }
	box<CONFIG> bbox() const;

	void serialize(ostream& stream) const { stream << center; }
	void deserialize(istream& stream) { stream >> center; }
private:
	coordinate_type center;
};

template<typename CONFIG>
class sphere: public shape<CONFIG> {
	typedef typename CONFIG::point_type coordinate_type;
	typedef typename CONFIG::REAL REAL;
public:
	sphere() = default;
	sphere(const coordinate_type& c__, REAL r__) : center(c__), radius(r__) {}

	coordinate_type& get_center() { return center; }
	coordinate_type get_center() const { return center; }
	REAL& get_radius() { return radius; }
	REAL get_radius() const { return radius; }

	shape<CONFIG>* clone() const { return static_cast<shape<CONFIG>*>(new sphere(*this)); }
	shape_type type() const noexcept { return shape_type::sphere; }
	box<CONFIG> bbox() const;

	void serialize(ostream& stream) const {
		stream << center;
		stream << radius;
	}
	void deserialize(istream& stream) {
		stream >> center;
		stream >> radius;
	}
private:
	coordinate_type center;
	REAL radius;
};

template<typename CONFIG>
class box: public shape<CONFIG> {
	typedef typename CONFIG::point_type coordinate_type;
	typedef typename CONFIG::REAL REAL;
public:
	box() = default;
	box(const coordinate_type& p1, const coordinate_type& p2){
		// there are four different (p1,p2) pair for same box.
		// box only allows min[i] <= max[i] for all i.
		for( uint i=0; i<CONFIG::D; ++i ){
			min[i] = std::min(p1[i],p2[i]);
			max[i] = std::max(p1[i],p2[i]);
		}
	}

	coordinate_type  get_min() const { return min; }
	coordinate_type& get_min() { return min; }
	coordinate_type  get_max() const { return max; }
	coordinate_type& get_max() { return max; }
	
	shape<CONFIG>* clone() const { return static_cast<shape<CONFIG>*>(new box(*this)); }
	shape_type type() const noexcept { return shape_type::box; }
	box<CONFIG> bbox() const;

	void serialize(ostream& stream) const {
		stream << min;
		stream << max;
	}
	void deserialize(istream& stream) {
		stream >> min;
		stream >> max;
	}
private:
	coordinate_type min, max; // only use the two vertex of box such that min[i] < max[i] for all i.
};

template<typename CONFIG>
class or_set: public shape<CONFIG> {
	typedef typename CONFIG::point_type coordinate_type;
	typedef typename CONFIG::REAL REAL;
public:
	or_set() = default;
	or_set(any_shape<CONFIG> obj1, any_shape<CONFIG>& obj2): lhs_(std::move(obj1)), rhs_(std::move(obj2)) {}

	const any_shape<CONFIG>& left() const { return lhs_; }
	any_shape<CONFIG>& left() { return lhs_; }
	const any_shape<CONFIG>& right() const { return rhs_; }
	any_shape<CONFIG>& right() { return rhs_; }
	
	shape<CONFIG>* clone() const { return static_cast<shape<CONFIG>*>(new or_set(*this)); }
	shape_type type() const noexcept { return shape_type::or_set; }
	box<CONFIG> bbox() const;

	void serialize(ostream& stream) const {
		stream << lhs_;
		stream << rhs_;
	}
	void deserialize(istream& stream) {
		stream >> lhs_;
		stream >> rhs_;
	}
private:
	any_shape<CONFIG> lhs_, rhs_;
};

template<typename CONFIG> std::unique_ptr<shape<CONFIG> > shape_factory(shape_type type)
{
	shape<CONFIG>* ptr;
	switch(type){
	case shape_type::or_set: ptr = new or_set<CONFIG>; break;
	case shape_type::point: ptr = new point<CONFIG>; break;
	case shape_type::sphere: ptr = new sphere<CONFIG>; break;
	case shape_type::box: ptr = new box<CONFIG>; break;
	default: ptr = 0; break;
	}
	return std::unique_ptr<shape<CONFIG> >(ptr);
}

namespace {
#define CASE(LHS, RHS) case shape_type::RHS : \
	return collide(static_cast<const LHS<CONFIG>&>(lhs), static_cast<const RHS<CONFIG>&>(rhs))
template<typename CONFIG> bool collide_impl_(const shape<CONFIG>& lhs, const shape<CONFIG>& rhs)
{
	switch(lhs.type()){
	case shape_type::universe:
		return true;
	case shape_type::or_set : {
		auto& obj = static_cast<const or_set<CONFIG>&>(lhs);
		return collide(obj.left(),rhs) || collide(obj.right(),rhs);
	}
	case shape_type::point :
		switch(rhs.type()){
		CASE(point,point);
		CASE(point,sphere);
		CASE(point,box);
		default: return true;
		}
	case shape_type::sphere :
		switch(rhs.type()){
		CASE(sphere,sphere);
		CASE(sphere,box);
		default: return true;
		}
	case shape_type::box :
		switch(rhs.type()){
		CASE(box,box);
		default: return true;
		}
	}
	return true; // cannot reach here
}
#undef CASE
}

template<typename CONFIG> bool collide( const shape<CONFIG>& lhs, const shape<CONFIG>& rhs)
{
	return static_cast<std::int8_t>(lhs.type()) <= static_cast<std::int8_t>(rhs.type()) ?
		collide_impl_(lhs,rhs):
		collide_impl_(rhs,lhs);
}

template<typename CONFIG> bool collide( const any_shape<CONFIG>& lhs, const any_shape<CONFIG>& rhs)
{
	if( lhs.empty() || rhs.empty() ) return true;
	else return collide(lhs.get(), rhs.get());
}

template<typename CONFIG> bool collide( const any_shape<CONFIG>& lhs, const shape<CONFIG>& rhs)
{
	if( lhs.empty() ) return true;
	else return collide(lhs.get(), rhs);
}

template<typename CONFIG> bool collide( const shape<CONFIG>& lhs, const any_shape<CONFIG>& rhs)
{
	if( rhs.empty() ) return true;
	else return collide(lhs, rhs.get());
}

template<typename CONFIG> bool collide( const point<CONFIG>& lhs, const point<CONFIG>& rhs)
{
	for( int i=0; i<CONFIG::D; ++i ) if( lhs.get_center()[i] != rhs.get_center()[i] ) return false;
	return true;
}

template<typename CONFIG> bool collide( const point<CONFIG>& lhs, const sphere<CONFIG>& rhs)
{
	return normsq( lhs.get_center()-rhs.get_center() ) <= (rhs.get_radius()*rhs.get_radius());
}

template<typename CONFIG> bool collide( const point<CONFIG>& lhs, const box<CONFIG>& rhs)
{
	for( uint i=0; i<CONFIG::D; ++i )
		if( lhs.get_center()[i] < rhs.get_min()[i] || rhs.get_max()[i] < lhs.get_center()[i] ) return false;
	return true;
}

template<typename CONFIG> bool collide( const sphere<CONFIG>& lhs, const sphere<CONFIG>& rhs)
{
	return normsq( lhs.get_center() - rhs.get_center() ) <=
		(lhs.get_radius()+rhs.get_radius())*(lhs.get_radius()+rhs.get_radius());
}

namespace {
template<typename REAL>
REAL minimum_of_quadratic(REAL b, REAL c, REAL x0, REAL x1 )
{
	// return the minimum of x^2+bx+c in [x0,x1]
	REAL p0 = -0.5*b;
	return std::max<REAL>( ( x0 <= p0 && p0 <= x1 )?
		c+0.5*b*p0 :
		std::min(x0*x0+b*x0, x1*x1+b*x1) + c, 0.0 );
}
}
template<typename CONFIG> bool collide( const sphere<CONFIG>& lhs, const box<CONFIG>& rhs)
{
	auto a = lhs.get_center()[0];
	auto dist = minimum_of_quadratic(-(a+a),a*a,rhs.get_min()[0],rhs.get_max()[0]);
	for( uint i=1; i<CONFIG::D; ++i ){
		auto a = lhs.get_center()[i];
		dist += minimum_of_quadratic(-(a+a),a*a,rhs.get_min()[i],rhs.get_max()[i]);
	}
	return dist <= lhs.get_radius()*lhs.get_radius();
}

template<typename CONFIG> bool collide( const box<CONFIG>& lhs, const box<CONFIG>& rhs)
{
    for( uint i=0; i<CONFIG::D; ++i ) {
    	if( lhs.get_max()[i] < rhs.get_min()[i] || rhs.get_max()[i] < lhs.get_min()[i] ) return false;
    }

	return true;
}

template<typename CONFIG> box<CONFIG> any_shape<CONFIG>::bbox() const
{
	if( content ) return content->bbox();
	else {
		box<CONFIG> bx;
		for( uint i = 0; i<CONFIG::D; ++i ) {
			bx.get_min()[i] = -std::numeric_limits<typename CONFIG::REAL>::infinity();
			bx.get_max()[i] = std::numeric_limits<typename CONFIG::REAL>::infinity();
		}
		return bx;
	}
}

template<typename CONFIG> box<CONFIG> point<CONFIG>::bbox() const
{
	return box<CONFIG>(center, center);
}

template<typename CONFIG> box<CONFIG> sphere<CONFIG>::bbox() const
{
	box<CONFIG> bx;
	for( uint i=0; i<CONFIG::D; ++i ){
		bx.get_min()[i] = center[i]-radius;
		bx.get_max()[i] = center[i]+radius;
	}
	return bx;
}

template<typename CONFIG>
box<CONFIG> box<CONFIG>::bbox() const
{
	return *this;
}

template<typename CONFIG> box<CONFIG> or_set<CONFIG>::bbox() const
{
	box<CONFIG> bx, lhs = lhs_.bbox(), rhs = rhs_.bbox();
	for( uint i=0; i<CONFIG::D; ++i ){
		bx.get_min()[i] = std::min(lhs.get_min()[i], rhs.get_min()[i]);
		bx.get_max()[i] = std::min(lhs.get_max()[i], rhs.get_max()[i]);
	}
	return bx;
}

template<typename CONFIG> ostream& operator<<( ostream& stream, const shape<CONFIG>& obj )
{
	obj.serialize(stream);
	return stream;
}

template<typename CONFIG> istream& operator>>( istream& stream, shape<CONFIG>& obj )
{
	obj.deserialize(stream);
	return stream;
}

template<typename CONFIG> or_set<CONFIG> operator|(any_shape<CONFIG> lhs, any_shape<CONFIG> rhs)
{
	return or_set<CONFIG>(std::move(lhs), std::move(rhs));
}

} // geometry
} // mui

#endif
