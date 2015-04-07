/*
 * span.h
 *
 *  Created on: Mar 15, 2014
 *      Author: ytang
 */

#ifndef SPAN_H_
#define SPAN_H_

#include "util.h"
#include "config.h"
#include "geometry.h"
#include "stream.h"
#include "stream_vector.h"

namespace mui {

template<class CONFIG>
class span {
protected:
	using container = std::vector<geometry::any_shape<CONFIG> >;
	container parts;
public:
	span() {}

	span operator || ( const geometry::any_shape<CONFIG> &s ) {
		parts.push_back(s);
		return *this;
	}
	void operator && ( const geometry::any_shape<CONFIG> &s ) {
		auto i = parts.begin();
		while( i != parts.end() ) {
			if ( !i->collide(s) ) i = parts.erase(i);
			else ++i;
		}
	}
	void reset() { parts.clear(); }

	bool collide( const geometry::any_shape<CONFIG>& other ) const {
		for( const auto& my: parts )
			if( mui::geometry::collide(my,other) ) return true;
		return false;
	}
	bool collide( const span& other ) const {
		for( const auto& r: other.parts )
			if( this->collide(r) ) return true;
		return false;
	}

	geometry::box<CONFIG> bbox() const {
		if( parts.empty() ) return geometry::box<CONFIG>();
		geometry::box<CONFIG> bx = parts.front().bbox();
		for( std::size_t i=1; i<parts.size(); ++i ) bx = mui::geometry::bbox(bx,parts[i]);
		return bx;
	}

	void serialize( ostream& stream ) const {
		stream << parts;
	}
	void deserialize( istream& stream ){
		stream >> parts;
	}
};

template<typename CONFIG>
ostream& operator<< ( ostream& stream, const span<CONFIG>& data )
{
	data.serialize(stream);
	return stream;
}

template<typename CONFIG>
istream& operator>> ( istream& stream, span<CONFIG>& data )
{
	data.deserialize(stream);
	return stream;
}

template<typename CONFIG>
bool collide( const span<CONFIG>& lhs, const span<CONFIG>& rhs )
{
	return lhs.collide(rhs);
}

template<typename CONFIG>
bool collide( const geometry::any_shape<CONFIG>& lhs, const span<CONFIG>& rhs )
{
	return rhs.collide(lhs);
}
template<typename CONFIG>
bool collide( const span<CONFIG>& lhs, const geometry::any_shape<CONFIG>& rhs )
{
	return lhs.collide(rhs);
}

}

#endif /* SPAN_H_ */
