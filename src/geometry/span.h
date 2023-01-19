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
 * @file span.h
 * @author Y. H. Tang
 * @date 15 March 2014
 * @brief Provides functions to determine whether geometries are colliding.
 */

#ifndef SPAN_H_
#define SPAN_H_

#include "../general/util.h"
#include "../config.h"
#include "geometry.h"
#include "../storage/stream.h"
#include "../storage/stream_vector.h"

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
