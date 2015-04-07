/* -*- c++ -*-
 * stream_vector.h
 *
 *  Created on: Mar 18, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_VECTOR_H
#define MUI_STREAM_VECTOR_H

#include <vector>
#include "stream.h"

namespace mui {
template<typename TYPE>
inline istream& operator>>(istream& stream, std::vector<TYPE>& ret)
{
	std::size_t size;
	stream >> size;
	std::vector<TYPE> vec(size);
	for( auto& a: vec ) stream >> a;
	ret.swap(vec);
	return stream;
}
// specialization is only for char because of endian problem
inline istream& operator>>(istream& stream, std::vector<char>& ret)
{
	std::size_t size;
	stream >> size;
	std::vector<char> vec(size);
	stream.read(vec.data(), size);
	ret.swap(vec);
	return stream;
}

template<typename TYPE>
inline ostream& operator<<(ostream& stream, const std::vector<TYPE>& vec)
{
	stream << vec.size();
	for( const auto& a: vec ) stream << a;
	return stream;
}
inline ostream& operator<<(ostream& stream, const std::vector<char>& vec)
{
	stream << vec.size();
	stream.write(vec.data(), vec.size());
	return stream;
}
}

#endif
