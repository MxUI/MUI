/* -*- c++ -*-
 * stream_string.h
 *
 *  Created on: Mar 18, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_STRING_H
#define MUI_STREAM_STRING_H

#include <string>
#include "stream.h"

namespace mui {
inline istream& operator>>(istream& stream, std::string& ret)
{
	std::size_t size;
	stream >> size;
	std::string str(size, '*');
	for( char& a: str ) stream >> a;
	ret.swap(str);
	return stream;
}
inline ostream& operator<< ( ostream& stream, const std::string& str )
{
	std::size_t size = {str.size()};
	stream << size;
	stream.write(str.c_str(), str.size());
	return stream;
}
}

#endif
