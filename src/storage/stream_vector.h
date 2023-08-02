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
 * @file stream_vector.h
 * @author S. Kudo
 * @date 18 March 2014
 * @brief Defines the stream in/out for std::vector data type.
 *
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
