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
 * @file stream_tuple.h
 * @author S. Kudo
 * @date 26 March 2014
 * @brief Defines the stream in/out for std::tuple data type.
 *
 * The order of elements in stream-in/out is first-argument to last-argument.
 */

#ifndef MUI_STREAM_TUPLE_H
#define MUI_STREAM_TUPLE_H

#include <tuple>

#include "stream.h"

namespace mui {

// build integer sequence 0, 1, ... in compile time
template<std::size_t... indexes>
struct index_sequence {
	typedef index_sequence<indexes...,sizeof...(indexes)> next;
};
template<std::size_t N>
struct make_index_sequence {
	typedef typename make_index_sequence<N-1>::type::next type;
};
template<>
struct make_index_sequence<0>{
	typedef index_sequence<> type;
};


namespace {
template<std::size_t N, std::size_t MAX, typename... Args>
struct input_tuple_impl_ {
	static void apply( istream& stream, std::tuple<Args...>& t){
		stream >> std::get<N>(t);
		input_tuple_impl_<N+1,MAX,Args...>::apply(stream,t);
	}
};
template<std::size_t N, typename... Args>
struct input_tuple_impl_<N,N,Args...> {
	static void apply( istream& , std::tuple<Args...>& ){}
};
}

template<typename... Args>
inline istream& operator>> ( istream& stream, std::tuple<Args...>& ret )
{
	std::tuple<Args...> t;
	input_tuple_impl_<0,sizeof...(Args),Args...>::apply(stream,t);
	ret.swap(t);
	return stream;
}


namespace{
template<std::size_t N, std::size_t MAX, typename... Args>
struct output_tuple_impl_ {
	static void apply( ostream& stream, const std::tuple<Args...>& t){
		stream << std::get<N>(t);
		output_tuple_impl_<N+1,MAX,Args...>::apply(stream,t);
	}
};
template<std::size_t N, typename... Args>
struct output_tuple_impl_<N,N,Args...> {
	static void apply( ostream& , const std::tuple<Args...>& ){}
};
}

template<typename... Args>
inline ostream& operator<< ( ostream& stream, const std::tuple<Args...>& t )
{
	output_tuple_impl_<0,sizeof...(Args),Args...>::apply(stream,t);
	return stream;
}

}

#endif
