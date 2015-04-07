/* -*- c++ -*-
 * stream_tuple.h
 *
 *  Created on: Mar 26, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_TUPLE_H
#define MUI_STREAM_TUPLE_H

#include <tuple>
#include "stream.h"

// stream-in/out for std::tuple.
// The order of elements in stream-in/out is first-argument to last-argument.
// So then you can stream-out by tuple 
//   stream << std::make_tuple(1,3.0,std::string("string"));
// and stream-in by elements
//   int i; double d; std::string s;
//   stream >> i >> d >> s;
// and vice-versa.

namespace mui {

// build interger sequence 0, 1, ... in compile time
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
