/*
 * reader_variable.h
 *
 *  Created on: Mar 11, 2014
 *      Author: ytang
 */

#ifndef READER_VARIABLE_H_
#define READER_VARIABLE_H_

#include <tuple>
#include <functional>
#include <utility>
#include "message.h"
#include "stream.h"
#include "stream_tuple.h"
#include "stream_string.h"

namespace mui
{

// parse message as variables and pass them to f as arguments of it
template<typename... Args>
struct reader_variables {
	typedef std::function<void(Args...)> function_type;
	typedef std::tuple<typename std::remove_reference<Args>::type...> tuple_type;

	reader_variables() = default;
	reader_variables( function_type f ) : f_(std::move(f)) {}
	void operator()( const message& msg ){
		// parse msg as tuple of variables
		auto stream = make_istream(msg.data());
		tuple_type t;
		stream >> t;
		// split tuple before applying
		apply(t, typename make_index_sequence<sizeof...(Args)>::type());
	}
private:
	template<std::size_t... indexes>
	void apply( tuple_type& t, index_sequence<indexes...> ){
		f_(std::get<indexes>(std::move(t))...);
	}
	function_type f_;
};

}


#endif /* READER_VARIABLE_H_ */
