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
 * @file reader_variable.h
 * @author Y. H. Tang
 * @date 11 March 2014
 * @brief Creates a structure to parse a message as variables and pass them
 * to a function as arguments.
 */

#ifndef READER_VARIABLE_H_
#define READER_VARIABLE_H_

#include <tuple>
#include <functional>
#include <utility>
#include "message.h"
#include "../../storage/stream.h"
#include "../../storage/stream_tuple.h"
#include "../../storage/stream_string.h"

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
