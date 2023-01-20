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
 * @file message.h
 * @author Y. H. Tang
 * @date 11 March 2014
 * @brief Structure to contain and manipulate data from internal data to MPI
 * message.
 */

#ifndef MESSAGE_H_
#define MESSAGE_H_

#include <memory>
#include "../../general/util.h"
#include "../../storage/stream.h"
#include "../../storage/stream_string.h"
#include "../../storage/stream_vector.h"
#include "../../storage/stream_tuple.h"

namespace mui {

struct message
{
public:
	using id_type = std::string;
private:
        id_type     id_;
        std::size_t id_size_;
        std::vector<char> data_;
public:
	message() : id_(), id_size_(0u), data_() {}

	template<typename... types>
	static message make( const id_type& id, types&&... data ) {
		message msg;
		msg.id_ = id;
		msg.id_size_ = streamed_size(id);
		std::size_t n = msg.id_size_ + streamed_size(data...);
		msg.data_.resize(n);
		auto stream = make_ostream(msg.data_.data());
		stream << id << std::forward_as_tuple(data...);
		return msg;
	}
	static message make( std::vector<char> data ){
		message m;
		m.data_.swap(data);
		auto in = make_istream(m.data_.begin());
		in >> m.id_;
		m.id_size_ = streamed_size(m.id_);
		return m;
	}

	bool has_id() const { return !id_.empty(); }

	const id_type& id() const { return id_; }
	const char* data() const { return data_.data() + id_size_; }
	std::size_t size() const { return data_.size() - id_size_; }
	std::vector<char> detach() {
		std::vector<char> data;
		data.swap(data_);
		id_.clear();
		id_size_ = 0;
		return data;
	}
private:
	friend ostream& operator<<( ostream&, const message&);
};

// serialization & deserialization method
inline ostream& operator<< ( ostream& stream, const message& m )
{
	stream << m.data_;
	return stream;
}

inline istream& operator>> ( istream& stream, message& m )
{
	std::vector<char> v;
	stream >> v;
	m = message::make(v);
	return stream;
}

}

#endif /* MESSAGE_H_ */
