/*
 * message.h
 *
 *  Created on: Mar 11, 2014
 *      Author: ytang
 */

#ifndef MESSAGE_H_
#define MESSAGE_H_

#include <memory>
#include "util.h"
#include "stream.h"
#include "stream_string.h"
#include "stream_vector.h"
#include "stream_tuple.h"

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






