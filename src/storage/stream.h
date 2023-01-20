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
 * @file stream.h
 * @author S. Kudo
 * @date 18 March 2014
 * @brief Defines base stream class container_stream and associated
 * functors.
 *
 * container_stream is akin to std::stringstream. Its functionality is to
 * serialize & deserialize data to container<char>.
 */

#ifndef MUI_STREAM_H_
#define MUI_STREAM_H_

#include <memory>
#include <algorithm>

#include "../general/endian_traits.h"

namespace mui {

class istream {
public:
	virtual ~istream() {}
	virtual void read( char* ptr, std::size_t size ) = 0;
};

class ostream {
public:
	virtual ~ostream() {}
	virtual void write( const char* ptr, std::size_t size ) = 0;
};

class iostream : public istream, public ostream {
public:
	virtual ~iostream() {}
};

template<template<typename T, typename =std::allocator<T> > class Seq,
         typename Alloc=std::allocator<char> >
class container_stream: public iostream {
public:
	~container_stream() {}
	void read( char* ptr, std::size_t size ){
		if( content.size() < size ) throw("");
		std::copy_n(content.begin(),size,ptr);
		content.erase(content.begin(),content.begin()+size);
	}
	void write( const char* ptr, std::size_t size ){
		content.insert(content.end(),ptr,ptr+size);
	}
	Seq<char,Alloc>& data() { return content; }
	const Seq<char,Alloc>& data() const { return content; }
private:
	Seq<char,Alloc> content;
};

/*
 * inputiterator_stream is a istream from InputIterator
 */
template<typename ConstInputIterator>
class iitr_stream: public istream {
public:
	iitr_stream( const iitr_stream& ) = default;
	iitr_stream( ConstInputIterator cur_ )
		: cur(cur_) {}
	~iitr_stream() {}
	void read( char* ptr, std::size_t sz ) {
		std::copy_n(cur,sz,ptr);
		std::advance(cur,sz);
	}

	ConstInputIterator current() const { return cur; }
private:
	ConstInputIterator cur;
};

template<typename ConstInputIterator>
iitr_stream<ConstInputIterator> make_istream(ConstInputIterator begin)
{
	return iitr_stream<ConstInputIterator>(begin);
}

/*
 * outputiterator_stream is a ouput_stream to OutputIterator
 */
template<typename OutputIterator>
class oitr_stream: public ostream {
public:
	oitr_stream( const oitr_stream& ) = default;
	oitr_stream( OutputIterator begin ) : cur(begin) {}
	~oitr_stream() {}
	void write( const char* ptr, std::size_t sz ) {
		cur = std::copy_n(ptr,sz,cur);
	}

	OutputIterator current() const { return cur; }
private:
	OutputIterator cur;
};

template<typename OutputIterator>
oitr_stream<OutputIterator> make_ostream(OutputIterator cur)
{
	return oitr_stream<OutputIterator>(cur);
}

/*
 * ocount_stream can be used for calculate the size of serialized data.
 */
class ocount_stream : public ostream {
public:
	ocount_stream(std::size_t off=0u) : sum(off) {}
	~ocount_stream() {}

	std::size_t size() const { return sum; }
	void write( const char*, std::size_t size ) { sum += size; }
private:
	std::size_t sum;
};

inline std::size_t streamed_size() { return 0u; }

template<typename T, typename... Args>
std::size_t streamed_size( const T& a, const Args&... args )
{
	ocount_stream stream;
	stream << a;
	return stream.size() + streamed_size(args...);
}

  // Use SFINAE to enable this only for types we marked as not
  // requiring endian conversion
  template<typename T,
	   typename std::enable_if<endian_traits<T>::convert == false>::type* = nullptr>
  istream& operator>>(istream& stream, T& dest) {
    stream.read(reinterpret_cast<char*>(&dest), sizeof(T));
    return stream;
  }

  // Use SFINAE to enable this only for types we marked as requiring
  // endian conversion
  template<typename T,
	   typename std::enable_if<endian_traits<T>::convert == true>::type* = nullptr>
  istream& operator>>(istream& stream, T& dest) {
    detail::endian_converter<sizeof(T)> conv;
    stream.read(conv.data.buf, sizeof(T));
    conv.betoh();
    std::memcpy(reinterpret_cast<char*>(&dest),
		conv.data.buf, sizeof(T));
    return stream;
  }

  // Use SFINAE to enable this only for types we marked as not
  // requiring endian conversion
  template<typename T,
	   typename std::enable_if<endian_traits<T>::convert == false>::type* = nullptr>
  ostream& operator<<(ostream& stream, const T& src) {
    stream.write(reinterpret_cast<const char*>(&src), sizeof(T));
    return stream;
  }

  // Use SFINAE to enable this only for types we marked as requiring
  // endian conversion
  template<typename T,
	   typename std::enable_if<endian_traits<T>::convert == true>::type* = nullptr>
  ostream& operator<<(ostream& stream, const T& src) {
    detail::endian_converter<sizeof(T)> conv;
    std::memcpy(conv.data.buf,
		reinterpret_cast<const char*>(&src), sizeof(T));
    conv.htobe();
    stream.write(conv.data.buf, sizeof(T));
    return stream;
  }

template<typename F, typename S>
istream& operator>> ( istream& stream, std::pair<F,S>& pair )
{
	stream >> pair.first >> pair.second;
	return stream;
}
template<typename F, typename S>
ostream& operator<< ( ostream& stream, const std::pair<F,S>& pair )
{
	stream << pair.first << pair.second;
	return stream;
}

template<typename T>
istream& operator>> ( istream& stream, std::complex<T>& cx )
{
	T r, i;
	stream >> r >> i;
	cx.real(r);
	cx.imag(i);
	return stream;
}
template<typename T>
ostream& operator<< ( ostream& stream, const std::complex<T>& cx )
{
	stream << cx.real() << cx.imag();
	return stream;
}

}

#endif
