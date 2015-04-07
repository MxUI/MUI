/* -*- c++ -*-
 * stream.h
 *
 *  Created on: Mar 18, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_H_
#define MUI_STREAM_H_

#include <memory>
#include <algorithm>

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

/*
 * container_steram is something like std::stringstream.
 * Its functionality is serialize & deserialize data to
 * container<char>.
 */
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

#define Makeopsh(TYPE,T)						\
	inline istream& operator>> ( istream& stream, TYPE& t )		\
	{								\
		stream.read(reinterpret_cast<char*>(&t),sizeof(t));	\
		return stream;						\
	}								\
	inline ostream& operator<< ( ostream& stream, TYPE t )		\
	{								\
		stream.write(reinterpret_cast<char*>(&t),sizeof(t));	\
		return stream;						\
	}

// use network-byte-order(big-endian)
#define Makeopshs(TYPE,SZ)						\
	inline istream& operator>> ( istream& stream, TYPE& t )		\
	{								\
		int##SZ##_t be;						\
		stream.read(reinterpret_cast<char*>(&be),sizeof(be));	\
		be = be##SZ##toh (be);					\
		std::memcpy(reinterpret_cast<char*>(&t),		\
		            reinterpret_cast<char*>(&be), sizeof(be));	\
		return stream;						\
	}								\
	inline ostream& operator<< ( ostream& stream, TYPE t )		\
	{								\
		int##SZ##_t be;						\
		std::memcpy(reinterpret_cast<char*>(&be),		\
		            reinterpret_cast<char*>(&t), sizeof(be));	\
		be = htobe##SZ (be);					\
		stream.write(reinterpret_cast<char*>(&be),sizeof(be));	\
		return stream;						\
	}

Makeopsh(char,8);
Makeopsh(signed char,8);
Makeopsh(unsigned char,8);

#ifndef MUI_IGNORE_ENDIAN
#if __BYTE_ORDER == __BIG_ENDIAN
#  define Makeopshi Makeopsh
#  if __FLOAT_WORD_ORDER == __BYTE_ORDER
#    define Makeopshf Makeopsh
#  else
#    define Makeopshf Makeopshs
#  endif
#else
#  define Makeopshi Makeopshs
#  if __FLOAT_WORD_ORDER == __BYTE_ORDER
#    define Makeopshf Makeopshs
#  else
#    define Makeopshf Makeopsh
#  endif
#endif
#else
#define Makeopshi Makeopsh
#define Makeopshf Makeopsh
#endif

Makeopshi(int16_t,16);
Makeopshi(uint16_t,16);
Makeopshi(int32_t,32);
Makeopshi(uint32_t,32);
Makeopshi(int64_t,64);
Makeopshi(uint64_t,64);
Makeopshf(float,32);
Makeopshf(double,64);

#undef Makeopshf
#undef Makeopshi
#undef Makeopshs
#undef Makeopsh

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
