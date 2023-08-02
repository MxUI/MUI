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
 * @file dynstorage.h
 * @author S. Kudo
 * @date 10 February 2014
 * @brief Implementation of a compound dynamic data structure used throughout
 * MUI.
 */

#ifndef DYNSTORAGE_H_
#define DYNSTORAGE_H_

#include <vector>
#include <utility>
#include <cstdint>
#include <typeinfo>
#include <type_traits>

#include "stream.h"
#include "../general/util.h"

namespace mui {

struct bad_storage_id: std::runtime_error {
	bad_storage_id( const char* err ): std::runtime_error(err) {}
};
struct bad_storage_cast: std::bad_cast {
	virtual const char* what() const noexcept { return "MUI Error [dynstorage.h]: Storage error, bad cast."; }
};


namespace {
// meta helper functions

// get_typeid_: this value is based on the order of Variadic template argument.
// get_typeid_<int, int, double, float, char, etc...>::value==0,
// get_typeid_<double, int, double, float, char, etc...>::value == 1, and so on
template<typename T, typename... Args> struct get_typeid_;
template<typename T, typename Head, typename... Tail> struct get_typeid_<T,Head,Tail...> {
	typedef typename std::remove_reference<typename std::remove_cv<T>::type>::type type;
	static constexpr std::int32_t value =
		std::conditional<std::is_same<type,Head>::value,
		                 std::integral_constant<std::int32_t,0>,
		                 std::integral_constant<std::int32_t,1+get_typeid_<T,Tail...>::value> >::type::value;
};
template<typename T> struct get_typeid_<T> {
	static constexpr std::int32_t value = 0;
};

// build i-th element and apply functor
template<std::int32_t i, typename R, typename... Args> struct make_value_;
template<std::int32_t i, typename R, typename Head, typename... Tail > struct make_value_<i,R,Head,Tail...> {
	template<typename F>
	static R apply( std::int32_t which, F f ) {
		if( i == which ) {
			Head head;
			return f(std::move(head));
		}
		else return make_value_<i+1,R,Tail...>::apply(which,f);
	}
};
template<std::int32_t i,typename R> struct make_value_<i,R> {
	template<typename F>
	static R apply( int, F ) { throw bad_storage_id("MUI Error [dynstorage.h]: Storage error, bad type id."); }
};


template<typename Head, typename... Tail> struct get_head_ { typedef Head type; };
template<typename... Types>
using get_head_t_ = typename get_head_<Types...>::type;

// type dispatcher. 
template<id_t i, typename R, typename... Args> struct apply_visitor_impl_;
template<id_t i, typename R, typename Head, typename... Tail>
struct apply_visitor_impl_<i,R,Head,Tail...> {
	template<typename F>
	static R apply( id_t which, void* content, F& f ) { // this applies the content to f
		if( i == which ) return f(*static_cast<Head*>(content));
		else return apply_visitor_impl_<i+1,R,Tail...>::apply(which,content,f);
	}
	template<typename F>
	static R apply( id_t which, const void* content, F& f ) {
		if( i == which ) return f(*static_cast<const Head*>(content));
		else return apply_visitor_impl_<i+1,R,Tail...>::apply(which,content,f);
	}
	template<typename F>
	static R applym( id_t which, void* content, F& f ) { // move the content to f
		if( i == which ) return f(std::move(*static_cast<Head*>(content)));
		else return apply_visitor_impl_<i+1,R,Tail...>::applym(which,content,f);
	}
};
template<id_t i, typename R> struct apply_visitor_impl_<i,R> {
	template<typename F>
	static R apply( id_t, const void*, F& ) { throw bad_storage_id("MUI Error [dynstorage.h]: Storage error, bad id."); }
	template<typename F>
	static R applym( id_t, void*, F& )      { throw bad_storage_id("MUI Error [dynstorage.h]: Storage error, bad id."); }
};
template<typename R,typename... Types> using applyer_ = apply_visitor_impl_<0,R,Types...>;
}

// storage can hold any type of Types... and it has value semantics.
// feature of storage
//   copy&construct as a value
//   member functions: swap, clear, empty
//   apply_visitor: user can get the actual type of storage by visiting it.
// ref: boost::any, boost::type_erasure, boost::variant
template<typename... Types>
struct storage {
private: // functors
	struct deleter_ { template<typename T> void operator()( T& t ){ delete std::addressof(t); } };
	struct cloner_ { template<typename T> void* operator()( const T& t ){ return static_cast<void*>(new T(t)); } };

private: // internal typedefs
	using Head_ = get_head_t_<Types...>;

public:
	using id_t = std::int32_t;
	storage() noexcept : which_(bad_id), content_(0) {}
	storage( const storage& rhs ) : which_(rhs.which_), content_(rhs.content_? rhs.clone_() : 0) {}
	storage( storage&& rhs ) noexcept : which_(rhs.which_), content_(rhs.content_) {
		rhs.which_ = bad_id;
		rhs.content_ = 0;
	}
	~storage() {
		if(!empty()) this->apply_visitor(deleter_());
		which_ = bad_id;
		content_ = 0;
	}
	storage& operator=( storage rhs ) {
		rhs.swap(*this);
		return *this;
	}

	template<typename ValueType>
	explicit storage( const ValueType& value)
		: which_(get_typeid_<ValueType,Types...>::value),
		  content_(new ValueType(value)) {
		static_assert(get_typeid_<ValueType,Types...>::value != bad_id,
		              "MUI Error [dynstorage.h]: Storage error, unsupported type. Please add type to type_list.");
	}
	template<typename ValueType>
	explicit storage(ValueType&& value,
	        typename std::enable_if<!std::is_same<storage&,ValueType>::value>::type* = 0,
	        typename std::enable_if<!std::is_const<ValueType>::value>::type* = 0 )
		: which_(get_typeid_<ValueType,Types...>::value),
		  content_(new typename std::decay<ValueType>::type(std::move(value))) {
		static_assert(get_typeid_<ValueType,Types...>::value != bad_id,
				"MUI Error [dynstorage.h]: Storage error, unsupported type. Please add type to type_list.");
	}

	template<typename ValueType>
	storage& operator=(ValueType&& rhs) {
		storage(std::forward<ValueType>(rhs)).swap(*this);
		return *this;
	}

	void swap( storage& rhs ) noexcept {
		std::swap(which_,rhs.which_);
		std::swap(content_,rhs.content_);
	}

private:
	void* clone_() const { return apply_visitor(cloner_()); }
public:
	void clear() noexcept { storage().swap(*this); }
	bool empty()    const noexcept { return content_ == 0; }
	explicit operator bool() const noexcept { return !empty(); }

	static const id_t bad_id = sizeof...(Types);
	id_t which() const { return which_; }

	template<typename F, typename R_=typename std::result_of<F(const Head_&)>::type>
	R_ apply_visitor( F f ) const& {
		return applyer_<R_,Types...>::apply(which_, content_, f);
	}
	template<typename F, typename R_=typename std::result_of<F(Head_&)>::type>
	R_ apply_visitor( F f ) & {
		return applyer_<R_,Types...>::apply(which_, content_, f);
	}
	template<typename F, typename R_=typename std::result_of<F(Head_&&)>::type>
	R_ apply_visitorm( F f ) && {
		return applyer_<R_,Types...>::applym(which_, content_, f);
	}


private:
	id_t which_ = bad_id;
	void* content_ = 0;

	template<typename ValueType, typename... Types2>
	friend ValueType* storage_cast(storage<Types2...>*);
};


template<typename... Args>
inline void swap(storage<Args...>& lhs, storage<Args...>& rhs )
{
	lhs.swap(rhs);
}

template<typename ValueType, typename... Args>
ValueType* storage_cast(storage<Args...>* obj )
{
	return obj && obj->which() == get_typeid_<ValueType,Args...>::value ?
		static_cast<ValueType*>(obj->content_) : 0;
}

template<typename ValueType, typename... Args>
const ValueType* storage_cast(const storage<Args...>* obj )
{
	return storage_cast<ValueType>(const_cast<storage<Args...>*>(obj));
}

template<typename ValueType, typename... Args>
ValueType storage_cast( storage<Args...>& obj ){
	typedef typename std::remove_reference<ValueType>::type nonref;
	nonref* result = storage_cast<nonref>(&obj);
	if(!result) throw bad_storage_cast();

	typedef typename std::add_lvalue_reference<ValueType>::type ref_type;
	return static_cast<ref_type>(*result);
}
template<typename ValueType, typename... Args>
ValueType storage_cast( const storage<Args...>& obj ){
	typedef typename std::remove_reference<ValueType>::type nonref;
	return storage_cast<const nonref&>(const_cast<storage<Args...>&>(obj));
}

template<typename ValueType, typename... Args>
inline ValueType&& storage_cast(storage<Args...>&& obj)
{
	static_assert(std::is_rvalue_reference<ValueType&&>::value
	              || std::is_const<typename std::remove_reference<ValueType>::type>::type,
		"MUI Error [dynstorage.h]: This type of cast is not supported.");
	return storage_cast<ValueType&&>(obj);
}


// (de)serialization
namespace {
struct serializer {
	template<typename T> void operator()( const T& t ){ ost << t; }
	ostream& ost;
};
}
template<typename... Args>
ostream& operator<<( ostream& stream, const storage<Args...>& st )
{
	stream << st.which();
	if( !st.empty() ) st.apply_visitor(serializer{stream});
	return stream;
}

namespace {
template<typename Storage>
struct deserializer {
	template<typename T> Storage operator()( T&& t ){
		ist >> t;
		return Storage(std::move(t));
	}
	istream& ist;
};
}
template<typename... Args>
istream& operator>>( istream& stream, storage<Args...>& st )
{
	using type=storage<Args...>;
	typename type::id_t which;
	stream >> which;
	st = (which == type::bad_id) ? type() :
		make_value_<0,type,Args...>::apply(which,deserializer<type>{stream});
	return stream;
}

}
#endif
