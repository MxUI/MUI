/* -*- c++ -*-
 * stream_ordered.h
 *
 *  Created on: Mar 26, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_ORDERED_H
#define MUI_STREAM_ORDERED_H

#include <map>
#include <set>
#include "stream.h"

namespace mui {
template<typename K, typename V>
inline istream& operator>>(istream& stream, std::map<K,V>& ret)
{
	std::size_t size;
	stream >> size;
	std::map<K,V> map;
	for( std::uint64_t i=0; i<size; ++i ) {
		std::pair<K,V> p;
		stream >> p;
		map.emplace(std::move(p));
	}
	ret.swap(map);
	return stream;
}
template<typename K, typename V>
inline ostream& operator<<(ostream& stream, const std::map<K,V>& map)
{
	stream << map.size();
	for( const auto& p : map ) stream << p;
	return stream;
}

template<typename K, typename V>
inline istream& operator>>(istream& stream, std::multimap<K,V>& ret)
{
	std::size_t size;
	stream >> size;
	std::multimap<K,V> map;
	for( std::uint64_t i=0; i<size; ++i ) {
		std::pair<K,V> p;
		stream >> p;
		map.emplace(std::ove(p));
	}
	ret.swap(map);
	return stream;
}
template<typename K, typename V>
inline ostream& operator<<(ostream& stream, const std::multimap<K,V>& map)
{
	stream << map.size();
	for( const auto& p : map ) stream << p;
	return stream;
}

template<typename K>
inline istream& operator>>(istream& stream, std::set<K>& ret)
{
	std::size_t size;
	stream >> size;
	std::set<K> set;
	for( std::uint64_t i=0; i<size; ++i ) {
		K k;
		stream >> k;
		set.emplace(std::move(k));
	}
	ret.swap(set);
	return stream;
}

template<typename K>
inline ostream& operator<<(ostream& stream, const std::set<K>& set)
{
	stream << set.size();
	for( auto& k : set ) stream << k;
	return stream;
}


template<typename K>
inline istream& operator>>(istream& stream, std::multiset<K>& ret)
{
	std::size_t size;
	stream >> size;
	std::multiset<K> set;
	for( std::uint64_t i=0; i<size; ++i ) {
		K k;
		stream >> k;
		set.emplace(std::move(k));
	}
	ret.swap(set);
	return stream;
}

template<typename K>
inline ostream& operator<< ( ostream& stream, const std::multiset<K>& set )
{
	stream << set.size();
	for( const auto& k : set ) stream << k;
	return stream;
}
}

#endif
