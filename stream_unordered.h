/* -*- c++ -*-
 * stream_unordered.h
 *
 *  Created on: Mar 18, 2014
 *      Author: skudo
 */

#ifndef MUI_STREAM_UNORDERED_H
#define MUI_STREAM_UNORDERED_H

#include <unordered_map>
#include <unordered_set>
#include "stream.h"

namespace mui {
template<typename K, typename V>
inline istream& operator>>(istream& stream, std::unordered_map<K,V>& ret)
{
	std::size_t size;
	stream >> size;
	std::unordered_map<K,V> map;
	map.reserve(size);
	for( std::uint64_t i=0; i<size; ++i ) {
		std::pair<K,V> p;
		stream >> p;
		map.emplace(std::move(p));
	}
	ret.swap(map);
	return stream;
}
template<typename K, typename V>
inline ostream& operator<<(ostream& stream, const std::unordered_map<K,V>& map)
{
	stream << map.size();
	for( const auto& p : map ) stream << p;
	return stream;
}

template<typename K, typename V>
inline istream& operator>>(istream& stream, std::unordered_multimap<K,V>& ret)
{
	std::size_t size;
	stream >> size;
	std::unordered_multimap<K,V> map;
	map.reserve(size);
	for( std::uint64_t i=0; i<size; ++i ) {
		std::pair<K,V> p;
		stream >> p;
		map.emplace(std::move(p));
	}
	ret.swap(map);
	return stream;
}
template<typename K, typename V>
inline ostream& operator<<(ostream& stream, const std::unordered_multimap<K,V>& map)
{
	stream << map.size();
	for( const auto& p : map ) stream << p;
	return stream;
}

template<typename K>
inline istream& operator>>(istream& stream, std::unordered_set<K>& ret)
{
	std::size_t size;
	stream >> size;
	std::unordered_set<K> set;
	set.reserve(size);
	for( std::uint64_t i=0; i<size; ++i ) {
		K k;
		stream >> k;
		set.emplace(std::move(k));
	}
	ret.swap(set);
	return stream;
}

template<typename K>
inline ostream& operator<<(ostream& stream, const std::unordered_set<K>& set)
{
	stream << set.size();
	for( const auto& k : set ) stream << k;
	return stream;
}


template<typename K>
inline istream& operator>>(istream& stream, std::unordered_multiset<K>& ret)
{
	std::size_t size;
	stream >> size;
	std::unordered_multiset<K> set;
	set.reserve(size);
	for( std::uint64_t i=0; i<size; ++i ) {
		K k;
		stream >> k;
		set.emplace(std::move(k));
	}
	ret.swap(set);
	return stream;
}

template<typename K>
inline ostream& operator<<(ostream& stream, const std::unordered_multiset<K>& set)
{
	stream << set.size();
	for( const auto& k : set ) stream << k;
	return stream;
}
}

#endif
