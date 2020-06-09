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
 * @file stream_unordered.h
 * @author S. Kudo
 * @date 18 March 2014
 * @brief Defines the stream in/out for the unordered std::unordered_map.
 * data type.
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
