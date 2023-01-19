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
 * @file virtual_container.h
 * @author Y. H. Tang
 * @date 16 March 2014
 * @brief Provides a virtual container interface that is used to wrap around
 * data structures created by dynstorage.h and allows casting between storage
 * types.
 *
 */

#ifndef VIRTUAL_CONTAINER_H_
#define VIRTUAL_CONTAINER_H_

#include "../general/util.h"

namespace mui {

template<typename V, typename ARRAY>
struct index_iterator: public std::iterator<std::forward_iterator_tag,V>
{
	typedef std::iterator<std::forward_iterator_tag,V> iterator;

	typedef typename iterator::iterator_category iterator_category;
	typedef typename iterator::value_type value_type;
	typedef typename iterator::difference_type difference_type;
	typedef typename iterator::pointer pointer;
	typedef typename iterator::reference reference;

	index_iterator( ARRAY& array, std::size_t index=0 ): array_(array), index_(index) {}
	index_iterator() = default;
	index_iterator( const index_iterator& ) = default;
	index_iterator& operator=( const index_iterator& ) = default;

	value_type operator*() const { return array_[index_]; }
	reference operator*() { return array_[index_]; }
	const pointer operator->() const { return &(array_[index_]); }
	pointer operator->() { return &(array_[index_]); }


	index_iterator& operator++(){
		++index_;
		return *this;
	}
	index_iterator operator++(int){
		index_iterator prev = *this;
		operator++();
		return prev;
	}
	bool operator==( const index_iterator& rhs ) const { return index_ == rhs.index_; }
	bool operator!=( const index_iterator& rhs ) const { return index_ != rhs.index_; }
private:
	ARRAY& array_;
	std::size_t index_;
};


template < typename TYPE, typename CONFIG>
class virtual_container {
public:
	using elem_type      = std::pair<typename CONFIG::point_type,TYPE>;
	using container_type = std::vector<elem_type>;
	typedef index_iterator<const elem_type,const virtual_container> iterator;

	inline virtual_container( const container_type &container, std::vector<size_t> map ) :
		container_(container), map_(std::move(map)) {
	}
	virtual_container( const container_type &container, const std::vector<bool> &pred ) :
		container_(container) {
		for( size_t i = 0 ; i < pred.size() ; i++ ) if (pred[i]) map_.push_back(i);
	}

	// operator [] does no bound check
	inline const elem_type& operator [] ( size_t i ) const {
		return container_[ map_[i] ];
	}

	// at performs bound check
	inline const elem_type& at ( size_t i ) const {
		if ( i >= map_.size() ) typename CONFIG::EXCEPTION(std::out_of_range("MUI Error [virtual_container.h]: Out of range."));
		return container_[ map_[i] ];
	}

	inline iterator begin() const {
		iterator iter(*this);
		return iter;
	}
	inline iterator cbegin() const { return begin(); }

	inline iterator end()const{
		iterator iter(*this,size());
		return iter;
	}
	inline iterator cend() const { return end(); }

	inline size_t size() const { return map_.size(); }
protected:
	container_type const & container_;
	std::vector<size_t> map_;

};

template < typename TYPE, typename CONFIG>
virtual_container<TYPE,CONFIG> make_vc( const std::vector<typename CONFIG::point_type,TYPE> &container,
                                        const std::vector<size_t> &map )
{
	return virtual_container<TYPE,CONFIG>(container,map);
}

template < typename TYPE, typename CONFIG>
virtual_container<TYPE,CONFIG> make_vc( const std::vector<std::pair<typename CONFIG::point_type,TYPE> > &container,
                                        const std::vector<bool> &pred )
{
	return virtual_container<TYPE,CONFIG>(container,pred);
}

}

#endif /* VIRTUAL_CONTAINER_H_ */
