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
******************************************************************************/

/**
 * @file bin.h
 * @author S. Kudo
 * @date 10 April 2014
 * @brief Structures and methods to create an underlying binning structure
 * for data received through an interface.
 */

#ifndef BIN_H
#define BIN_H

#include "../geometry/geometry.h"
#include "../config.h"

namespace mui {

namespace {
template<typename INT>
struct count_iterator : std::iterator<std::random_access_iterator_tag,INT,INT,INT,INT&>{
	INT cur_;
	count_iterator() {}
	count_iterator(const count_iterator&) = default;
	count_iterator( INT cur ): cur_(cur) {}

	count_iterator& operator++(){
		++cur_;
		return *this;
	}
	count_iterator& operator--(){
		--cur_;
		return *this;
	}
	count_iterator operator++(int){
		count_iterator ret = *this;
		operator++();
		return ret;
	}
	count_iterator operator--(int){
		count_iterator ret = *this;
		operator--();
		return ret;
	}
	count_iterator  operator+(INT n) const { return count_iterator(cur_+n); }
	count_iterator& operator+=(INT n) { cur_ +=n; return *this; }
	count_iterator  operator-(INT n) const { return count_iterator(cur_-n); }
	count_iterator& operator-(INT n) { cur_ -= n; return *this; }
	INT operator-(count_iterator&rhs) const { return cur_-rhs.cur_; }
	
	INT& operator*() { return cur_; }
	const INT& operator*() const { return cur_; }

	bool operator!=( count_iterator rhs ) const {
		return cur_ != rhs.cur_;
	}
	bool operator==( count_iterator rhs ) const {
		return !operator!=(rhs);
	}
};

template<typename INT>
count_iterator<INT>  operator+(INT n, count_iterator<INT> rhs) { return count_iterator<INT>(rhs.cur_+n); }
template<int depth> struct set_map_{
	static void apply( int offset, int lda[], int lh[][2],
	                   const std::vector<std::size_t>& displs, std::vector<std::size_t>& map ) {
		for( int i=lh[depth][0]; i<lh[depth][1]; ++i )
			set_map_<depth-1>::apply( offset+i*lda[depth], lda, lh, displs, map );
	}
};
template<> struct set_map_<0>{
	static void apply( int offset, int[] , int lh[][2],
	                   const std::vector<std::size_t>& displs, std::vector<std::size_t>& map ) {
		typedef count_iterator<std::size_t> iter;
		map.insert(map.end(), iter(displs[offset+lh[0][0]]), iter(displs[offset+lh[0][1]]));
	}
};
}

template<typename T, typename CONFIG>
struct bin_range;

template<typename T,typename CONFIG>
struct bin_iterator : std::iterator<std::forward_iterator_tag,std::pair<typename CONFIG::point_type,T> > {
	using P = std::pair<typename CONFIG::point_type,T>;
	static const int D = CONFIG::D;
	bin_iterator( const bin_range<T,CONFIG>& range_ );
	explicit bin_iterator( const bin_range<T,CONFIG>& range_, int) : range(range_) { invalidate(); }
	bin_iterator( const bin_iterator& ) = default;
	bin_iterator& operator=( const bin_iterator& ) = default;

	inline bin_iterator& operator++();
	inline bin_iterator operator++( int ){
		bin_iterator other = *this;
		this->operator++();
		return other;
	}

	inline void invalidate();
	inline void validate();
	inline bool operator!=( const bin_iterator& rhs ) const { return index != rhs.index; }
	inline bool operator==( const bin_iterator& rhs ) const { return !(this->operator!=(rhs)); }

	inline const P& operator*() const;
	inline const P* operator->() const;

	const bin_range<T,CONFIG>& range;
	int count[D-1];
	std::size_t high;
	std::size_t index;
};

template<typename T, typename CONFIG>
struct bin_range {
	using P = std::pair<typename CONFIG::point_type,T>;
	static const int D = CONFIG::D;
	using iterator = bin_iterator<T,CONFIG>;
	iterator begin() const { return iterator(*this); }
	iterator end() const { return iterator(*this,0); }
	
	bin_range( int lda_[], int lh_[][2], const std::vector<std::size_t>& d_, const std::vector<P>& v_ )
		: displs(d_), value(v_) {
		for( int i=0; i<D; ++i ){
			lda[i] = lda_[i];
			lh[i][0] = lh_[i][0];
			lh[i][1] = lh_[i][1];
		}
	}
	bin_range(const bin_range&) = default;
	bin_range& operator=(const bin_range&) = default;

	int lda[D];
	int lh[D][2];
	const std::vector<std::size_t>& displs;
	const std::vector<P>& value;
};

template<typename T, typename CONFIG>
bin_iterator<T,CONFIG>::bin_iterator( const bin_range<T,CONFIG>& range_ ): range(range_)
{
	for( int i=0; i<D-1; ++i ) count[i] = range.lh[i+1][0];
	validate();
}
template<typename T, typename CONFIG>
bin_iterator<T,CONFIG>& bin_iterator<T,CONFIG>::operator++() 
{
	if( ++index == high ){
		int i=0;
		for( ; i<D-2; ++i ){
			if( ++count[i] != range.lh[i+1][1] ) goto VALIDATE;
			count[i] = range.lh[i+1][0];
		}
		++count[i];
	VALIDATE:
		validate();
	}
	return *this;
}
template<typename T, typename CONFIG>
void bin_iterator<T,CONFIG>::invalidate()
{ 
	std::size_t offset = 0;
	for( int i=0; i<D-1; ++i ) offset += range.lda[i+1]*(range.lh[i+1][1]-1);
	index = high  = range.displs[offset+range.lh[0][1]];
}
template<typename T, typename CONFIG>
inline void bin_iterator<T,CONFIG>::validate()
{
	if( count[D-2] == range.lh[D-1][1] ) return;
	std::size_t offset = 0;
	for( int i=0; i<D-1; ++i ) offset += range.lda[i+1]*count[i];
	index = range.displs[offset+range.lh[0][0]];
	high  = range.displs[offset+range.lh[0][1]];
}

template<typename T, typename CONFIG>
inline const typename bin_iterator<T,CONFIG>::P& bin_iterator<T,CONFIG>::operator*() const { return range.value[index]; }
template<typename T, typename CONFIG>
inline const typename bin_iterator<T,CONFIG>::P* bin_iterator<T,CONFIG>::operator->() const { return std::addressof(range.value[index]); }


template<typename CONFIG>
struct bin_t {
private:
	using point_type = typename CONFIG::point_type;
	static const int D = CONFIG::D;
	static const bool QUIET = CONFIG::QUIET;

	std::vector<std::size_t> displs;
	// sorted[displs[i]] is the first element of the i-th bin

	point_type min, max;
	std::size_t n[CONFIG::D];
	typename CONFIG::REAL h;
	using REAL = typename CONFIG::REAL;
public:
	template<typename T>
	bin_t( std::vector<std::pair<point_type,T> >& val ){
		if( val.empty() ){
			displs.resize(2,0);
			h = 1.0;
			return;
		}

		// calculate h & n
		min = max = val[0].first;
		for( std::size_t i=1; i<val.size(); ++i ){
			point_type p = val[i].first;
			for( std::size_t i=0; i<D; ++i ) {
				min[i] = std::min(min[i],p[i]);
				max[i] = std::max(max[i],p[i]);
			}
		}

		size_t zero_count=0;
		REAL vol = std::abs(max[0]-min[0]);
		if(almost_equal(vol, static_cast<REAL>(0))) { // check if first dimension is zero size, if so set to 1
			vol = 1.0;
			zero_count++;
		}

		REAL vol_multi = vol;
		for( int i=1; i<D; ++i ){
			vol_multi = max[i]-min[i];
			if(almost_equal(vol_multi, static_cast<REAL>(0))) { // check if other dimensions are zero size, if so set them to 1
				vol_multi = static_cast<REAL>(1);
				zero_count++;
			}
			vol *= vol_multi;
		}

		if (zero_count == D) // if each dimension was actually zero (rather than just a subset) then set vol to zero
			vol = static_cast<REAL>(0);

		h = std::pow(static_cast<REAL>(6)*vol/static_cast<REAL>(val.size()),1.0/D); // about 6 points per bin

		if(almost_equal(h, static_cast<REAL>(0))){ // if h is still zero (only in the case of all dimensions being zero) then warn the user as this may be a problem
			h = static_cast<REAL>(1); // in this special case set h to 1 arbitrarily so bins work numerically
			if(val.size() > 1 && !QUIET)
				std::cout << "MUI Warning [bin.h]: Bin support size fixed to 1.0, check interface dimensionality or problem decomposition." << std::endl;
		}

		std::size_t nn=1;
		for( int i=0; i<D; ++i ) {
			n[i] = static_cast<size_t>(std::ceil((max[i]-min[i])/h));
			n[i] = static_cast<size_t>(std::max( n[i], std::size_t(1) ));
			nn *= n[i];
		}

		// make index
		std::vector<std::size_t> index(val.size()+1,0);
		std::vector<std::size_t> counts(nn,0);
		for( std::size_t i=0; i<val.size(); ++i ) {
			index[i] = get_index_(val[i].first);
			++counts[index[i]];
		}
		displs.resize(nn+1,0); // add 1 for sentinel
		std::partial_sum(counts.begin(),counts.end(), displs.begin()+1);

		counts = displs;
		std::vector<std::pair<point_type,T> > v(val.size());
		for( std::size_t i=0; i<val.size(); ++i ) v[counts[index[i]]++] = val[i];
		v.swap(val);
	}

	std::vector<std::size_t> query( const geometry::box<CONFIG>& bx ) const {
		std::vector<std::size_t> map;
		int lda[D];
		int lh[D][2];
		if( initialize_query_(bx,lda,lh) ) return map;
		map.reserve(lda[D-1]*12);
		set_map_<D-1>::apply( 0, lda, lh, displs, map );
		return map;
	}

	template<typename T>
	bin_range<T,CONFIG> query2( const geometry::box<CONFIG>& bx, const std::vector<std::pair<point_type,T> >& v ) const {
		int lda[D];
		int lh[D][2];
		initialize_query_(bx,lda,lh);
		return bin_range<T,CONFIG>{lda,lh,displs,v};
	}

	REAL domain_size() {
		REAL dim_size = norm(max-min);
		// Special case if domain only contains a single point
		if(dim_size == 0) dim_size = 1.0;
		return dim_size;
	}

private:
	bool initialize_query_( const geometry::box<CONFIG>& bx, int lda[], int lh[][2] ) const {
		bool broken = false;
		lda[0] = 1;
		for( int i=1; i<CONFIG::D; ++i ) lda[i] = lda[i-1]*n[i-1];
		for( int i=0; i<CONFIG::D; ++i ) {
			lh[i][0] = static_cast<int>(std::max(std::floor((bx.get_min()[i]-min[i])/h),typename CONFIG::REAL(0)));
			lh[i][1] = static_cast<int>(std::min(std::ceil((bx.get_max()[i]-min[i])/h),typename CONFIG::REAL(n[i])));
			if( lh[i][0] >= lh[i][1] ){
				lh[i][0] = lh[i][1];
				broken = true;
			}
		}
		return broken;
	}

	inline std::size_t get_index_( const point_type& pt ) const {
		std::size_t m = 1, ret=0;
		for( int i=0; i<D; ++i ) {
			std::size_t d = static_cast<size_t>((std::floor((pt[i]-min[i])/h)));
			d = static_cast<size_t>(std::min(d,n[i]-1)); // the values may change here if (max[i]-min[i])/h is equal to a integer value.
			ret += m*d;
			m *= n[i];
		}
		return ret;
	}
};

}
#endif
