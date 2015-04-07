/*
 * bin.h
 *
 *  Created on: Apr 10, 2014
 *      Author: skudo
 */
#ifndef BIN_H
#define BIN_H

#include "geometry.h"
#include "config.h"

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

	std::vector<std::size_t> displs;
	// sorted[displs[i]] is the first element of the i-th bin

	point_type min, max;
	std::size_t n[CONFIG::D];
	typename CONFIG::REAL h;
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

		auto vol = max[0]-min[0];
		for( int i=1; i<D; ++i ) vol *= max[i]-min[i];

		h = std::pow(6*vol/val.size(),1.0/D); // about 6 points per bin
		h = h == 0.0 ? 1.0: h;
		std::size_t nn=1;
		for( int i=0; i<D; ++i ) {
			n[i] = std::ceil((max[i]-min[i])/h);
			n[i] = std::max( n[i], std::size_t(1) );
			nn *= n[i];
		}

		// make index
		std::vector<std::size_t> index(val.size()+1);
		std::vector<std::size_t> counts(nn);
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

private:
	bool initialize_query_( const geometry::box<CONFIG>& bx, int lda[], int lh[][2] ) const {
		bool broken = false;
		lda[0] = 1;
		for( int i=1; i<CONFIG::D; ++i ) lda[i] = lda[i-1]*n[i-1];
		for( int i=0; i<CONFIG::D; ++i ) {
			lh[i][0] = std::max(std::floor((bx.get_min()[i]-min[i])/h),typename CONFIG::REAL(0));
			lh[i][1] = std::min(std::ceil((bx.get_max()[i]-min[i])/h),typename CONFIG::REAL(n[i]));
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
			std::size_t d = int(std::floor((pt[i]-min[i])/h));
			d = std::min(d,n[i]-1); // the values may change here if (max[i]-min[i])/h is equal to a integer value.
			ret += m*d;
			m *= n[i];
		}
		return ret;
	}
};

}
#endif
