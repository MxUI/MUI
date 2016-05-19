/* -*- c++ -*-
 * spatial_storage.h
 *
 *  Created on: Apr 21, 2014
 *      Author: skudo
 */

#define SPATIAL_STORAGE_H

#include <exception>
#include <mutex>
#include "dynstorage.h"
#include "virtual_container.h"

namespace mui {

template<typename BIN, typename STORAGE, typename CONFIG>
class spatial_storage {
private: // type definitions
	using storage_t = STORAGE;
public:
	using point_type = typename CONFIG::point_type;
	using EXCEPTION = typename CONFIG::EXCEPTION;

private: // functors
	struct insert_ { 
		template<typename T> void operator()( T& t ){
			T& rhs = storage_cast<T&>(st);
			t.insert(t.end(), std::make_move_iterator(rhs.begin()), std::make_move_iterator(rhs.end()));
		}
		storage_t st;
	};
	struct construct_ {
		template<typename T> void operator()(T& t){ ::new(ptr) BIN(t); }
		void* ptr;
	};
public:
	spatial_storage() noexcept {
		is_bin_ = false;
	}
	~spatial_storage() {
		destroy_if_bin_();
	}
	spatial_storage( const spatial_storage& rhs ): data_(rhs.data_) {
		is_bin_.store(rhs.is_bin_.load());
		if( rhs.is_built() ) ::new(&bin_) BIN(rhs.bin_);
	}
	spatial_storage( spatial_storage&& rhs ) noexcept : data_(std::move(rhs.data_)) {
		is_bin_.store(rhs.is_bin_.load());
		if( rhs.is_built() ) ::new(&bin_) BIN(std::move(rhs.bin_));
		rhs.destroy_if_bin_();
	}
	spatial_storage( storage_t rhs ): data_(std::move(rhs)) {
		is_bin_ = false;
	}
	spatial_storage& operator=( spatial_storage rhs ) {
		rhs.swap(*this);
		return *this;
	}
	
	void swap( spatial_storage& rhs ) noexcept(noexcept(BIN(std::move(rhs.bin_)))) {
		data_.swap(rhs.data_);
		if( is_bin_ && rhs.is_bin_ ) bin_.swap(rhs.bin_);
		else if( is_bin_ && !rhs.is_bin_ ) {
			::new(&(rhs.bin_)) BIN(std::move(bin_));
			rhs.is_bin_ = true;
			destroy_if_bin_();
		} else if( !is_bin_ && rhs.is_bin_ ) {
			::new(&(bin_)) BIN(std::move(rhs.bin_));
			is_bin_ = true;
			rhs.destroy_if_bin_();
		}
	}
	

	// no lock internally; it's job of uniface because
	// we need more large lock if there is no rw-lock.
	template<typename REGION, typename FOCUS, typename SAMPLER>
	typename SAMPLER::OTYPE 
	query(const REGION& reg, const FOCUS& f, const SAMPLER& s) const {
		using vec = std::vector<std::pair<point_type,typename SAMPLER::ITYPE> >;
		if( data_.empty() ) 
			return s.filter(f, virtual_container<typename SAMPLER::ITYPE,CONFIG>(vec(),std::vector<bool>()));
		if( !is_built() ) EXCEPTION(std::logic_error("spatial storage: query error. "
		                                              "not builded bin. Internal data was corrupsed."));
		const vec& st = storage_cast<const vec&>(data_);
		return s.filter(f,virtual_container<typename SAMPLER::ITYPE,CONFIG>(st,bin_.query(reg)));
	}
	template<typename REGION, typename FOCUS, typename SAMPLER>
	typename SAMPLER::OTYPE 
	query2(const REGION& reg, const FOCUS& f, const SAMPLER& s) const {
		using vec = std::vector<std::pair<point_type,typename SAMPLER::ITYPE> >;
		if( data_.empty() ) 
			return s.filter(f, vec());
		if( !is_built() ) EXCEPTION( std::logic_error("spatial storage: query error. "
		                                               "not builded bin. Internal data was corrupsed.") );
		const vec& st= storage_cast<const vec&>(data_);
		return s.filter(f,bin_.query2(reg,st));
	}
	void build() {
		if( is_built() ) EXCEPTION(std::logic_error("spatial storage: build error. cannot build twice."));
		if( !data_.empty() ){
			data_.apply_visitor(construct_{(void*)&bin_});
			is_bin_ = true;
		}
	}
	template<typename REGION, typename FOCUS, typename SAMPLER>
	typename SAMPLER::OTYPE 
	build_and_query_ts(const REGION& reg, const FOCUS& f, const SAMPLER& s) {
		// this method is thread-safe. other methods are not.
		{
			std::unique_lock<std::mutex> lock(mutex_);
			if( !is_built() ) build();
		}
		return query(reg, f, s);
	}
	void insert( storage_t storage ) {
		destroy_if_bin_();
		if( !storage ) return;
		if( !data_ ) data_ = std::move(storage);
		else if( data_.which() == storage.which() ) data_.apply_visitor(insert_{storage});
		else EXCEPTION(bad_storage_id("spatial storage: insert error. Type doesn't match."));
	}
	
	bool is_built() const { return is_bin_; }
	bool empty() const { return data_.empty(); }
private:
	void destroy_if_bin_() { 
		if( is_built() ) {
			bin_.~BIN(); 
			is_bin_ = false;
		}
	}

	storage_t data_;

	std::atomic<bool> is_bin_;
	union {
		char nodata_ = '\0';
		BIN bin_;
	};
	mutable std::mutex mutex_;
};
}

#define SPATIAL_STORAGE_H
