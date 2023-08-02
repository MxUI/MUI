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
 * @file spatial_storage.h
 * @author S. Kudo
 * @date 21 April 2014
 * @brief Defines the spatial_storage data type.
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
	template<typename REGION, typename FOCUS, typename SAMPLER, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	query(const REGION& reg, const FOCUS& f, SAMPLER& s, ADDITIONAL && ... additional) const {
		using vec = std::vector<std::pair<point_type,typename SAMPLER::ITYPE> >;

		if( data_.empty() )
			return s.filter( f, virtual_container<typename SAMPLER::ITYPE,CONFIG>(vec(), std::vector<bool>()), additional... );
		if( !is_built() ) EXCEPTION(std::logic_error("MUI Error [spatial_storage.h]: Query error, "
				                                         "bin structure not built yet."));

		const vec& st = storage_cast<const vec&>(data_);

		return s.filter( f, virtual_container<typename SAMPLER::ITYPE,CONFIG>(st, bin_.query(reg)), additional...);
	}

	void build() {
		if( is_built() ) EXCEPTION(std::logic_error("MUI Error [spatial_storage.h]: Build error, cannot build bin structure twice."));
		if( !data_.empty() ) {
		  data_.apply_visitor(construct_{static_cast<void*>(&bin_)});
			is_bin_ = true;
		}
	}

	template<typename FOCUS, typename SAMPLER, typename ...ADDITIONAL>
	typename SAMPLER::OTYPE
	build_and_query_ts(const FOCUS& f, SAMPLER& s, ADDITIONAL && ... additional) {
		// this method is thread-safe. other methods are not.
		{
			std::unique_lock<std::mutex> lock(mutex_);
			if( !is_built() ) build();
		}

		return query(s.support(f, bin_.domain_size()).bbox(), f, s, additional...);
	}

	void insert( storage_t storage ) {
		destroy_if_bin_();
		if( !storage ) return;
		if( !data_ ) data_ = std::move(storage);
		else if( data_.which() == storage.which() ) data_.apply_visitor(insert_{storage});
		else EXCEPTION(bad_storage_id("MUI Error [spatial_storage.h]: Insert error. Type doesn't match."));
	}

	template<typename TYPE>
	const std::vector<std::pair<point_type,TYPE> >& return_data() {
		return storage_cast<const std::vector<std::pair<point_type,TYPE> >& >(data_);
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
