/*****************************************************************************
* Multiscale Universal Interface Code Coupling Library                       *
*                                                                            *
* Copyright (C) 2019 Y. H. Tang, S. Kudo, X. Bian, Z. Li, G. E. Karniadakis, *
*                    S. M. Longshaw, A. Skillen                              *
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
 * @file uniface.h
 * @author S. Kudo
 * @date 11 February 2014
 * @brief Provides the majority of the useful functionality for MUI,
 * including all fetch, commit and push functions.
 *
 * The majority of the user interface for MUI is defined here, as are all of
 * the data structures.
 */

#ifndef UNIFACE_H_
#define UNIFACE_H_

#include "general/util.h"
#include "communication/comm.h"
#include "communication/comm_factory.h"
#include "config.h"
#include "storage/dynstorage.h"
#include "storage/spatial_storage.h"
#include "communication/lib_dispatcher.h"
#include "communication/message/message.h"
#include "communication/message/reader_variable.h"
#include "storage/stream_vector.h"
#include "storage/stream_unordered.h"
#include "storage/stream_string.h"
#include "storage/bin.h"
#include "storage/stream.h"

#ifdef PYTHON_BINDINGS
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

namespace mui {

template<typename CONFIG = default_config>
class uniface {
public:
	// public typedefs (see config.h for descriptions)
	static const int  D     = CONFIG::D;
	static const bool FIXEDPOINTS = CONFIG::FIXEDPOINTS;
	static const bool QUIET = CONFIG::QUIET;

	using REAL = typename CONFIG::REAL;
	using point_type = typename CONFIG::point_type;
	using time_type = typename CONFIG::time_type;
	using iterator_type = typename CONFIG::iterator_type;
	using data_types = typename CONFIG::data_types;
	using span_t = geometry::any_shape<CONFIG>;
private:
	// meta functions to split tuple and add vector<pair<point_type,_1> >
	template<typename T> struct add_vp_ { using type = std::vector<std::pair<point_type,T> >; };
		template<typename T> struct def_storage_;
	template<typename... TYPES> struct def_storage_<type_list<TYPES...> >{
		using type = storage<typename add_vp_<TYPES>::type...>;
	};

	template<typename T> struct add_vi_ { using type = std::vector<std::pair<size_t, T> >; };
	template<typename T> struct def_storage_raw_;
	template<typename... TYPES> struct def_storage_raw_<type_list<TYPES...> >{
		using type = storage<typename add_vi_<TYPES>::type...>;
	};

	template<typename T> struct def_storage_single_;
	template<typename... TYPES> struct def_storage_single_<type_list<TYPES...> >{
		using type = storage<TYPES...>;
	};

	// internal typedefinitions for full frame
	using storage_t = typename def_storage_<data_types>::type;
	using spatial_t = spatial_storage<bin_t<CONFIG>,storage_t,CONFIG>;
	using frame_type = std::unordered_map<std::string, storage_t>;
	using bin_frame_type = std::unordered_map<std::string, spatial_t>;
	// internal typdefinitions for data values only (static points)
	using storage_raw_t = typename def_storage_raw_<data_types>::type;
	using frame_raw_type = std::unordered_map<std::string, storage_raw_t>;
	// internal typedefinitions for single value
	using storage_single_t = typename def_storage_single_<data_types>::type;

	struct peer_state {
		peer_state() : disable_send(false), disable_recv(false), ss_stat_send(false), ss_stat_recv(false) {}

		using spans_type = std::map<std::pair<time_type,time_type>,span_t>;

		bool is_recving(time_type t, const span_t& s) const {
			return scan_spans_(t,s,recving_spans);
		}
		
		void set_recving( time_type start, time_type end, span_t s ) {
			recving_spans.emplace(std::make_pair(start,end),std::move(s));
		}
		
		bool is_sending(time_type t, const span_t& s) const {
			return scan_spans_(t,s,sending_spans);
		}
		
		void set_sending(time_type start, time_type end, span_t s) {
			sending_spans.emplace(std::make_pair(start,end), std::move(s));
		}

		void set_pts(std::vector<point_type>& pts) {
			pts_ = pts;
		}

		const std::vector<point_type>& pts() const {
			return pts_;
		}
		
		void set_send_disable() {
			disable_send = true;
		}
		
		void set_recv_disable() {
			disable_recv = true;
		}
		
		bool is_send_disabled() const {
			return disable_send;
		}
		
		bool is_recv_disabled() const {
			return disable_recv;
		}

		void set_ss_send_status(bool status) {
			ss_stat_send = status;
		}

		bool ss_send_status() const {
			return ss_stat_send;
		}

		void set_ss_recv_status(bool status) {
			ss_stat_recv = status;
		}

		bool ss_recv_status() const {
			return ss_stat_recv;
		}

		time_type current_t() const { return latest_timestamp; }
		iterator_type current_it() const { return latest_subiter; }
		time_type next_t() const { return next_timestamp; }
		iterator_type next_it() const { return next_subiter; }
		void set_current_t( time_type t ) { latest_timestamp = t; }
		void set_current_sub( iterator_type i ) { latest_subiter = i; }
		void set_next_t( time_type t ) { next_timestamp = t; }
		void set_next_sub( iterator_type i ) { next_subiter = i; }
	private:
		bool scan_spans_(time_type t, const span_t& s, const spans_type& spans ) const {
			bool prefetched = false;
			auto end = spans.lower_bound(std::make_pair(t,t));
			if( spans.size() == 1 ) end = spans.end();

			for( auto itr = spans.begin(); itr != end; ++itr ) {
				if( t < itr->first.second || almost_equal(t, itr->first.second) ) {
					prefetched = true;
					if( collide(s,itr->second) ) return true;
				}
			}

			// if prefetched at t, but no overlap region, then return false;
			// otherwise return true;
			return !prefetched;
		}

		time_type latest_timestamp = std::numeric_limits<time_type>::lowest();
		iterator_type latest_subiter = std::numeric_limits<iterator_type>::lowest();
		time_type next_timestamp = std::numeric_limits<time_type>::lowest();
		iterator_type next_subiter = std::numeric_limits<iterator_type>::lowest();
		spans_type recving_spans;
		spans_type sending_spans;
		std::vector<point_type> pts_;
		std::unordered_map<std::string, storage_single_t> assigned_vals_;
		bool disable_send;
		bool disable_recv;
		bool ss_stat_send;
		bool ss_stat_recv;
	};

private: // data members
	std::unique_ptr<communicator> comm;
	dispatcher<message::id_type, std::function<void(message)> > readers;

	std::map<std::pair<time_type, iterator_type>, bin_frame_type> log;

	frame_type push_buffer;
	frame_raw_type push_buffer_raw;
	std::vector<point_type> push_buffer_pts;

	std::unordered_map<std::string, storage_single_t > assigned_values;

	std::vector<peer_state> peers;
	std::vector<bool> peer_is_sending;
	bool smart_send_set_ = true;
	time_type span_start = std::numeric_limits<time_type>::lowest();
	time_type span_timeout = std::numeric_limits<time_type>::lowest();
	span_t current_span;
	time_type recv_start = std::numeric_limits<time_type>::lowest();
	time_type recv_timeout = std::numeric_limits<time_type>::lowest();
	span_t recv_span;
	time_type memory_length = std::numeric_limits<time_type>::max();
	std::mutex mutex;
	bool initialized_pts_;
	size_t fixedPointCount_;
	time_type fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
	iterator_type fetch_i_hist_ = std::numeric_limits<iterator_type>::lowest();

public:
	uniface( const char URI[] ) : uniface( comm_factory::create_comm(URI, QUIET) ) {}
	uniface( std::string const &URI ) : uniface( comm_factory::create_comm(URI.c_str(), QUIET) ) {}
	uniface( communicator* comm_ ) : comm(comm_), initialized_pts_(false), fixedPointCount_(0) {
		using namespace std::placeholders;

		peers.resize(comm->remote_size());
		peer_is_sending.resize(comm->remote_size(), true);

		readers.link("timestamp", reader_variables<int32_t, std::pair<time_type,iterator_type> >(
					 std::bind(&uniface::on_recv_confirm, this, _1, _2)));
		readers.link("forecast", reader_variables<int32_t, std::pair<time_type,iterator_type>>(
					 std::bind(&uniface::on_recv_forecast, this, _1, _2)));
		readers.link("data", reader_variables<std::pair<time_type,iterator_type>, frame_type>(
					 std::bind(&uniface::on_recv_data, this, _1, _2)));
		readers.link("rawdata", reader_variables<int32_t, std::pair<time_type,iterator_type>, frame_raw_type>(
					 std::bind(&uniface::on_recv_rawdata, this, _1, _2, _3)));
		readers.link("points", reader_variables<int32_t, std::vector<point_type>>(
					 std::bind(&uniface::on_recv_points, this, _1, _2)));
		readers.link("assignedVals", reader_variables<std::string, storage_single_t>(
					 std::bind(&uniface::on_recv_assignedVals, this, _1, _2)));
		readers.link("receivingSpan", reader_variables<int32_t, time_type,time_type, span_t>(
		       std::bind(&uniface::on_recv_span, this, _1, _2, _3, _4)));
		readers.link("sendingSpan", reader_variables<int32_t, time_type,time_type, span_t>(
		       std::bind(&uniface::on_send_span, this, _1, _2, _3, _4)));
		readers.link("receivingDisable", reader_variables<int32_t>(
					 std::bind(&uniface::on_send_disable, this, _1)));
		readers.link("sendingDisable", reader_variables<int32_t>(
					 std::bind(&uniface::on_recv_disable, this, _1)));
	}

	uniface( const uniface& ) = delete;
	uniface& operator=( const uniface& ) = delete;

	/** \brief Announce the value \c value with the parameter \c attr
	* Useful if, for example, you wish to pass a parameter
	* rather than a field without an associated timestamp
	*/
	template<typename TYPE>
	void push( const std::string& attr, const TYPE& value ) {
		comm->send(message::make("assignedVals", attr, storage_single_t(TYPE(value))));
	}

	/** \brief Push data with tag "attr" to buffer
	* Push data with tag "attr" to bcuffer. If using CONFIG::FIXEDPOINTS=true,
	* data must be pushed in the same order that the points were previously pushed.
   	*/
	template<typename TYPE>
	void push( const std::string& attr, const point_type& loc, const TYPE& value ) {
		if( FIXEDPOINTS ) {
			// If this push is before first commit then build local points list
			if( !initialized_pts_ ) push_buffer_pts.emplace_back( loc );

			storage_raw_t& n = push_buffer_raw[attr];
			if( !n ) n = storage_raw_t(std::vector<std::pair<size_t,TYPE> >());
			storage_cast<std::vector<std::pair<size_t,TYPE> >&>(n).emplace_back( fixedPointCount_, value );

			// Increment counter for flat fixed point list
			fixedPointCount_++;
		} 
		else {
			storage_t& n = push_buffer[attr];
			if( !n ) n = storage_t(std::vector<std::pair<point_type,TYPE> >());
			storage_cast<std::vector<std::pair<point_type,TYPE> >&>(n).emplace_back( loc, value );
		}
	}

#ifdef PYTHON_BINDINGS
    template<typename TYPE>
    void push_many(const std::string& attr, const class py::array_t<REAL>& points,
                   const class py::array_t<TYPE>& values) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        auto values_arr = values.template unchecked<1>();
        assert(points_arr.shape(0) == values_arr.shape(0));
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            push<TYPE>(attr, p, values_arr(i));
        }
    }

    template<class SAMPLER, class TIME_SAMPLER>
    py::array_t<typename SAMPLER::OTYPE,py::array::c_style>
    fetch_many(const std::string& attr,const py::array_t<REAL,py::array::c_style> points, const time_type t,
        const SAMPLER &sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        py::array_t<typename SAMPLER::OTYPE,py::array::c_style> values(points_arr.shape(0));
        auto values_arr = values.template mutable_unchecked<1>();
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            values_arr(i)  = fetch(attr, p, t, sampler, t_sampler, barrier_enabled);
        }
        return values;
    }

    template<class SAMPLER, class TIME_SAMPLER>
    py::array_t<typename SAMPLER::OTYPE,py::array::c_style>
    fetch_many(const std::string& attr,const py::array_t<REAL,py::array::c_style> points, const time_type t,
        const iterator_type it, const SAMPLER &sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        py::array_t<typename SAMPLER::OTYPE,py::array::c_style> values(points_arr.shape(0));
        auto values_arr = values.template mutable_unchecked<1>();
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            values_arr(i)  = fetch(attr, p, t, it, sampler, t_sampler, barrier_enabled);
        }
        return values;
    }

    template<class SAMPLER, class TIME_SAMPLER, class ALGORITHM>
    py::array_t<typename SAMPLER::OTYPE,py::array::c_style>
    fetch_many(const std::string& attr,const py::array_t<REAL,py::array::c_style> points, const time_type t,
        const SAMPLER &sampler, const TIME_SAMPLER &t_sampler, const ALGORITHM &algorithm, bool barrier_enabled = true) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        py::array_t<typename SAMPLER::OTYPE,py::array::c_style> values(points_arr.shape(0));
        auto values_arr = values.template mutable_unchecked<1>();
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            values_arr(i)  = fetch(attr, p, t, sampler, t_sampler, algorithm, barrier_enabled);
        }
        return values;
    }

    template<class SAMPLER, class TIME_SAMPLER, class ALGORITHM>
    py::array_t<typename SAMPLER::OTYPE,py::array::c_style>
    fetch_many(const std::string& attr,const py::array_t<REAL,py::array::c_style> points, const time_type t,
        const iterator_type it, const SAMPLER &sampler, const TIME_SAMPLER &t_sampler, const ALGORITHM &algorithm, bool barrier_enabled = true) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        py::array_t<typename SAMPLER::OTYPE,py::array::c_style> values(points_arr.shape(0));
        auto values_arr = values.template mutable_unchecked<1>();
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            values_arr(i)  = fetch(attr, p, t, it, sampler, t_sampler, algorithm, barrier_enabled);
        }
        return values;
    }


    template<typename TYPE, class TIME_SAMPLER>
    py::array_t<REAL, py::array::c_style>
    fetch_points_np(const std::string& attr, const time_type t, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, TYPE test_value = static_cast<TYPE>(0)) {
        std::vector<point_type> points = fetch_points<TYPE>(attr, t, t_sampler, barrier_enabled);
        size_t n = points.size();
        test_value += 1;
        py::array_t<REAL, py::array::c_style> points_np({n, static_cast<size_t>(D)});
        auto points_np_arr = points_np.template mutable_unchecked<2>();
        for (std::size_t i = 0; i < n; i++)
            for (std::size_t j = 0; j < D; j++)
                points_np_arr(i,j) = (points[i].data())[j];
        return points_np;
    }

    template<typename TYPE, class TIME_SAMPLER>
    py::array_t<REAL, py::array::c_style>
    fetch_points_np(const std::string& attr, const time_type t, const iterator_type it, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, TYPE test_value = static_cast<TYPE>(0)) {
        std::vector<point_type> points = fetch_points<TYPE>(attr, t, it, t_sampler, barrier_enabled);
        size_t n = points.size();
        test_value += 1;
        py::array_t<REAL, py::array::c_style> points_np({n, static_cast<size_t>(D)});
        auto points_np_arr = points_np.template mutable_unchecked<2>();
        for (std::size_t i = 0; i < n; i++)
            for (std::size_t j = 0; j < D; j++)
                points_np_arr(i,j) = (points[i].data())[j];
        return points_np;
    }
#endif

    /** \brief Fetch a single parameter from the interface
	* Overloaded \c fetch to fetch a single parameter of name \c attr.
	* There is no barrier on this fetch as there is no time associated
	* with the value.
	*/
	template<typename TYPE>
	TYPE fetch( const std::string& attr ) {
		storage_single_t& n = assigned_values[attr];
		if( !n ) return TYPE();
		return storage_cast<TYPE&>(n);
	}

	/** \brief Fetch from the interface, blocking with barrier at time=t
    */
	template<class SAMPLER, class TIME_SAMPLER, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr,const point_type& focus, const time_type t,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true,
		   ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t
		if( fetch_t_hist_ != t && barrier_enabled )
			barrier(t_sampler.get_upper_bound(t));

		fetch_t_hist_ = t;

		std::vector<std::pair<std::pair<time_type,iterator_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
			  	  	  	  	  	  	  	  	  	  	  	   std::numeric_limits<iterator_type>::lowest());

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
			  	  	  	  	  	  	  	  	  	  	  	   std::numeric_limits<iterator_type>::lowest());
		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( start->first, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return t_sampler.filter(t, v);
	}

	/** \brief Fetch from the interface, blocking with barrier at time=t,it
	*/
	template<class SAMPLER, class TIME_SAMPLER, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr,const point_type& focus, const time_type t, const iterator_type it,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true,
		   ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t,iteration=it
		if((fetch_t_hist_ != t || fetch_i_hist_ != it) && barrier_enabled)
			barrier(t_sampler.get_upper_bound(t),t_sampler.get_upper_bound(it));

		fetch_t_hist_ = t;
		fetch_i_hist_ = it;

		std::vector<std::pair<std::pair<time_type,iterator_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
														   t_sampler.get_lower_bound(it)-threshold(it));

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
														   t_sampler.get_upper_bound(it)+threshold(it));
		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( start->first, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return t_sampler.filter(std::make_pair(t,it), v);
	}

	/** \brief Fetch from the interface with coupling algorithms, blocking with barrier at time=t
	*/
	template<class SAMPLER, class TIME_SAMPLER, class COUPLING_ALGO, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr,const point_type& focus, const time_type t,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, const COUPLING_ALGO &cpl_algo, 
		   bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t
		if( fetch_t_hist_ != t && barrier_enabled )
			barrier(t_sampler.get_upper_bound(t));

		fetch_t_hist_ = t;

		std::vector<std::pair<std::pair<time_type,iterator_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
														std::numeric_limits<iterator_type>::lowest());

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
														std::numeric_limits<iterator_type>::lowest());
		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( start->first, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return cpl_algo.relaxation(std::make_pair(std::numeric_limits<time_type>::lowest(), static_cast<iterator_type>(t)), focus, t_sampler.filter(t, v));
	}

	/** \brief Fetch from the interface with coupling algorithms, blocking with barrier at time=t,it
	 */
	template<class SAMPLER, class TIME_SAMPLER, class COUPLING_ALGO, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr,const point_type& focus, const time_type t, const iterator_type it,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, const COUPLING_ALGO &cpl_algo, 
		   bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t,iteration=it
		if((fetch_t_hist_ != t || fetch_i_hist_ != it) && barrier_enabled)
			barrier(t_sampler.get_upper_bound(t),t_sampler.get_upper_bound(it));

		fetch_t_hist_ = t;
		fetch_i_hist_ = it;

		std::vector<std::pair<std::pair<time_type,iterator_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
														   t_sampler.get_lower_bound(it)-threshold(it));

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
														   t_sampler.get_upper_bound(it)+threshold(it));
		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( start->first, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return cpl_algo.relaxation(std::make_pair(t,it), focus, t_sampler.filter(std::make_pair(t,it), v));
	}

	/** \brief Fetch points currently stored in the interface, blocking with barrier at time=t
	*/
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<point_type>
	fetch_points( const std::string& attr, const time_type t,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t
		if( fetch_t_hist_ != t && barrier_enabled )
			barrier(t_sampler.get_upper_bound(t));

		fetch_t_hist_ = t;

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector <point_type> return_points;

		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
													       std::numeric_limits<iterator_type>::lowest());

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
														   std::numeric_limits<iterator_type>::lowest());
		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ){
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			const vec& ds = iter->second.template return_data<TYPE>();
			return_points.reserve(ds.size());
			for( size_t i=0; i<ds.size(); i++ ) {
				return_points.emplace_back(ds[i].first);
			}
		}

		return return_points;
	}

	/** \brief Fetch points currently stored in the interface, blocking with barrier at time=t,it
	*/
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<point_type>
	fetch_points( const std::string& attr, const time_type t, const iterator_type it,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t,iteration=it
		if( fetch_t_hist_ != t && fetch_i_hist_ != it && barrier_enabled)
			barrier(t_sampler.get_upper_bound(t),t_sampler.get_upper_bound(it));

		fetch_t_hist_ = t;
		fetch_i_hist_ = it;

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector <point_type> return_points;

		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
													       t_sampler.get_lower_bound(it)-threshold(it));

		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
														   t_sampler.get_upper_bound(it)+threshold(it));

		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ){
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			const vec& ds = iter->second.template return_data<TYPE>();
			return_points.reserve(ds.size());
			for( size_t i=0; i<ds.size(); i++ ) {
				return_points.emplace_back(ds[i].first);
			}
		}

		return return_points;
	}

	/** \brief Fetch values currently stored in the interface, blocking with barrier at time=t
	*/
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<TYPE>
	fetch_values( const std::string& attr, const time_type t,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t,iteration=it
		if( fetch_t_hist_ != t && barrier_enabled )
			barrier(t_sampler.get_upper_bound(t));

		fetch_t_hist_ = t;

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector<TYPE> return_values;

		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
													   	   std::numeric_limits<iterator_type>::lowest());
		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
													   	   std::numeric_limits<iterator_type>::lowest());

		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ){
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			const vec& ds = iter->second.template return_data<TYPE>();
			return_values.reserve(ds.size());
			for( size_t i=0; i<ds.size(); i++ ) {
				return_values.emplace_back(ds[i].second);
			}
		}

		return return_values;
	}

	/** \brief Fetch values currently stored in the interface, blocking with barrier at time=t,it
	*/
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<TYPE>
	fetch_values( const std::string& attr, const time_type t, const iterator_type it,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		// Only enter barrier on first fetch for time=t,iteration=it
		if( fetch_t_hist_ != t && fetch_i_hist_ != it && barrier_enabled)
			barrier(t_sampler.get_upper_bound(t),t_sampler.get_upper_bound(it));

		fetch_t_hist_ = t;
		fetch_i_hist_ = it;

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector<TYPE> return_values;

		std::pair<time_type,iterator_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),
													   	   t_sampler.get_lower_bound(it)-threshold(it));
		std::pair<time_type,iterator_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),
													   	   t_sampler.get_upper_bound(it)+threshold(it));

		auto end = log.upper_bound(curr_time_upper);

		if( log.size() == 1 ) end = log.end();

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ){
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			const vec& ds = iter->second.template return_data<TYPE>();
			return_values.reserve(ds.size());
			for( size_t i=0; i<ds.size(); i++ ) {
				return_values.emplace_back(ds[i].second);
			}
		}

		return return_values;
	}

	/** \brief Serializes pushed data and sends it to remote nodes
	* Serializes pushed data and sends it to remote nodes.
	* Returns the actual number of peers contacted
	*/
	int commit( time_type t, iterator_type it = std::numeric_limits<iterator_type>::lowest() ) {
		std::pair<time_type, iterator_type> time(t, it);

    	// Check Smart Send if announcement made
    	if ( !smart_send_set_ ) {
    		// Reset all peers to default of enabled
    		std::fill(peer_is_sending.begin(), peer_is_sending.end(), true);
    		update_smart_send(t);
    		smart_send_set_ = true;
    	}

		if( FIXEDPOINTS ) {
			// This only happens during the first commit
			if( push_buffer_pts.size() > 0 ) {
				comm->send( message::make("points",comm->local_rank(),std::move(push_buffer_pts)),peer_is_sending );
				initialized_pts_ = true;
				push_buffer_pts.clear();
			}

			// Reset counter for flat point structure
			fixedPointCount_ = 0;

			if( push_buffer_raw.size() > 0 ) {
				comm->send( message::make("rawdata",comm->local_rank(),time,std::move(push_buffer_raw)),peer_is_sending );
				push_buffer_raw.clear();
			}
		}
		else {
			if( push_buffer.size() > 0 ) {
				comm->send( message::make("data",time,std::move(push_buffer)),peer_is_sending );
				push_buffer.clear();
			}
		}

		comm->send( message::make("timestamp",comm->local_rank(),time),peer_is_sending );

		return std::count( peer_is_sending.begin(),peer_is_sending.end(),true );
	}

	/** \brief Updates Smart Send locality data
	* Creates a new comm rank mapping for Smart Send functionality
	*/
	void update_smart_send( time_type t ) {
		if( (((span_start < t) || almost_equal(span_start, t)) &&
	        ((t < span_timeout) || almost_equal(t, span_timeout))) ) {
			for( size_t i=0; i < peers.size(); i++ ) {
				// Check if peer is explicitly disabled
				if( peers[i].is_recv_disabled() ) {
					peer_is_sending[i] = false;
					continue;
				}

				// Perform geometric check against defined regions
				peer_is_sending[i] = peers[i].is_recving( t, current_span );
			}
	  }
	  else { // Ensure explicitly disabled peers are taken into account if outside Smart Send time bounds
		  for( size_t i=0; i < peers.size(); i++ ) {
			  if( peers[i].is_recv_disabled() ) peer_is_sending[i] = false;
		  }
	  }
	}

	/** \brief Sends a forecast of an upcoming time to remote nodes
	*/
	void forecast( time_type t, iterator_type it = std::numeric_limits<iterator_type>::lowest()) {
		std::pair<time_type,iterator_type> time(t,it);
		comm->send(message::make("forecast", comm->local_rank(), time));
	}

	/** \brief Tests whether data is available at time=t
	*/
	bool is_ready( const std::string& attr, time_type t ) const {
		using logitem_ref_t = typename decltype(log)::const_reference;
		return std::any_of(log.begin(), log.end(), [=](logitem_ref_t time_frame) {
			return time_frame.second.find(attr) != time_frame.second.end(); }) // return false for attributes that don't exist.
			&& std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
			return (p.is_send_disabled()) || (!p.is_sending(t, recv_span)) ||
				   ((((p.current_t() > t) || almost_equal(p.current_t(), t)) || (p.next_t() > t))); });
	}

	/** \brief Tests whether data is available at time=t,it
	*/
	bool is_ready( const std::string& attr, time_type t, iterator_type it ) const {
		using logitem_ref_t = typename decltype(log)::const_reference;
		return std::any_of(log.begin(), log.end(), [=](logitem_ref_t time_frame) {
			return time_frame.second.find(attr) != time_frame.second.end(); }) // return false for attributes that don't exist.
			&& std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
			return (p.is_send_disabled()) || (!p.is_sending(t, recv_span)) ||
				   ((((p.current_t() > t) || almost_equal(p.current_t(), t)) || (p.next_t() > t)) &&
					(((p.current_it() > it) || almost_equal(p.current_it(), it)) || (p.current_it() > it))); });
	}

	/** \brief Blocking barrier at time=t. Initiates receive from remote nodes.
	*/
	void barrier( time_type t ) {
		// barrier must be thread-safe because it is called in fetch()
		std::lock_guard<std::mutex> lock(mutex);

		auto start = std::chrono::system_clock::now();

		for(;;) {
			size_t peers_unblocked = 0;
			for( size_t p = 0; p < peers.size(); p++ ) {
				if( peers[p].is_send_disabled() ) { peers_unblocked++; continue; } // Rank disabled, immediate break
				if( !peers[p].is_sending(t, recv_span) ) { peers_unblocked++; continue; } // Rank disabled due to Smart Send geometry check
				if( (peers[p].current_t() > t || almost_equal(peers[p].current_t(), t)) || peers[p].next_t() > t ) { // Final time check
					peers_unblocked++;
					continue;
				}
			}
			// All peers unblocked, break loop
			if( peers_unblocked == peers.size() )
				break;
			else // Acquire messages
				acquire();
		}

		if( !QUIET ) {
			if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
				std::cout << "MUI Warning [uniface.h]: Communication barrier spent over 5 seconds" << std::endl;
			}
		}
	}

	/** \brief Blocking barrier at time=t,it. Initiates receive from remote nodes.
	*/
	void barrier( time_type t, iterator_type it ) {
		// barrier must be thread-safe because it is called in fetch()
		std::lock_guard<std::mutex> lock(mutex);

		auto start = std::chrono::system_clock::now();

		for(;;) {
			size_t peers_unblocked = 0;
			for( size_t p = 0; p < peers.size(); p++ ) {
				if( peers[p].is_send_disabled() ) { peers_unblocked++; continue; } // Rank disabled, immediate break
				if( !peers[p].is_sending(t, recv_span) ) { peers_unblocked++; continue; } // Rank disabled due to Smart Send geometry check
				if( ((peers[p].current_t() > t || almost_equal(peers[p].current_t(), t)) || peers[p].next_t() > t) && // Final time check
					((peers[p].current_it() > it || almost_equal(peers[p].current_it(), it)) || peers[p].next_it() > it) ) {
					peers_unblocked++;
					continue;
				}
			}
			// All peers unblocked, break loop
			if( peers_unblocked == peers.size() )
				break;
			else // Acquire messages
				acquire();
		}

		if( !QUIET ) {
			if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
				std::cout << "MUI Warning [uniface.h]: Communication barrier spent over 5 seconds" << std::endl;
			}
		}
	}

	/** \brief Blocking barrier for Smart Send send values. Initiates receive from remote nodes.
    */
	void barrier_ss_send( ) {
		// barrier must be thread-safe because it is called in fetch()
		std::lock_guard<std::mutex> lock(mutex);

		auto start = std::chrono::system_clock::now();

		for(;;) {    // barrier must be thread-safe because it is called in fetch()
			if( std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
				return (p.ss_send_status()); }) ) break;
			acquire(); // To avoid infinite-loop when synchronous communication
		}

		for(size_t i=0; i<peers.size(); i++) {
			peers[i].set_ss_send_status(false);
		}

		if( !QUIET ) {
			if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
				std::cout << "MUI Warning [uniface.h]: Smart Send communication barrier spent over 5 seconds" << std::endl;
			}
		}
	}

	/** \brief Blocking barrier for Smart Send receive values. Initiates receive from remote nodes.
    */
	void barrier_ss_recv( ) {
		// barrier must be thread-safe because it is called in fetch()
		std::lock_guard<std::mutex> lock(mutex);

		auto start = std::chrono::system_clock::now();

		for(;;) {    // barrier must be thread-safe because it is called in fetch()
			if( std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
				return (p.ss_recv_status()); }) ) break;
			acquire(); // To avoid infinite-loop when synchronous communication
		}

		for(size_t i=0; i<peers.size(); i++) {
			peers[i].set_ss_recv_status(false);
		}

		if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
			if( !QUIET )
				std::cout << "MUI Warning [uniface.h]: Smart Send communication barrier spent over 5 seconds" << std::endl;
		}
	}

	/** \brief Announces to all remote nodes using non-blocking peer-to-peer approach "I'll send this span"
	*/
	void announce_send_span( time_type start, time_type timeout, span_t s, bool synchronised = false) {
		span_start = start;
		span_timeout = timeout;
		current_span.swap(s);
		comm->send(message::make("sendingSpan", comm->local_rank(), start, timeout, std::move(current_span)));
		if( synchronised ) barrier_ss_send();
		smart_send_set_ = false;
	}

	/** \brief Announces to all remote nodes "I'm disabled for send"
	*/
	void announce_send_disable( bool synchronised = false ) {
		comm->send(message::make("sendingDisable", comm->local_rank()));
		if( synchronised ) barrier_ss_send();
	}

	/** \brief Announces to all remote nodes using non-blocking peer-to-peer approach "I'm receiving this span"
	*/
	void announce_recv_span( time_type start, time_type timeout, span_t s, bool synchronised = false ) {
		recv_start = start;
		recv_timeout = timeout;
		recv_span.swap(s);
		comm->send(message::make("receivingSpan", comm->local_rank(), start, timeout, std::move(recv_span)));
		if( synchronised ) barrier_ss_recv();
		smart_send_set_ = false;
	}

	/** \brief Announces to all remote nodes "I'm disabled for receive"
	*/
	void announce_recv_disable( bool synchronised = false ) {
		comm->send(message::make("receivingDisable", comm->local_rank()));
		if( synchronised ) barrier_ss_recv();
	}

	/** \brief Removes log between (-inf, @last]
	*/
	void forget( time_type last, bool reset_log = false ) {
		std::pair<time_type,iterator_type> upper_limit(last+threshold(last),
													   std::numeric_limits<iterator_type>::lowest());

		log.erase(log.begin(), log.upper_bound(upper_limit));

		if( reset_log ) {
			std::pair<time_type,iterator_type> curr_time(std::numeric_limits<time_type>::lowest(),
														 std::numeric_limits<iterator_type>::lowest());

			if( !log.empty() ) curr_time = log.rbegin()->first;

			for( size_t i=0; i < peers.size(); i++ ) {
				peers[i].set_current_t(curr_time.first);
				peers[i].set_current_sub(curr_time.second);
			}
		}

		fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
	}

	/** \brief Removes log between ([-inf,-inf], [@last.first,@last.second]]
	*/
	void forget( std::pair<time_type,iterator_type> last, bool reset_log = false ) {
		std::pair<time_type,iterator_type> upper_limit(last.first+threshold(last.first),
												   last.second+threshold(last.second));

		log.erase(log.begin(), log.upper_bound(upper_limit));

		if( reset_log ) {
			std::pair<time_type,iterator_type> curr_time(std::numeric_limits<time_type>::lowest(),
														 std::numeric_limits<iterator_type>::lowest());

			if( !log.empty() ) curr_time = log.rbegin()->first;

			for( size_t i=0; i < peers.size(); i++ ) {
				peers[i].set_current_t(curr_time.first);
				peers[i].set_current_sub(curr_time.second);
			}
		}

		fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
		fetch_i_hist_ = std::numeric_limits<iterator_type>::lowest();
	}

	/** \brief Removes log between [@first, @last]
	*/
	void forget( time_type first, time_type last, bool reset_log = false ) {
		std::pair<time_type,iterator_type> lower_limit(first-threshold(first),
													   std::numeric_limits<iterator_type>::lowest());
		std::pair<time_type,iterator_type> upper_limit(last+threshold(last),
				         	 	 	 	 	 	 	   std::numeric_limits<iterator_type>::lowest());

		log.erase(log.lower_bound(lower_limit), log.upper_bound(upper_limit));

		if( reset_log ) {
			std::pair<time_type,iterator_type> curr_time(std::numeric_limits<time_type>::lowest(),
														 std::numeric_limits<iterator_type>::lowest());

			if( !log.empty() ) curr_time = log.rbegin()->first;

			for( size_t i=0; i < peers.size(); i++ ) {
				peers[i].set_current_t(curr_time.first);
				peers[i].set_current_sub(curr_time.second);
			}
		}

		fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
	}

	/** \brief Removes log between [[@first.first,@first.second], [@last.first,@last.second]]
	*/
	void forget( std::pair<time_type,iterator_type> first, std::pair<time_type,iterator_type> last, bool reset_log = false ) {
		std::pair<time_type,iterator_type> lower_limit(first.first-threshold(first.first),
													   first.second-threshold(first.second));
		std::pair<time_type,iterator_type> upper_limit(last.first+threshold(last.first),
													   last.second+threshold(last.second));

		log.erase(log.lower_bound(lower_limit), log.upper_bound(upper_limit));

		if( reset_log ) {
			std::pair<time_type,iterator_type> curr_time(std::numeric_limits<time_type>::lowest(),
														 std::numeric_limits<iterator_type>::lowest());

			if( !log.empty() ) curr_time = log.rbegin()->first;

			for( size_t i=0; i<peers.size(); i++ ) {
				peers[i].set_current_t(curr_time.first);
				peers[i].set_current_sub(curr_time.second);
			}
		}

		fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
		fetch_i_hist_ = std::numeric_limits<iterator_type>::lowest();
	}

	/** \brief Removes log between (-inf, current-@length] automatically.
	*/
	void set_memory( time_type length ) {
		memory_length = length;

		fetch_t_hist_ = std::numeric_limits<time_type>::lowest();
		fetch_i_hist_ = std::numeric_limits<iterator_type>::lowest();
	}
	
	/** \brief Returns the URI host (domain) for the created uniface
	*/
	std::string uri_host() {
		return comm->uri_host();
	}

	/** \brief Returns the URI path (name) for the created uniface
	*/
	std::string uri_path() {
		return comm->uri_path();
	}

	/** \brief Returns the URI protocol for the created uniface
	*/
	std::string uri_protocol() {
		return comm->uri_protocol();
	}

private:
	/** \brief Triggers communication
	*/
	void acquire() {
		message m = comm->recv();
		if( m.has_id() ) readers[m.id()](m);
	}

	/** \brief Handles "timestamp" messages
	*/
	void on_recv_confirm( int32_t sender, std::pair<time_type,iterator_type> timestamp ) {
		peers[sender].set_current_t(timestamp.first);
		peers[sender].set_current_sub(timestamp.second);
	}

	/** \brief Handles "forecast" messages
	*/
	void on_recv_forecast( int32_t sender, std::pair<time_type,iterator_type> timestamp ) {
		peers[sender].set_next_t(timestamp.first);
		peers[sender].set_next_sub(timestamp.second);
	}

	/** \brief Handles "data" messages
	*/
	void on_recv_data( std::pair<time_type,iterator_type> timestamp, frame_type frame ) {
		auto itr = log.find(timestamp);

		if( itr == log.end() )
			std::tie(itr,std::ignore) = log.insert(std::make_pair(timestamp,bin_frame_type()));

		auto& cur = itr->second;

		for( auto& p: frame ){
			auto pstr = cur.find(p.first);
			if( pstr == cur.end() ) cur.insert(std::make_pair(std::move(p.first),spatial_t(std::move(p.second))));
			else pstr->second.insert(p.second);
		}

		log.erase(log.begin(), log.upper_bound({timestamp.first-memory_length, timestamp.second}));
	}

	/** \brief Handles "data" messages
	*/
	void on_recv_rawdata( int32_t sender, std::pair<time_type,iterator_type> timestamp, frame_raw_type frame ) {
		on_recv_data( timestamp, associate( sender, frame ) );
	}

	/** \brief Handles "receivingSpan" messages
	*/
	void on_recv_span( int32_t sender, time_type start, time_type timeout, span_t s ) {
		peers[sender].set_recving(start,timeout,std::move(s));
		peers[sender].set_ss_recv_status(true);
	}

	/** \brief Handles "sendingSpan" messages
	*/
	void on_send_span( int32_t sender, time_type start, time_type timeout, span_t s ) {
		peers[sender].set_sending(start,timeout,std::move(s));
		peers[sender].set_ss_send_status(true);
	}

	/** \brief Handles "sendingDisable" messages
	*/
	void on_recv_disable( int32_t sender ) {
		peers[sender].set_recv_disable();
		peers[sender].set_ss_recv_status(true);
		peer_is_sending[sender] = false;
	}

	/** \brief Handles "receivingDisable" messages
	*/
	void on_send_disable( int32_t sender ) {
		peers[sender].set_send_disable();
		peers[sender].set_ss_send_status(true);
	}

	/** \brief Handles "points" messages
	*/
	void on_recv_points( int32_t sender, std::vector<point_type> points ) {
		peers[sender].set_pts(points);
	}

	/** \brief Handles "assignedVals" messages
	*/
	void on_recv_assignedVals( std::string attr, storage_single_t data ) {
		typename std::unordered_map<std::string, storage_single_t >::iterator it = assigned_values.find(attr);
		if (it != assigned_values.end())
			it->second = data;
		else 
			assigned_values.insert( std::pair<std::string, storage_single_t>( attr, data ) );
	}

	/** \brief Associates raw data and stored point data together
	*/
	inline frame_type associate( int32_t sender, frame_raw_type& frame ) {
		frame_type buf;
		const auto& pts = peers[sender].pts();

		for( auto& p: frame ) {
			const auto& data = storage_cast<const std::vector<std::pair<size_t,REAL> >&>(p.second);

			buf.insert(std::make_pair(p.first, storage_t(std::vector<std::pair<point_type,REAL> >())));
			std::vector<std::pair<point_type,REAL> >& data_store = storage_cast<std::vector<std::pair<point_type,REAL> >&>(buf[p.first]);

		    data_store.resize(data.size());

			for( size_t i=0; i<data.size(); i++ ) {
				data_store[i].first = pts[data[i].first];
				data_store[i].second = data[i].second;
			}
		}

		return buf;
	}
};

}

#endif // _UNIFACE_H
