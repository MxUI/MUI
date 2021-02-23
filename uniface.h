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

#include "util.h"
#include "comm.h"
#include "comm_factory.h"
#include "config.h"
#include "dynstorage.h"
#include "spatial_storage.h"
#include "lib_dispatcher.h"
#include "message.h"
#include "reader_variable.h"
#include "stream_vector.h"
#include "stream_unordered.h"
#include "stream_string.h"
#include "bin.h"
#include "stream.h"

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
	using time_type  = typename CONFIG::time_type;
	using data_types = typename CONFIG::data_types;
	using span_t = geometry::any_shape<CONFIG>;
private:
	// meta functions to split tuple and add vector<pair<point_type,_1> >
	template<typename T> struct add_vp_  { using type = std::vector<std::pair<point_type,T> >; };
	template<typename T> struct add_v_   { using type = std::vector<T>; };
	template<typename T> struct def_storage_;
	template<typename... TYPES> struct def_storage_<type_list<TYPES...> >{
		using type = storage<typename add_vp_<TYPES>::type...>;
	};

	template<typename T> struct def_storage_raw_;
	template<typename... TYPES> struct def_storage_raw_<type_list<TYPES...> >{
		using type = storage<typename add_v_<TYPES>::type...>;
	};

	template<typename T> struct def_storage_single_;
	template<typename... TYPES> struct def_storage_single_<type_list<TYPES...> >{
		using type = storage<TYPES...>;
	};

	// internal typedefinitions
	using storage_t = typename def_storage_<data_types>::type;
	using spatial_t = spatial_storage<bin_t<CONFIG>,storage_t,CONFIG>;
	using frame_type = std::unordered_map<std::string, storage_t>;
	using bin_frame_type = std::unordered_map<std::string, spatial_t>;

	using storage_raw_t = typename def_storage_raw_<data_types>::type;
	using frame_raw_type = std::unordered_map<std::string, storage_raw_t>;

	using storage_single_t = typename def_storage_single_<data_types>::type;

	struct peer_state {
		peer_state() : disable_send(false), disable_recv(false) {}

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

		void set_pts(std::vector<point_type> pts) {
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

		time_type current_t() const { return latest_timestamp; }
		time_type current_sub() const { return latest_subiter; }
		time_type next_t() const { return next_timestamp; }
		time_type next_sub() const { return next_subiter; }
		void set_current_t( time_type t ) { latest_timestamp = t; }
		void set_current_sub( time_type t ) { latest_subiter = t; }
		void set_next_t( time_type t ) { next_timestamp = t; }
		void set_next_sub( time_type t ) { next_subiter = t; }
	private:
		bool scan_spans_(time_type t, const span_t& s, const spans_type& spans ) const {
			bool prefetched = false;
			auto end = spans.end();
			if(spans.size() > 1)
				end = spans.lower_bound(std::make_pair(t,t));

			for( auto itr = spans.begin(); itr != end; ++itr ) {
			    if( (t < itr->first.second) || almost_equal(t, itr->first.second) ) {
					prefetched = true;
					if( collide(s,itr->second) ) return true;
				}
			}
			// if prefetched at t, but no overlap region, then return false;
			// otherwise return true;
			return !prefetched;
		}
		time_type latest_timestamp = std::numeric_limits<time_type>::lowest();
		time_type latest_subiter = std::numeric_limits<time_type>::lowest();
		time_type next_timestamp = std::numeric_limits<time_type>::lowest();
		time_type next_subiter = std::numeric_limits<time_type>::lowest();
		spans_type recving_spans;
		spans_type sending_spans;
		std::vector<point_type> pts_;
		std::unordered_map<std::string, storage_single_t> assigned_vals_;
		bool disable_send;
		bool disable_recv;
	};

private: // data members
	std::unique_ptr<communicator> comm;
	dispatcher<message::id_type, std::function<void(message)> > readers;

	std::map<std::pair<time_type, time_type>, bin_frame_type> log;

	frame_type push_buffer;
	frame_raw_type push_buffer_raw;
	std::vector<point_type> push_buffer_pts;

	std::unordered_map<std::string, storage_single_t > assigned_values;

	std::vector<peer_state> peers;
	time_type span_start = std::numeric_limits<time_type>::lowest();
	time_type span_timeout = std::numeric_limits<time_type>::lowest();
	span_t current_span;
	time_type recv_start = std::numeric_limits<time_type>::lowest();
	time_type recv_timeout = std::numeric_limits<time_type>::lowest();
	span_t recv_span;
	time_type memory_length = std::numeric_limits<time_type>::max();
	std::mutex mutex;
	bool initialized_pts_;

public:
	uniface( const char URI[] ) : uniface( comm_factory::create_comm(URI) ) {}
	uniface( std::string const &URI ) : uniface( comm_factory::create_comm(URI.c_str()) ) {}
	uniface( communicator* comm_ ) : comm(comm_), initialized_pts_(false) {
		using namespace std::placeholders;

		peers.resize(comm->remote_size());

		readers.link("timestamp", reader_variables<int32_t, std::pair<time_type,time_type> >(
					 std::bind(&uniface::on_recv_confirm, this, _1, _2)));
		readers.link("forecast", reader_variables<int32_t, std::pair<time_type,time_type>>(
					 std::bind(&uniface::on_recv_forecast, this, _1, _2)));
		readers.link("data", reader_variables<std::pair<time_type,time_type>, frame_type>(
					 std::bind(&uniface::on_recv_data, this, _1, _2)));
		readers.link("rawdata", reader_variables<int32_t, std::pair<time_type,time_type>, frame_raw_type>(
					 std::bind(&uniface::on_recv_rawdata, this, _1, _2, _3)));
		readers.link("points", reader_variables<int32_t, std::vector<point_type>>(
					 std::bind(&uniface::on_recv_points, this, _1, _2)));
		readers.link("assignedVals", reader_variables<std::string, storage_single_t>(
					 std::bind(&uniface::on_recv_assignedVals, this, _1, _2)));
		readers.link("receivingSpan", reader_variables<int32_t, time_type, time_type, span_t>(
		             std::bind(&uniface::on_recv_span, this, _1, _2, _3, _4)));
		readers.link("sendingSpan", reader_variables<int32_t, time_type, time_type, span_t>(
		             std::bind(&uniface::on_send_span, this, _1, _2, _3, _4)));
		readers.link("receivingDisable", reader_variables<int32_t>(
					 std::bind(&uniface::on_send_disable, this, _1)));
		readers.link("sendingDisable", reader_variables<int32_t>(
					 std::bind(&uniface::on_recv_disable, this, _1)));
	}

	uniface( const uniface& ) = delete;
	uniface& operator=( const uniface& ) = delete;
    
    /** \brief Push data with tag "attr" to buffer
     * Push data with tag "attr" to buffer. If using CONFIG::FIXEDPOINTS=true,
     * data must be pushed in the same order that the points were previously pushed.
     */
    template<typename TYPE>
	void push( const std::string& attr, const point_type loc, const TYPE value ) {
		if( FIXEDPOINTS ) {
			if( !initialized_pts_ ) push_buffer_pts.emplace_back( loc );
			storage_raw_t& n = push_buffer_raw[attr];
			if( !n ) n = storage_raw_t(std::vector<TYPE>());
			storage_cast<std::vector<TYPE>&>(n).emplace_back( value );
		} 
		else {
			storage_t& n = push_buffer[attr];
			if( !n ) n = storage_t(std::vector<std::pair<point_type,TYPE> >());
			storage_cast<std::vector<std::pair<point_type,TYPE> >&>(n).emplace_back( loc, value );
		}
	}

	/** \brief Push the value \c value to the parameter \c attr
	  * Useful if, for example, you wish to pass a parameter
	  * rather than a field without an associated timestamp
	  */
	template<typename TYPE>
	void push( const std::string& attr, const TYPE value ) {
		comm->send(message::make("assignedVals", attr, storage_single_t(TYPE(value))));
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
           const SAMPLER& sampler, const TIME_SAMPLER &t_sampler) {
        // Arrays must have ndim = d; can be non-writeable
        point_type p = 0;
        auto points_arr = points.template unchecked<2>();
        py::array_t<typename SAMPLER::OTYPE,py::array::c_style> values(points_arr.shape(0));
        auto values_arr = values.template mutable_unchecked<1>();
        for (ssize_t i = 0; i < points_arr.shape(0); i++) {
            for (ssize_t j = 0; j < points_arr.shape(1); j++)
                p[j] = points_arr(i,j);
            values_arr(i)  = fetch(attr, p, t, sampler, t_sampler);
        }
        return values;
    }

    template<typename TYPE>
    py::array_t<REAL, py::array::c_style>
    fetch_points_np( const std::string& attr, const time_type t ) {
        std::vector<point_type> points = fetch_points<TYPE>(attr, t);
        size_t n = points.size();
        py::array_t<REAL, py::array::c_style> points_np({n, n});
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
	fetch( const std::string& attr,const point_type focus, const time_type t,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true,
		   ADDITIONAL && ... additional ) {
		if(barrier_enabled)
			barrier(t_sampler.get_upper_bound(t));

		std::vector<std::pair<std::pair<time_type,time_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),std::numeric_limits<time_type>::lowest());
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),std::numeric_limits<time_type>::lowest());

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			std::pair<time_type,time_type> time(start->first.first,std::numeric_limits<time_type>::lowest());
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( time, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return t_sampler.filter(t, v);
	}

	/** \brief Fetch from the interface, blocking with barrier at time=t1,t2
	  */
	template<class SAMPLER, class TIME_SAMPLER, typename ... ADDITIONAL>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr,const point_type focus, const time_type t1, const time_type t2,
		   SAMPLER& sampler, const TIME_SAMPLER &t_sampler, bool barrier_enabled = true,
		   ADDITIONAL && ... additional ) {
		if(barrier_enabled)
			barrier(t_sampler.get_upper_bound(t1),t_sampler.get_upper_bound(t2));

		std::vector<std::pair<std::pair<time_type,time_type>,typename SAMPLER::OTYPE> > v;
		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t1)-threshold(t1),t_sampler.get_lower_bound(t2)-threshold(t2));
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t1)+threshold(t1),t_sampler.get_upper_bound(t2)+threshold(t2));

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

		for( auto start = log.lower_bound(curr_time_lower); start != end; ++start ) {
			std::pair<time_type,time_type> time(start->first);
			const auto& iter = start->second.find(attr);
			if( iter == start->second.end() ) continue;
			v.emplace_back( time, iter->second.build_and_query_ts( focus, sampler, additional... ) );
		}

		return t_sampler.filter(std::make_pair(t1,t2), v);
	}

	/** \brief Fetch points currently stored in the interface, blocking with barrier at time=t
	 */
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<point_type>
	fetch_points( const std::string& attr, const time_type t,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
	  if(barrier_enabled)
	    barrier(t_sampler.get_upper_bound(t));

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector <point_type> return_points;

		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),std::numeric_limits<time_type>::lowest());
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),std::numeric_limits<time_type>::lowest());

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

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

	/** \brief Fetch points currently stored in the interface, blocking with barrier at time=t1,t2
	 */
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<point_type>
	fetch_points( const std::string& attr, const time_type t1, const time_type t2,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		if(barrier_enabled)
			barrier(t_sampler.get_upper_bound(t1),t_sampler.get_upper_bound(t2));

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector <point_type> return_points;

		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t1)-threshold(t1),t_sampler.get_lower_bound(t2)-threshold(t2));
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t1)+threshold(t1),t_sampler.get_upper_bound(t2)+threshold(t2));

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

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
		if(barrier_enabled)
			barrier(t_sampler.get_upper_bound(t));

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector<TYPE> return_values;

		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t)-threshold(t),std::numeric_limits<time_type>::lowest());
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t)+threshold(t),std::numeric_limits<time_type>::lowest());

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

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

	/** \brief Fetch values currently stored in the interface, blocking with barrier at time=t1,t2
	 */
	template<typename TYPE, class TIME_SAMPLER, typename ... ADDITIONAL>
	std::vector<TYPE>
	fetch_values( const std::string& attr, const time_type t1, const time_type t2,
				  const TIME_SAMPLER &t_sampler, bool barrier_enabled = true, ADDITIONAL && ... additional ) {
		if(barrier_enabled)
			barrier(t_sampler.get_upper_bound(t1),t_sampler.get_upper_bound(t2));

		using vec = std::vector<std::pair<point_type,TYPE> >;
		std::vector<TYPE> return_values;

		std::pair<time_type,time_type> curr_time_lower(t_sampler.get_lower_bound(t1)-threshold(t1),t_sampler.get_lower_bound(t2)-threshold(t2));
		std::pair<time_type,time_type> curr_time_upper(t_sampler.get_upper_bound(t1)+threshold(t1),t_sampler.get_upper_bound(t2)+threshold(t2));

		auto end = log.end();
		if(log.size() > 1)
			end = log.upper_bound(curr_time_upper);

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
	int commit( time_type t1, time_type t2 = std::numeric_limits<time_type>::lowest() ) {
	    std::vector<bool> is_sending(comm->remote_size(), true);
	    std::vector<bool> is_enabled(comm->remote_size(), true);
	    std::pair<time_type,time_type> time(t1,t2);

	    // Check if peer set to disabled (not linked to time span)
	    for( std::size_t i=0; i<peers.size(); ++i ) {
	    	if(peers[i].is_recv_disabled())
	    		is_enabled[i] = false;
		}

	    // Check for smart send based on t1
	    if( (((span_start < t1) || almost_equal(span_start, t1)) &&
	    	((t1 < span_timeout) || almost_equal(t1, span_timeout))) ) {
			for( std::size_t i=0; i<peers.size(); ++i ) {
				if(!is_enabled[i]) // Peer is completely disabled
					is_sending[i] = false;
				else // Check peer using typical smart send procedure
					is_sending[i] = peers[i].is_recving( t1, current_span );
			}
		}

		if( FIXEDPOINTS ) {
			if( push_buffer_pts.size() > 0 ) {
				comm->send(message::make("points",comm->local_rank(),std::move(push_buffer_pts)), is_sending);
				initialized_pts_ = true;
			}
			comm->send(message::make("rawdata",comm->local_rank(),time,std::move(push_buffer_raw)), is_sending);
			push_buffer_raw.clear();
			push_buffer_pts.clear();
		}
		else {
			comm->send(message::make("data",time,std::move(push_buffer)), is_sending);
			push_buffer.clear();
		}

		comm->send(message::make("timestamp",comm->local_rank(),time), is_sending);

		return std::count( is_sending.begin(), is_sending.end(), true );
	}

	/** \brief Sends a forecast of an upcoming time to remote nodes
	  */
	void forecast( time_type t1, time_type t2 = std::numeric_limits<time_type>::lowest()) {
		std::pair<time_type,time_type> time(t1,t2);
		comm->send(message::make("forecast", comm->local_rank(), time));
	}

	/** \brief Tests whether data is available at time t1
	  */
	bool is_ready( const std::string& attr, time_type t1 ) const {
		using logitem_ref_t = typename decltype(log)::const_reference;
		return std::any_of(log.begin(), log.end(), [=](logitem_ref_t time_frame) {
			return time_frame.second.find(attr) != time_frame.second.end(); }) // return false for attributes that don't exist.
			&& std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
			return (p.is_send_disabled()) || (!p.is_sending(t1, recv_span)) ||
				   ((((p.current_t() > t1) || almost_equal(p.current_t(), t1)) || (p.next_t() > t1))); });
	}

	/** \brief Tests whether data is available at time t1 (or t1,t2)
	  */
	bool is_ready( const std::string& attr, time_type t1, time_type t2 ) const {
		using logitem_ref_t = typename decltype(log)::const_reference;
		return std::any_of(log.begin(), log.end(), [=](logitem_ref_t time_frame) {
			return time_frame.second.find(attr) != time_frame.second.end(); }) // return false for attributes that don't exist.
			&& std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
			return (p.is_send_disabled()) || (!p.is_sending(t1, recv_span)) ||
				   ((((p.current_t() > t1) || almost_equal(p.current_t(), t1)) || (p.next_t() > t1)) &&
					(((p.current_sub() > t2) || almost_equal(p.current_sub(), t2)) || (p.current_sub() > t2))); });
	}

	/** \brief Blocking barrier at t1. Initiates receive from remote nodes.
	  */
	void barrier( time_type t1 ) {
		auto start = std::chrono::system_clock::now();
		for(;;) {    // barrier must be thread-safe because it is called in fetch()
			std::lock_guard<std::mutex> lock(mutex);
			if( std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
				return (p.is_send_disabled()) || (!p.is_sending(t1, recv_span)) ||
					   ((((p.current_t() > t1) || almost_equal(p.current_t(), t1)) || (p.next_t() > t1))); }) ) break;
			acquire(); // To avoid infinite-loop when synchronous communication
		}
		if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
			if( !QUIET )
				std::cout << "MUI Warning [uniface.h]: Communication barrier spends over 5 seconds" << std::endl;
		}
	}

	/** \brief Blocking barrier at t1 (or t1,t2). Initiates receive from remote nodes.
	  */
	void barrier( time_type t1, time_type t2 ) {
		auto start = std::chrono::system_clock::now();
		for(;;) {    // barrier must be thread-safe because it is called in fetch()
			std::lock_guard<std::mutex> lock(mutex);
			if( std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
				return (p.is_send_disabled()) || (!p.is_sending(t1, recv_span)) ||
					   ((((p.current_t() > t1) || almost_equal(p.current_t(), t1)) || (p.next_t() > t1)) &&
					    (((p.current_sub() > t2) || almost_equal(p.current_sub(), t2)) || (p.next_sub() > t2))); }) ) break;
			acquire(); // To avoid infinite-loop when synchronous communication
		}
		if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
			if( !QUIET )
				std::cout << "MUI Warning [uniface.h]: Communication barrier spends over 5 seconds" << std::endl;
		}
	}

	/** \brief Announces to all remote nodes "I'll send this span"
	  */
	void announce_send_span( time_type start, time_type timeout, span_t s ) {
	    comm->send(message::make("sendingSpan", comm->local_rank(), start, timeout, std::move(s)));
	    span_start = start;
	    span_timeout = timeout;
	    current_span.swap(s);
	}

	/** \brief Announces to all remote nodes "I'm disabled for send"
	  */
	void announce_send_disable() {
		comm->send(message::make("sendingDisable", comm->local_rank()));
	}

	/** \brief Announces to all remote nodes "I'm receiving this span"
	  */
	void announce_recv_span( time_type start, time_type timeout, span_t s ) {
		comm->send(message::make("receivingSpan", comm->local_rank(), start, timeout, std::move(s)));
		recv_start = start;
		recv_timeout = timeout;
		recv_span.swap(s);
	}

	/** \brief Announces to all remote nodes "I'm disabled for receive"
	  */
	void announce_recv_disable() {
		comm->send(message::make("receivingDisable", comm->local_rank()));
	}

	/** \brief Removes log between (-inf, @last]
	  */
	void forget( time_type last, bool reset_log = false ) {
		std::pair<time_type,time_type> upper_limit(last+threshold(last),std::numeric_limits<time_type>::lowest());
		log.erase(log.begin(), log.upper_bound(upper_limit));

		if(reset_log) {
			std::pair<time_type,time_type> curr_time(std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest());
			if(!log.empty()) curr_time = log.rbegin()->first;
			for(size_t i=0; i<peers.size(); i++) {
				peers.at(i).set_current_t(curr_time.first);
				peers.at(i).set_current_sub(curr_time.second);
			}
		}
	}

	/** \brief Removes log between ([-inf,-inf], [@last.first,@last.second]]
	  */
	void forget( std::pair<time_type, time_type> last, bool reset_log = false ) {
		std::pair<time_type,time_type> upper_limit(last.first+threshold(last.first),last.second+threshold(last.second));
		log.erase(log.begin(), log.upper_bound(upper_limit));

		if(reset_log) {
			std::pair<time_type,time_type> curr_time(std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest());
			if(!log.empty()) curr_time = log.rbegin()->first;
			for(size_t i=0; i<peers.size(); i++) {
				peers.at(i).set_current_t(curr_time.first);
				peers.at(i).set_current_sub(curr_time.second);
			}
		}
	}

	/** \brief Removes log between [@first, @last]
	  */
	void forget( time_type first, time_type last, bool reset_log = false ) {
		std::pair<time_type,time_type> lower_limit(first-threshold(first),std::numeric_limits<time_type>::lowest());
		std::pair<time_type,time_type> upper_limit(last+threshold(last),std::numeric_limits<time_type>::lowest());
		log.erase(log.lower_bound(lower_limit), log.upper_bound(upper_limit));

		if(reset_log) {
			std::pair<time_type,time_type> curr_time(std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest());
			if(!log.empty()) curr_time = log.rbegin()->first;
			for(size_t i=0; i<peers.size(); i++) {
				peers.at(i).set_current_t(curr_time.first);
				peers.at(i).set_current_sub(curr_time.second);
			}
		}
	}

	/** \brief Removes log between [[@first.first,@first.second], [@last.first,@last.second]]
	  */
	void forget( std::pair<time_type, time_type> first, std::pair<time_type, time_type> last, bool reset_log = false ) {
		std::pair<time_type,time_type> lower_limit(first.first-threshold(first.first),first.second-threshold(first.second));
		std::pair<time_type,time_type> upper_limit(last.first+threshold(last.first),last.second+threshold(last.second));
		log.erase(log.lower_bound(lower_limit), log.upper_bound(upper_limit));

		if(reset_log) {
			std::pair<time_type,time_type> curr_time(std::numeric_limits<time_type>::lowest(),std::numeric_limits<time_type>::lowest());
			if(!log.empty()) curr_time = log.rbegin()->first;
			for(size_t i=0; i<peers.size(); i++) {
				peers.at(i).set_current_t(curr_time.first);
				peers.at(i).set_current_sub(curr_time.second);
			}
		}
	}

	/** \brief Removes log between (-inf, current-@length] automatically.
	  */
	void set_memory( time_type length ) {
		memory_length = length;
	}
	
private:
	/** \brief Triggers communication
	  */
	void acquire() {
		message m = comm->recv();
		if( m.has_id() ){
			readers[m.id()](m);
		}
	}

	/** \brief Handles "timestamp" messages
	  */
	void on_recv_confirm( int32_t sender, std::pair<time_type,time_type> timestamp ) {
		peers.at(sender).set_current_t(timestamp.first);
		peers.at(sender).set_current_sub(timestamp.second);
	}

	/** \brief Handles "forecast" messages
	  */
	void on_recv_forecast( int32_t sender, std::pair<time_type,time_type> timestamp ) {
		peers.at(sender).set_next_t(timestamp.first);
		peers.at(sender).set_next_sub(timestamp.second);
	}

	/** \brief Handles "data" messages
	  */
	void on_recv_data( std::pair<time_type,time_type> timestamp, frame_type frame ) {
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
	void on_recv_rawdata( int32_t sender, std::pair<time_type,time_type> timestamp, frame_raw_type frame ) {
		frame_type buf = associate( sender, frame );
		on_recv_data( timestamp, buf );
	}

	/** \brief Handles "receivingSpan" messages
	  */
	void on_recv_span( int32_t sender, time_type start, time_type timeout, span_t s ) {
		peers.at(sender).set_recving(start,timeout,std::move(s));
	}

	/** \brief Handles "sendingSpan" messages
	  */
	void on_send_span( int32_t sender, time_type start, time_type timeout, span_t s ) {
		peers.at(sender).set_sending(start,timeout,std::move(s));
	}

	/** \brief Handles "sendingDisable" messages
	  */
	void on_recv_disable( int32_t sender ) {
		peers.at(sender).set_recv_disable();
	}

	/** \brief Handles "receivingDisable" messages
	  */
	void on_send_disable( int32_t sender ) {
		peers.at(sender).set_send_disable();
	}

	/** \brief Handles "points" messages
	  */
	void on_recv_points( int32_t sender, std::vector<point_type> points ) {
		peers.at(sender).set_pts(points);
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
    
	/** \brief Associates raw data and point data together after both received
	  */
	frame_type associate( int32_t sender, frame_raw_type& frame ) {
		frame_type buf;
		for( auto& p: frame ){
			storage_t& n = buf[p.first];
			if( !n ) n = storage_t(std::vector<std::pair<point_type,REAL> >());

			const std::vector<REAL>& data = storage_cast<const std::vector<REAL>&>(p.second);
			const std::vector<point_type>& pts = peers.at(sender).pts();

			assert( pts.size() == data.size() );

			std::vector<std::pair<point_type,REAL> >& data_store = storage_cast<std::vector<std::pair<point_type,REAL> >&>(n);

			data_store.reserve(data.size());

			for( std::size_t i=0; i<data.size(); i++ ){
				data_store.emplace_back( pts[i], data[i] );
			}
		}

		return buf;
	}
};

}

#endif // _UNIFACE_H
