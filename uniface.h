/* -*- c++ -*-
 * uniface.h
 *
 *  Created on: Feb 11, 2014
 *      Author: skudo
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
#include "stream.h"
#include "stream_vector.h"
#include "stream_unordered.h"
#include "stream_string.h"
#include "bin.h"

namespace mui {

/*! \brief Briefly what is uniface
 *
 *  An expanded description of uniface. 
 */
template<typename CONFIG = default_config>
class uniface {
public:
	// public typedefs.
	static const int  D     = CONFIG::D; //!< dimensionality of the domains
	static const bool DEBUG = CONFIG::DEBUG;
	using REAL       = typename CONFIG::REAL;
	using point_type = typename CONFIG::point_type;
	using time_type  = typename CONFIG::time_type;
	using data_types = typename CONFIG::data_types;

	using span_t = geometry::any_shape<CONFIG>;
private:
	// meta functions to split tuple and add vector<pair<point_type,_1> >
	template<typename T> struct add_vp_  { using type = std::vector<std::pair<point_type,T> >; };
	template<typename T> struct def_storage_;
	template<typename... TYPES> struct def_storage_<type_list<TYPES...> >{
		using type = storage<typename add_vp_<TYPES>::type...>;
	};

	// internal typedefinitions
	using storage_t = typename def_storage_<data_types>::type;
	using spatial_t = spatial_storage<bin_t<CONFIG>,storage_t,CONFIG>;
	using frame_type = std::unordered_map<std::string, storage_t>;
	using bin_frame_type = std::unordered_map<std::string, spatial_t>;

	struct peer_state {
		using spans_type = std::map<std::pair<time_type,time_type>,span_t>;

		bool is_recving(time_type t, const span_t& s) const {
			return scan_spans_(t,s,recving_spans);
		}
		void set_recving( time_type start, time_type timeout, span_t s ) {
			recving_spans.emplace(std::make_pair(start,timeout),std::move(s));
		}
		bool is_sending(time_type t, const span_t& s) const {
			return scan_spans_(t,s,sending_spans);
		}
		void set_sending(time_type start, time_type timeout, span_t s) {
			sending_spans.emplace(std::make_pair(start,timeout), std::move(s));
		}	

		time_type current_t() const { return latest_timestamp; }
		time_type next_t() const { return next_timestamp; }
		void set_current_t( time_type t ) { latest_timestamp = t; }
		void set_next_t( time_type t ) { next_timestamp = t; }
	private:
		bool scan_spans_(time_type t, const span_t& s, const spans_type& spans ) const {
			auto p = std::make_pair(t,t);
			auto end = spans.upper_bound(p);
			bool prefetched = false;
			
			for( auto itr = spans.begin(); itr != end; ++itr )
				if( t <= itr->first.second ){
					prefetched = true;
					if( collide(s,itr->second) ) return true;
				}
			// if prefetched at t, but no overlap region, then return false;
			// otherwise return true;
			return !prefetched;
		}
		time_type latest_timestamp = std::numeric_limits<time_type>::lowest();
		time_type next_timestamp = std::numeric_limits<time_type>::lowest();
		spans_type recving_spans;
		spans_type sending_spans;
	};

private: // data members
	std::unique_ptr<communicator> comm;
	dispatcher<message::id_type, std::function<void(message)> > readers;

	frame_type push_buffer;
	std::map<time_type,bin_frame_type> log;

	std::vector<peer_state> peers;
	time_type span_start   = std::numeric_limits<time_type>::lowest();
	time_type span_timeout = std::numeric_limits<time_type>::lowest();
	span_t current_span;
	time_type recv_start   = std::numeric_limits<time_type>::lowest();
	time_type recv_timeout = std::numeric_limits<time_type>::lowest();
	span_t recv_span;
	time_type memory_length = std::numeric_limits<time_type>::max();
	std::mutex mutex;
public:
	uniface( const char URI[] ) : uniface( comm_factory::create_comm(URI) ) {}
	uniface( std::string const &URI ) : uniface( comm_factory::create_comm(URI.c_str()) ) {}
	uniface( communicator* comm_ ) : comm(comm_) {
		using namespace std::placeholders;

		peers.resize(comm->remote_size());
		readers.link("timestamp", reader_variables<int32_t, time_type>(
		             std::bind(&uniface::on_recv_confirm, this, _1, _2)));
		readers.link("forecast", reader_variables<int32_t, time_type>(
		             std::bind(&uniface::on_recv_forecast, this, _1, _2)));
		readers.link("receiving span", reader_variables<int32_t, time_type, time_type, span_t>(
		             std::bind(&uniface::on_recv_span, this, _1, _2, _3, _4)));
		readers.link("sending span", reader_variables<int32_t, time_type, time_type, span_t>(
				std::bind(&uniface::on_send_span, this, _1, _2, _3, _4)));
		readers.link("data", reader_variables<time_type, frame_type>(
		             std::bind(&uniface::on_recv_data, this, _1, _2)));
	}

	uniface( const uniface& ) = delete;
	uniface& operator=( const uniface& ) = delete;

	template<typename TYPE>
	void push( const std::string& attr, const point_type loc, const TYPE value ) {
		storage_t& n = push_buffer[attr];
		if( !n ) n = storage_t(std::vector<std::pair<point_type,TYPE> >());
		storage_cast<std::vector<std::pair<point_type,TYPE> >&>(n).emplace_back( loc, value );
	}
	
	template<class SAMPLER, class TIME_SAMPLER>
	typename SAMPLER::OTYPE
	fetch( const std::string& attr, const point_type focus, const time_type t,
	       const SAMPLER& sampler, const TIME_SAMPLER &t_sampler ) {
		barrier(t_sampler.get_upper_bound(t));
		std::vector<std::pair<time_type,typename SAMPLER::OTYPE> > v;

		for( auto first=log.lower_bound(t_sampler.get_lower_bound(t)),
		     last = log.upper_bound(t_sampler.get_upper_bound(t)); first!= last; ++first ){
			time_type time = first->first;
			auto iter = first->second.find(attr);
			if( iter == first->second.end() ) continue;
			v.emplace_back(time,iter->second.build_and_query_ts(sampler.support(focus).bbox(),focus,sampler));
		}
		return t_sampler.filter(t, v);
	}

	// commit() serializes pushed data and send it to remote nodes and, after that,
	// send "confirm" message.
	// return actual number of peers contacted
	int commit( time_type timestamp ) {
		std::vector<bool> is_sending(comm->remote_size(), true);
		if( span_start <= timestamp  && timestamp <= span_timeout )
			for( std::size_t i=0; i<peers.size(); ++i ) {
				is_sending[i] = peers[i].is_recving( timestamp, current_span );
			}
		comm->send(message::make("data",timestamp,std::move(push_buffer)), is_sending);
		push_buffer.clear();
		comm->send(message::make("timestamp",comm->local_rank(),timestamp));
		return std::count( is_sending.begin(), is_sending.end(), true );
	}
	void forecast( time_type timestamp ) {
		comm->send(message::make("forecast", comm->local_rank(), timestamp));
	}

	bool is_ready( const std::string& attr, time_type t ) const {
		return  std::any_of(log.begin(), log.end(), [=](const bin_frame_type& a) {
				return a.find(attr) != a.end(); }) // return false for random @attr. 
			&& std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) {
				return (!p.is_sending(t,recv_span)) || (p.current_t() >= t || p.next_t() > t); });
	}
	void barrier( time_type t ) {
		auto start = std::chrono::system_clock::now();
		for(;;) {    // barrier must be thread-safe because it is called in fetch()
			std::lock_guard<std::mutex> lock(mutex);
			if( std::all_of(peers.begin(), peers.end(), [=](const peer_state& p) { 
				return (!p.is_sending(t,recv_span)) || (p.current_t() >= t || p.next_t() > t); }) ) break;
			acquire(); // To avoid infinite-loop when synchronous communication
		}
		if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) )
			std::cerr << "MUI Warning. Communication spends over 5 sec." << std::endl;
	}

	void announce_send_span( time_type start, time_type timeout, span_t s ){
		// say remote nodes that "I'll send this span."
		// not implemented yet. just save arguments.
		comm->send(message::make("sending span", comm->local_rank(), start, timeout, std::move(s)));
		span_start = start;
		span_timeout = timeout;
		current_span.swap(s);
	}
	void announce_recv_span( time_type start, time_type timeout, span_t s ){
		// say remote nodes that "I'm recving this span."
		comm->send(message::make("receiving span", comm->local_rank(), start, timeout, std::move(s)));
		recv_start = start;
		recv_timeout = timeout;
		recv_span.swap(s);
	}

	// remove log betwewn (-inf, @end]
	void forget( time_type end ) {
		log.erase(log.begin(), log.upper_bound(end));
	}
	// remove log between [@first, @last]
	void forget( time_type first, time_type last ) {
		log.erase(log.lower_bound(first), log.upper_bound(last));
	}
	// remove log between (-inf, curent-@length] automatically. 
	void set_memory( time_type length ) {
		memory_length = length;
	}
	
private:
	// triggers communication
	void acquire() {
		message m = comm->recv();
		if( m.has_id() ) readers[m.id()](m);
	}

	void on_recv_confirm( int32_t sender, time_type timestamp ) {
		peers.at(sender).set_current_t(timestamp);
	}
	void on_recv_forecast( int32_t sender, time_type timestamp ) {
		peers.at(sender).set_next_t(timestamp);
	}

	void on_recv_data( time_type timestamp, frame_type frame ) {
		// when message.id_ == "data"
		auto itr = log.find(timestamp);
		if( itr == log.end() ) std::tie(itr,std::ignore) = log.insert(std::make_pair(timestamp,bin_frame_type()));
		auto& cur = itr->second;
		for( auto& p: frame ){
			auto pstr = cur.find(p.first);
			if( pstr == cur.end() ) cur.insert(std::make_pair(std::move(p.first),spatial_t(std::move(p.second))));
			else pstr->second.insert(p.second);
		}
		log.erase(log.begin(), log.upper_bound(timestamp-memory_length));
	}
	void on_recv_span( int32_t sender, time_type start, time_type timeout, span_t s ) {
		peers.at(sender).set_recving(start,timeout,std::move(s));
	}
	void on_send_span( int32_t sender, time_type start, time_type timeout, span_t s ){
		peers.at(sender).set_sending(start,timeout,std::move(s));
	}
};

}

#endif // _UNIFACE_H

// TO-DO
// [x] in uniface::barrier: give time-out warning  <- requires try_recv
// [x] void uniface::forget( time_type time );
// [x] void uniface::forget( time_type begin, time_type end );
// [x] void uniface::set_memory( time_type length );
// [x] bool uniface::is_ready( std::string attr );
// [x] void uniface::forecast( time_type stamp );
// [ ] linear solver
// [ ] config generator
// [ ] logger
// [ ] python wrapper
// [ ] periodicity policy

// [ ] shm://domain/interface, create pipes in /dev/shm/interface
// [ ] void uniface::recommit( time_type stamp, INT version );
