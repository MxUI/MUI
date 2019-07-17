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
 * @file comm_tcp.h
 * @author S. Kudo
 * @date 11 February 2014
 * @brief File containing class definition of base TCP communicator
 */

#ifndef COMM_TCP_H
#define COMM_TCP_H


#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netdb.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <sys/epoll.h>
#include <sys/eventfd.h>
#include <fcntl.h>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <vector>
#include <deque>
#include <list>
#include <exception>
#include <atomic>

#include "message.h"
#include "comm.h"
#include "lib_uri.h"
#include "comm_factory.h"

namespace mui {

/*
 * smalluint is a integer which has variable length depend on
 * its value. smalluint requires small bytes if the value is smaller.
 * For example if the value is less than equal to 127, it only requires
 * one byte.
 * Format: the byte length is determined by the number of leading ones
 * of first byte.
 */
struct smalluint {
	uint64_t num;
	static int detect_size( unsigned char byte ){
		if(      byte <= 0x7fu ) return 1;
		else if( byte <= 0xbfu ) return 2;
		else if( byte <= 0xdfu ) return 3;
		else if( byte <= 0xefu ) return 4;
		else if( byte <= 0xf7u ) return 5;
		else if( byte <= 0xfbu ) return 6;
		else if( byte <= 0xfdu ) return 7;
		else if( byte <= 0xfeu ) return 8;
		else return 9;
	}
};

inline istream& operator>> ( istream& stream, smalluint& sml )
{
	unsigned char byte;
	stream >> byte;
	uint64_t num = 0u;
	int size;
	if(      byte <= 0x7fu ) { size = 1; num = byte&0xffu; }
	else if( byte <= 0xbfu ) { size = 2; num = byte&0x3fu; }
	else if( byte <= 0xdfu ) { size = 3; num = byte&0x1fu; }
	else if( byte <= 0xefu ) { size = 4; num = byte&0x0fu; }
	else if( byte <= 0xf7u ) { size = 5; num = byte&0x07u; }
	else if( byte <= 0xfbu ) { size = 6; num = byte&0x03u; }
	else if( byte <= 0xfdu ) { size = 7; num = byte&0x01u; }
	else if( byte <= 0xfeu ) size = 8;
	else size = 9;

	for( int i=0; i<size-1; ++i ){
		unsigned char cur;
		stream >> cur;
		num = (num<<8u)|cur;
	}
	sml.num = num;
	return stream;
}
inline ostream& operator<< ( ostream& stream, const smalluint& sml )
{
	uint64_t num = sml.num;
	unsigned char bytes[9];
	int size;
	if(      num <= 0x000000000000007full ) size = 1;
	else if( num <= 0x0000000000003fffull ) size = 2;
	else if( num <= 0x00000000001fffffull ) size = 3;
	else if( num <= 0x000000000fffffffull ) size = 4;
	else if( num <= 0x00000007ffffffffull ) size = 5;
	else if( num <= 0x000003ffffffffffull ) size = 6;
	else if( num <= 0x0001ffffffffffffull ) size = 7;
	else if( num <= 0x00ffffffffffffffull ) size = 8;
	else size = 9;

	for( int j=size-1; j>=0; --j ){
		bytes[j] = (num&0xffull);
		num >>= 8u;
	}
	bytes[0] |= (0xfff00ull>>(size-1));
	stream.write((char*)bytes,size);

	return stream;
}

namespace {
static const std::size_t MUI_TCP_BUF_SIZE = 4000;
int SYSCHECK(int value){ // carefully use this. just like "SYSCHECK(close(fd));"
	int err = errno;
	if(value == -1) throw std::system_error(err, std::system_category());
	return value;
}
}

struct mutex_timeout: std::runtime_error {
	mutex_timeout(): std::runtime_error("tcp error: lock time out. May be dead locked.") {}
};

/* unique_fd_
 * RAII class for file-descripter.
 */
class unique_fd_ {
public:
	unique_fd_(int fd=0): fd_(fd) {}
	unique_fd_(const unique_fd_&) = delete;
	unique_fd_(unique_fd_&& rhs): fd_(rhs.fd_) {
		rhs.fd_ = 0;
	}
	unique_fd_& operator=(const unique_fd_&) = delete;
	unique_fd_& operator=(unique_fd_&& rhs){
		checked_close_();
		fd_ = rhs.fd_;
		rhs.fd_ = 0;
		return *this;
	}
	~unique_fd_(){ checked_close_(); }
	operator int () const { return fd_; }
	void reset(int fd) {
		checked_close_();
		fd_ = fd;
	}
	void swap(unique_fd_& rhs) { std::swap(fd_, rhs.fd_); }
private:
	void checked_close_(){ if( fd_ ) close(fd_); }
	int fd_;
};

/* poll_scheduler
 * A wrapper for epoll. It calls callback functions if the fd binded to callback becomes ready.
 * port to select? kqueue?
 */
class poll_scheduler {
public:
	static const uint64_t IN = 1;
	static const uint64_t OUT = 4;
	static const uint64_t CALL = 32;
	static const uint64_t ET = 1u<<31; // edge triggered
	
	typedef std::function<void(uint64_t)> callback_t;
	
	poll_scheduler(int max){
		buf_.reserve(max+1);
		callbacks_.reserve(max+1);

		epfd_ = SYSCHECK(epoll_create(max));

		int pipes[2];
		SYSCHECK(pipe(pipes));
		pipe_read_.reset(pipes[0]);
		pipe_write_.reset(pipes[1]);
		add(pipe_read_, IN, std::bind(&poll_scheduler::pop_event_, this));
	}
	
	~poll_scheduler() {
		for( auto& a: callbacks_ )
			SYSCHECK(epoll_ctl(epfd_, EPOLL_CTL_DEL, a.first, NULL));
	}

	poll_scheduler(const poll_scheduler&) = delete;
	poll_scheduler& operator=(const poll_scheduler&) = delete;
	
	void run() {
		buf_.resize(callbacks_.size());
		int nfd = epoll_wait(epfd_, buf_.data(), buf_.size(), 100);
		int err = errno;
		if( nfd == 0 ) return;;
		if( nfd == -1 && err == EAGAIN ) return;
		else if( nfd == -1 ) throw std::system_error(err,std::system_category());
		for( int i=0; i<nfd; ++i )
			callbacks_[buf_[i].data.fd](translate_(buf_[i].events));
	}

	void add(int fd, uint64_t default_events, callback_t callback) {
		callbacks_.emplace(fd, std::move(callback));
		epoll_event ev;
		ev.data.fd = fd;
		ev.events = translate_from_(default_events) | EPOLLET;
		try {
			SYSCHECK(epoll_ctl(epfd_, EPOLL_CTL_ADD, fd, &ev));
		} catch(...) {
			callbacks_.erase(callbacks_.find(fd));
			throw;
		}
	}
	void schedule(int fd) {
		SYSCHECK(write(pipe_write_, &fd, sizeof(int)));
	}
private:
	void pop_event_() {
		int fd;
		SYSCHECK(read(pipe_read_, &fd, sizeof(int)));
		if( fd != -1 )
			callbacks_[fd](CALL);
	}
	uint64_t translate_(uint64_t ev) {
		return (ev&EPOLLIN?IN:0u) | (ev&EPOLLOUT?OUT:0u) | (ev&EPOLLET?ET:0u);
	}
	uint64_t translate_from_(uint64_t ps_ev) {
		return (ps_ev&IN?EPOLLIN:0u) | (ps_ev&OUT?EPOLLOUT:0u) |(ps_ev&ET?EPOLLET:0u);
	}
	
	std::unordered_map<int, callback_t> callbacks_;
	std::vector<epoll_event> buf_;
	unique_fd_ epfd_;
	unique_fd_ pipe_read_, pipe_write_;
};

namespace {
/* read_buffer
 * This manages buffers which will be used by read(3).
 * Internal use only.
 */
struct read_buffer : istream {
	static const int BUFSIZE = MUI_TCP_BUF_SIZE;
	struct node_t { char buf[BUFSIZE]; };
	
	std::pair<int,char*> get_buffer() {
		if(nodes_.empty() || tail_ == BUFSIZE) {
			nodes_.emplace_back();
			tail_ = 0u;
		}
		return std::make_pair(BUFSIZE-tail_, nodes_.back().buf+tail_);
	}
	void move_tail(int size) { // consume bytes. size must be <= get_buffer().first
		psize_ += size;
		size_ += size;
		tail_ += size;
	}

	void read(char* dest, std::size_t size) { // read bytes and seek
		size_ -= size;
		while(size) {
			std::size_t sz = std::min<std::size_t>(size, (pcur_!=nodes_.size()-1?BUFSIZE-ccur_:tail_-ccur_));
			dest = std::copy_n(nodes_[pcur_].buf+ccur_, sz, dest);
			size -= sz;
			ccur_ += sz;
			if( size ) {
				++pcur_;
				ccur_ = 0u;
			}
		}
	}
	void revert() { // undo seeks
		pcur_ = 0u;
		ccur_ = head_;
		size_ = psize_;
	}
	void detach() { // erase the used bytes
		nodes_.erase(nodes_.begin(), nodes_.begin()+pcur_);
		head_ = ccur_;
		psize_ = size_;
	}
	std::size_t size() const { return size_; }
	
	std::deque<node_t> nodes_;
	std::size_t psize_ = 0u;
	std::size_t size_ = 0u;
	std::size_t pcur_ = 0u;
	unsigned ccur_ = 0u;
	unsigned head_ = 0u;
	unsigned tail_ = 0u;
};
}

class read_que {
public:
	read_que( unique_fd_&& fd, std::function<void(message)> callback )
	: fd_(std::move(fd)), callback_(std::move(callback)) {}
	read_que( read_que&& ) = default;
	read_que& operator=( read_que&& ) = default;
	
	void try_recv(int){
		while(true){
			std::pair<int,char*> p = buf_.get_buffer();
			ssize_t r = read(fd_, p.second, p.first);
			int err = errno;
			if( r > 0 ){
				buf_.move_tail(r);
				char c;
				buf_ >> c;
				buf_.revert();
				std::uint64_t smlsize = smalluint::detect_size(c);
				if( buf_.size() < smlsize) continue;

				smalluint vecsize;
				buf_ >> vecsize;
				buf_.revert();
				if( buf_.size() < smlsize + vecsize.num ) continue;

				message msg;
				buf_ >> msg;
				callback_(std::move(msg));
				buf_.detach();
				continue;
			} else if( r == -1 && (err == EAGAIN || err == EWOULDBLOCK) ){
				break;
			}else if( r == -1 ) throw std::system_error(err,std::system_category());
			else throw std::runtime_error("tcp connection is broken.\n");
		}
	}
	int get_fd() const { return fd_; }
private:
	unique_fd_ fd_;
	std::function<void(message)> callback_;
	read_buffer buf_;
};

class write_que {
public:
	explicit write_que(unique_fd_&& fd): fd_(std::move(fd)) {}
	write_que(write_que&& rhs): fd_(std::move(rhs.fd_)), send_(std::move(send_)) {}
	write_que& operator=(write_que&& rhs) {
		std::swap(fd_, rhs.fd_);
		send_.swap(rhs.send_);
		return *this;
	}
	
	void try_send(int){
		while(true) {
			mutex_.lock();
			auto iter = send_.begin();
			auto end  = send_.end();
			mutex_.unlock();
			if(iter == end) break;
			auto& p = *iter;

			std::size_t sz = p.second.size() - p.first;
			char* buf = p.second.data() + p.first;
			ssize_t r = write(fd_, buf, sz);
			int err = errno;
			if( r == sz ) {
				std::lock_guard<std::mutex> lock(mutex_);
				send_.pop_front();
			} else if( r > 0 ) {
				p.first += r;
			} else if( r == -1 && (err == EAGAIN || err == EWOULDBLOCK) ){
				break;
			} else if( r == -1 ) throw std::system_error(err,std::system_category());
			else throw std::runtime_error("tcp connection is broken.\n");
		}
	}
	void push( std::vector<char> data ){
		if( data.size() > std::numeric_limits<ssize_t>::max() ) 
			throw std::range_error("cannot send data larger than ssize_t");
		if( !data.empty() ) {
			std::lock_guard<std::mutex> lock(mutex_);
			send_.emplace_back(0u, std::move(data));
		}
	}
	int get_fd() const { return fd_; }
private:
	unique_fd_ fd_;
	std::mutex mutex_;
	std::list<std::pair<std::size_t,std::vector<char> > > send_;
};

class comm_fd: public communicator {
public:
	comm_fd(int local_rank, int local_size, std::vector<unique_fd_>&& wfds, std::vector<unique_fd_>&& rfds)
	: local_rank_(local_rank), local_size_(local_size), remote_size_(wfds.size()), poll_(rfds.size()+wfds.size()) {
		assert(wfds.size() == rfds.size());
		die_.store(false);

		rques_.reserve(rfds.size());
		wques_.reserve(wfds.size());

		for(auto& r: rfds) rques_.emplace_back(std::move(r),std::bind(&comm_fd::push_msg_,this,std::placeholders::_1));
		for(auto& q: rques_)
			poll_.add(q.get_fd(), poll_scheduler::IN | poll_scheduler::ET, std::bind(&read_que::try_recv,&q,std::placeholders::_1));
		for(auto& w: wfds) wques_.emplace_back(std::move(w));
		for(auto& q: wques_)
			poll_.add(q.get_fd(), poll_scheduler::OUT | poll_scheduler::ET, std::bind(&write_que::try_send,&q,std::placeholders::_1));

		std::thread th = std::thread(std::bind(&comm_fd::run_, this));
		thread_.swap(th);
	}
	
	~comm_fd() {
		die_.store(true);
		poll_.schedule(-1);
		thread_.join();
	}
	
	int local_rank() const { return local_rank_; }
	int local_size() const { return local_size_; }
	int remote_size() const { return remote_size_; }
protected:
	void send_impl_(message msg, const std::vector<bool>& dest) {
		std::vector<char> v(streamed_size(msg));
		auto stream = make_ostream(v.begin());
		stream << msg;
		for( int i=0; i<remote_size(); ++i ) {
			if(dest[i]){
				wques_[i].push(v);
				poll_.schedule(wques_[i].get_fd());
			}
		}
	}
	message recv_impl_() {
		std::unique_lock<std::mutex> lock(recv_mutex_);
		recv_cv_.wait(lock, [=](){ return !mesgs_.empty(); });
		message msg = std::move(mesgs_.front());
		mesgs_.pop_front();
		return msg;
	}
private:
	void push_msg_(message msg) {
		std::unique_lock<std::mutex> lock(recv_mutex_);
		mesgs_.emplace_back(std::move(msg));
		recv_cv_.notify_one(); // only ONE thrad can cann comm_fd.recv_impl_();
	}

	void run_() {
		while(true){
			if( die_ ) break;
			poll_.run();
		}
	}

	int local_rank_, local_size_, remote_size_;

	std::thread thread_;	
	std::atomic_bool die_;

	poll_scheduler poll_;
	std::vector<read_que> rques_;
	std::vector<write_que> wques_;

	std::mutex recv_mutex_;
	std::condition_variable recv_cv_;
	std::list<message> mesgs_;
};

inline communicator* create_comm_tcp( const char* str ){
	uri u(str);
	bool is_server = u.host().empty();

	struct free_addrinfo_ { void operator()( addrinfo* ptr ) const { freeaddrinfo(ptr); } };
	unique_fd_ sock;
	std::unique_ptr<addrinfo,free_addrinfo_> info;
	addrinfo hint, *tmp=NULL;

	std::memset((char*)&hint, 0, sizeof(hint));
	hint.ai_family = AF_INET;	// use ip v4
	hint.ai_socktype = SOCK_STREAM;	// TCP/IP
	hint.ai_flags = AI_NUMERICSERV|(is_server?AI_PASSIVE:0);
	
	if(int err = getaddrinfo(is_server?0:u.host().c_str(), u.path().c_str(), &hint, &tmp))
		throw std::runtime_error(gai_strerror(err));
	info.reset(tmp);
	
	for(; tmp!=NULL; tmp = tmp->ai_next){
		int s = socket(tmp->ai_family, tmp->ai_socktype, tmp->ai_protocol);
		if(s == -1 ) continue;
		sock.reset(s);
		
		if( is_server ){
			if( bind(sock, tmp->ai_addr, tmp->ai_addrlen) != -1 ) break;
		} else {
			if( connect(sock, tmp->ai_addr, tmp->ai_addrlen) != -1 ) break;
		}
	}
	if(tmp == NULL) throw std::runtime_error("tcp connection error.");
	
	unique_fd_ wfd, rfd;
	
	if( is_server ){
		SYSCHECK(listen(sock,1));
		wfd.reset(SYSCHECK(accept(sock,0,0)));
		SYSCHECK(fcntl(wfd,F_SETFL,O_NONBLOCK));
		rfd.reset(SYSCHECK(dup(wfd)));
	} else {
		SYSCHECK(fcntl(sock,F_SETFL,O_NONBLOCK));
		wfd.swap(sock);
		rfd.reset(SYSCHECK(dup(wfd)));
	}
	std::vector<mui::unique_fd_> w; w.emplace_back(std::move(wfd));
	std::vector<mui::unique_fd_> r; r.emplace_back(std::move(rfd));
	return static_cast<communicator*>(new comm_fd(0,1,std::move(w),std::move(r)));
}

const static bool comm_tcp_registered_ = comm_factory::instance().link( "tcp", create_comm_tcp );

inline communicator* create_comm_shm( const char* str ){
	uri u(str);
	std::string patha = std::string("/dev/shm/") + u.host() + "_A";
	std::string pathb = std::string("/dev/shm/") + u.host() + "_B";
	bool is_wr = false;
	
	if(mkfifo(patha.c_str(), S_IRUSR | S_IWUSR)) {
		int err = errno;
		if(err != EEXIST) throw std::system_error(err, std::system_category());
		is_wr = true;
	}
	if(!is_wr)
		if(int ret = mkfifo(pathb.c_str(), S_IRUSR | S_IWUSR))
			throw std::system_error(errno, std::system_category());
	
	unique_fd_ afd = SYSCHECK(open(patha.c_str(), is_wr?O_RDONLY:O_WRONLY));
	SYSCHECK(fcntl(afd,F_SETFL,O_NONBLOCK));
	unique_fd_ bfd = SYSCHECK(open(pathb.c_str(), is_wr?O_WRONLY:O_RDONLY));
	SYSCHECK(fcntl(bfd,F_SETFL,O_NONBLOCK));
	
	std::vector<mui::unique_fd_> w; w.emplace_back(std::move(is_wr?bfd:afd));
	std::vector<mui::unique_fd_> r; r.emplace_back(std::move(is_wr?afd:bfd));
	return static_cast<communicator*>(new comm_fd(0,1,std::move(w),std::move(r)));
}

const static bool comm_shm_registered_ = comm_factory::instance().link( "shm", create_comm_shm);

}

#if 0
void print( mui::message msg )
{
	auto stream = mui::make_istream(msg.data());
	std::string str;
	stream >> str;
	std::cerr<< "GET MESSAGE!!\t" << str << std::endl;
}

int main(int argc, char **argv)
{
	//std::unique_ptr<mui::communicator> c(mui::create_comm_tcp(argc<2?"tcp:///37129":"tcp://localhost/37129"));
	std::unique_ptr<mui::communicator> c(mui::create_comm_shm("shm://my_pipe/"));

	c->send(mui::message::make("test",std::string("Hello, ")));
	mui::message msg = c->recv();
	print(msg);
	sleep(3);
	return 0;
}
#endif

#endif
