/* -*- c++ -*-
 * comm_mpi_nxn.h
 *
 *  Created on: Mar 4, 2014
 *      Author: skudo
 */

#ifndef COMM_MPI_NXN_H
#define COMM_MPI_NXN_H

#include <algorithm>
#include <list>
#include <deque>
#include <vector>
#include <mpi.h>
#include "comm.h"

namespace mui {
class mpicomm_nxn : public communicator {
private:
	int tag;
	MPI_Comm local_comm;
	MPI_Comm intercomm;
	std::list<std::pair<MPI_Request,std::vector<char> > > bufs;
	std::deque<char> recv_buf;
public:
	mpicomm_nxn( int tag_, MPI_Comm local_comm_, MPI_Comm remote_comm )
		: tag(tag_), local_comm(local_comm_)
	{
		MPI_Intercomm_create( local_comm_, 0, remote_comm, 0, tag, &intercomm );
	}
	virtual ~mpicomm_nxn(){}
	virtual int local_size() const {
		int size;
		MPI_Comm_size(local_comm,&size);
		return size;
	}
	virtual int remote_size() const {
		int size;
		MPI_Comm_remote_size(inter_comm,&size);
		return size;
	}

protected:
	std::size_t send_impl_( const char* ptr, std::size_t size, std::vector<bool> is_sending ) {
		_M_test();
		if( local_rank == 0 ){
			std::vector<int> recvsizes(local_size());
			MPI_Gather((void*)&size, 1, MPI_INT,
			           (void*)recvsizes.data(), 1, MPI_INT,
			           0, local_comm);
			std::vector<int> displ(local_size());
			for( std::size_t i=0; i<displ.size(); ++i ) displ[i] = (!i ? 0: displ[i-1]+recvsizes[i-1]);
			int total=0;
			for( int v : recvsizes ) total += v;

			bufs.push_back(std::make_pair(MPI_Request(),std::vector<char>(total)));
			MPI_Gatherv(const_cast<char*>(ptr), size, MPI_BYTE,
			            bufs.back().second.data(), recvsizes.data(), displ.data(), MPI_BYTE, 0, local_comm );
			MPI_Isend(bufs.back().second.data(), total, MPI_BYTE, 0, tag, intercomm, &(bufs.back().first));
		} else {
			MPI_Gather((void*)&size, 1, MPI_INT,
			           NULL, 0, MPI_INT, 0, local_comm);
			MPI_Gatherv(const_cast<char*>(ptr), size, MPI_BYTE,
			            0, 0, 0, MPI_BYTE, 0, local_comm );
		}
		return size;
	}
	std::size_t recv_impl_( char* ptr, std::size_t size ) {
		_M_test();
		int count;
		std::vector<char> temp;
		if( local_rank() ==  0 ){
			MPI_Status status;
			MPI_Probe(0, tag, intercomm, &status);
			MPI_Get_count(&statuc,MPI_BYTE,&count);
			temp.resize(count);
			MPI_Recv(temp.data(),count,MPI_BYTE,0,tag,intercomm,MPI_IGNORE);
		}
		MPI_Bcast(&count,1,MPI_INT,0,local_comm);
		if( local_rank() != 0 ) temp.resize(count);
		MPI_Bcast(temp.data(),count,MPI_BYTE,0,local_comm);
		recv_buf.insert(recv_buf.end(), temp.begin(), temp.end());

		std::size_t ret_size = std::min(size,recv_buf.size());
		std::copy_n(recv_buf.begin(),ret_size,ptr);
		return ret_size;
	}
	int try_recv_impl_( char* ptr, std::size_t size, std::size_t* received_size ){
		// not implemented yet
		throw("mui::mpicomm::try_recv_impl_ is not implemented yet.");
	}
private:
	int local_rank() const {
		int rank;
		MPI_Comm_rank(local_comm,&rank);
		return rank;
	}
	
	void _M_test() {
		if( local_rank() == 0 ){
			int test = true;
			while( test&&!bufs.empty() ) {
				MPI_Test(&(bufs.front().first),&test,MPI_STATUS_IGNORE);
				if( test ) bufs.pop_front();
			}
		}
	}

};
}
#endif
