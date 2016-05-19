/* -*- c++ -*-
 * comm_mpi_smart.h
 *
 *  Created on: Mar 12, 2014
 *      Author: ytang
 */

#ifndef COMM_MPI_SMART_H
#define COMM_MPI_SMART_H

#include "util.h"
#include "comm.h"
#include "comm_mpi.h"
#include "comm_factory.h"
#include "stream.h"
#include "message.h"

namespace mui {

class comm_mpi_smart : public comm_mpi {
private:
	std::list<std::pair<MPI_Request,std::shared_ptr<std::vector<char> > > > send_buf;
public:
	comm_mpi_smart( const char URI[], MPI_Comm world = MPI_COMM_WORLD ) : comm_mpi(URI, world) {}
	virtual ~comm_mpi_smart() {}

private:
	void send_impl_( message msg, const std::vector<bool> &is_sending ) {
		test_completion();
		auto bytes = std::make_shared<std::vector<char> >(msg.detach());
		for(int i = 0 ; i < remote_size_ ; i++)
			if ( is_sending[i] ) {
				send_buf.emplace_back(MPI_Request(), bytes);
				MPI_Isend(bytes->data(), bytes->size(), MPI_BYTE, i, 0, 
				          domain_remote_, &(send_buf.back().first));
			}
	}

	message recv_impl_() {
		test_completion();
		std::vector<char> temp;
		MPI_Status status;
		MPI_Probe(MPI_ANY_SOURCE, 0, domain_remote_, &status);
		int count;
		MPI_Get_count(&status,MPI_BYTE,&count);
		temp.resize(count);
		MPI_Recv( temp.data(), count, MPI_BYTE, status.MPI_SOURCE, status.MPI_TAG, domain_remote_, MPI_STATUS_IGNORE );
		return message::make(std::move(temp));
	}

	void test_completion() {
		for( auto itr=send_buf.begin(), end=send_buf.end(); itr != end; ){
			int test = false;
			MPI_Test(&(itr->first),&test,MPI_STATUS_IGNORE);
			if( test ) itr = send_buf.erase(itr);
			else ++itr;
		}
	}
};

inline communicator *create_comm_mpi_smart( const char URI[] ) {
	return new comm_mpi_smart(URI);
}

const static bool registered = comm_factory::instance().link( "mpi", create_comm_mpi_smart );

}
#endif
