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
 * @file comm_mpi_nxn.h
 * @author S. Kudo
 * @date 4 March 2014
 * @brief Structures and methods for a many-to-many (nxn) communicator type.
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
			for( size_t v : recvsizes ) total += v;

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
