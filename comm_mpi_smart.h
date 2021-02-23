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
 * @file comm_mpi_smart.h
 * @author Y. H. Tang
 * @date 12 March 2014
 * @brief Structures and methods for a smart (communication reducing)
 * communicator type.
 */

#ifndef COMM_MPI_SMART_H
#define COMM_MPI_SMART_H

#include <cstdint>
#include <cstdlib>
#include "util.h"
#include "comm.h"
#include "comm_mpi.h"
#include "comm_factory.h"
#include "message.h"
#include "stream.h"

namespace mui {

class comm_mpi_smart : public comm_mpi {
private:
	std::list<std::pair<MPI_Request,std::shared_ptr<std::vector<char> > > > send_buf;
public:
	comm_mpi_smart( const char URI[], MPI_Comm world = MPI_COMM_WORLD ) : comm_mpi(URI, world) {}
	virtual ~comm_mpi_smart() {}

private:
	void send_impl_( message msg, const std::vector<bool> &is_sending ) {
		auto bytes = std::make_shared<std::vector<char> >(msg.detach());
		for( int i = 0 ; i < remote_size_ ; i++ ){
			if( is_sending[i] ){
				if(bytes->size() > INT_MAX){
					std::cerr << "MUI Error [comm_mpi_smart.h]: Trying to send more data than is possible with MPI_Isend." << std::endl
							  << "This is likely because there is too much data per MPI rank." << std::endl
							  << "The program will now abort. Try increasing the number of MPI ranks." << std::endl;
					std::abort();
				}
				send_buf.emplace_back(MPI_Request(), bytes);
				MPI_Isend(bytes->data(), bytes->size(), MPI_BYTE, i, 0, 
				          domain_remote_, &(send_buf.back().first));
			}
		}
		// Block until MPI_Isend calls completed successfully, maximum 60 seconds with printed warnings
		test_completion_blocking();
	}

	message recv_impl_() {
		// Catch any unsent MPI_Isend calls before fetch starts, non-blocking
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

	/** \brief Non-blocking check for complete MPI_Isend calls
     */
	void test_completion() {
		for( auto itr=send_buf.begin(), end=send_buf.end(); itr != end; ){
			int test = false;
			MPI_Test(&(itr->first),&test,MPI_STATUS_IGNORE);
			if( test ) itr = send_buf.erase(itr);
			else ++itr;
		}
	}

	/** \brief Time-limited blocking check for complete MPI_Isend calls
     */
	void test_completion_blocking() {
		int timeout = 0;
		auto start = std::chrono::system_clock::now();
		while (send_buf.size() > 0) {
			for( auto itr=send_buf.begin(), end=send_buf.end(); itr != end; ){
				int test = false;
				MPI_Test(&(itr->first),&test,MPI_STATUS_IGNORE);
				if( test ) itr = send_buf.erase(itr);
				else ++itr;
			}
			if( (std::chrono::system_clock::now() - start) > std::chrono::seconds(5) ) {
				timeout++;
				int time_left = 60 - (timeout*5);
				if( time_left == 0) break;
				std::cout << "MUI Warning [comm_mpi_smart.h]: MPI_Isend calls spent over 5 seconds completing, all operations will timeout in "
						  << time_left << "s" << std::endl;
			}
		}
	}
};

inline communicator *create_comm_mpi_smart( const char URI[] ) {
	return new comm_mpi_smart(URI);
}

const static bool registered = comm_factory::instance().link( "mpi", create_comm_mpi_smart );

}
#endif
