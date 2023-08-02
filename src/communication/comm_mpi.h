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
 * @file comm_mpi.h
 * @author S. Kudo
 * @date 11 February 2014
 * @brief Class definition of base MPI communicator.
 */

#ifndef COMM_MPI_H_
#define COMM_MPI_H_

#include <mpi.h>
#include "comm.h"
#include "../general/util.h"
#include "lib_uri.h"

namespace mui {

/*
 * MPI base implementation of Comm
 */

class comm_mpi : public communicator {
public:
	comm_mpi( const char URI[], const bool quiet, MPI_Comm world ) :
		domain_local_(0),
		domain_remote_(0),
		local_size_(0),
		local_rank_(0),
		remote_size_(0),
		global_size_(0),
		global_rank_(0),
		uri_host_(std::string()),
		uri_path_(std::string()),
		uri_protocol_(std::string()),
		initialized(false),
		init_by_me(false) {
		init(URI, quiet, world);
	}
	virtual ~comm_mpi() {
		finalize();
	}

	void init( const char URI[], const bool quiet, MPI_Comm world ) {
		if ( initialized )
			throw( std::runtime_error("MUI Error [comm_mpi.h]: Duplicate MUI communicator initialization") );

		// check MPI initialization status
		int init;
		MPI_Initialized( &init );
		if (!init) {
			init_by_me = true;
			MPI_Init( NULL, NULL );
		}

		// duplicate world to avoid potential conflict with host program
		if (world == MPI_COMM_WORLD) MPI_Comm_dup( MPI_COMM_WORLD, &world );
		MPI_Comm_size( world, &global_size_ );
		MPI_Comm_rank( world, &global_rank_ );
		if (global_size_ < 2)
			std::cout <<  "MUI Warning [comm_mpi.h]: Number of global ranks less than 2" << std::endl;

		// get upper bond for tag hash
		int prime;
		{
			int* v;
			int flag;
			MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_TAG_UB,&v,&flag);
			prime = *v;
		}

		// parse URI, split the world using domain tag
		uri desc(URI);
		uri_host_ = desc.host();
		uri_path_ = desc.path();
		uri_protocol_ = desc.protocol();
		int domain_hash = std::hash<std::string>()( desc.host() );
		domain_hash=std::abs(domain_hash);
		int ifs_hash = std::hash<std::string>()( desc.path() );
		ifs_hash=std::abs(ifs_hash);
		MPI_Comm_split( world, domain_hash % prime, global_rank_, &domain_local_ );
		MPI_Comm_size( domain_local_, &local_size_ );
		MPI_Comm_rank( domain_local_, &local_rank_ );

		// set intercommunication using interface tag
		std::vector<int> all_domain( global_size_, 0 );
		std::vector<int> all_ifs   ( global_size_, 0 );
		MPI_Allgather( &domain_hash, 1, MPI_INT, &all_domain[0], 1, MPI_INT, world );
		MPI_Allgather( &ifs_hash,    1, MPI_INT, &all_ifs[0],    1, MPI_INT, world );
		for( int i = 0 ; i < global_size_ ; i++ ) {
			if ( i == global_rank_ ) continue;
			if ( all_domain[i] != domain_hash && all_ifs[i] == ifs_hash ) {
				MPI_Intercomm_create( domain_local_, 0, world, i, ifs_hash % prime, &domain_remote_ );
				break;
			}
		}
		MPI_Comm_remote_size( domain_remote_, &remote_size_ );

		// output for debugging
		if( !quiet ) {
      std::cout << "MUI [comm_mpi.h]: Rank: " << global_rank_ << ", "
                << "Identifier: " << URI << ", "
                << "Domain size: " << local_size_ << ", "
                << "Peers: " << remote_size_
                << std::endl;
		}
		else {
		  if( local_rank_ == 0 ) {
		    std::cout << "MUI [comm_mpi.h]: " << "Identifier: " << URI << ", "
                  << "Domain size: " << local_size_ << ", "
                  << "Peers: " << remote_size_
                  << std::endl;
		  }
		}
		MPI_Barrier( world );

		initialized = true;
	}

	void finalize() {
		if ( initialized && init_by_me ) MPI_Finalize();
	}

	virtual int local_size() const { return local_size_; }
	virtual int local_rank() const { return local_rank_; }
	virtual int remote_size() const { return remote_size_; }
	virtual int global_size() const { return global_size_; }
	virtual int global_rank() const { return global_rank_; }
	virtual std::string uri_host() const { return uri_host_; }
	virtual std::string uri_path() const { return uri_path_; }
	virtual std::string uri_protocol() const { return uri_protocol_; }

protected:
	MPI_Comm domain_local_;
	MPI_Comm domain_remote_;
	int local_size_;
	int local_rank_;
	int remote_size_;
	int global_size_;
	int global_rank_;
	std::string uri_host_;
	std::string uri_path_;
	std::string uri_protocol_;
private:
	bool initialized;
	bool init_by_me;
};

}

#endif
