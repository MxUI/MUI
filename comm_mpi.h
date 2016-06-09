/* -*- c++ -*-
 * comm_mpi.h
 *
 *  Created on: Feb 11, 2014
 *      Author: skudo
 */

#ifndef COMM_MPI_H_
#define COMM_MPI_H_

#include <mpi.h>
#include "comm.h"
#include "util.h"
#include "lib_uri.h"

namespace mui {

/*
 * MPI base implementation of Comm
 */

class comm_mpi : public communicator {
public:
	comm_mpi(const char URI[], MPI_Comm world ) :
		domain_local_(0), domain_remote_(0),
		local_size_(0), local_rank_(0), remote_size_(0), global_size_(0), global_rank_(0),
		initialized(false), init_by_me(false) {
		init(URI, world);
	}
	virtual ~comm_mpi() {
		finalize();
	}

	void init( const char URI[], MPI_Comm world ) {
		if ( initialized ) {
			throw( std::runtime_error("duplicate MUI communicator initialization") );
		}

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
		if (global_size_ < 2) {
			fprintf( stderr, "WARNING: # of global ranks less than 2\n" );
		}

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
		int domain_hash = std::hash<std::string>()( desc.host() );
		domain_hash=std::abs(domain_hash);
		int ifs_hash    = std::hash<std::string>()( desc.path() );
		ifs_hash=std::abs(ifs_hash);
		MPI_Comm_split( world, domain_hash % prime, global_rank_, &domain_local_ );
		MPI_Comm_size( domain_local_, &local_size_ );
		MPI_Comm_rank( domain_local_, &local_rank_ );

		// set intercommunication using interface tag
		std::vector<int> all_domain( global_size_, 0 );
		std::vector<int> all_ifs   ( global_size_, 0 );
		MPI_Allgather( &domain_hash, 1, MPI_INT, &all_domain[0], 1, MPI_INT, world );
		MPI_Allgather( &ifs_hash,    1, MPI_INT, &all_ifs[0],    1, MPI_INT, world );
		for(int i = 0 ; i < global_size_ ; i++) {
			if ( i == global_rank_ ) continue;
			if ( all_domain[i] != domain_hash && all_ifs[i] == ifs_hash ) {
				MPI_Intercomm_create( domain_local_, 0, world, i, ifs_hash % prime, &domain_remote_ );
				break;
			}
		}
		MPI_Comm_remote_size( domain_remote_, &remote_size_ );

		// output for debugging
		for(int i=0;i<global_size_;i++) {
			if ( i == global_rank_ )
				std::cout	<<"rank "<<global_rank_<<'\t'
							<<"identifier "<<URI<<'\t'
							<<"domain size "<<local_size_<<'\t'
							<<"peer number "<<remote_size_<<'\t'
							<<std::endl;
			MPI_Barrier( world );
		}

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

protected:
	MPI_Comm domain_local_;
	MPI_Comm domain_remote_;
	int local_size_;
	int local_rank_;
	int remote_size_;
	int global_size_;
	int global_rank_;
private:
	bool initialized, init_by_me;
};

}

#endif
