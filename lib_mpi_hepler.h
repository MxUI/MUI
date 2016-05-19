/*
 * lib_mpi_hepler.h
 *
 *  Created on: Jun 1, 2015
 *      Author: ytang
 */

#ifndef LIB_MPI_HEPLER_H_
#define LIB_MPI_HEPLER_H_

#include <mpi.h>

namespace mui {

namespace mpi {

template<typename T>
inline MPI_Datatype mpi_type( T              const& t );
inline MPI_Datatype mpi_type( int            const &t ) { return MPI_INT; }
inline MPI_Datatype mpi_type( long           const &t ) { return MPI_LONG; }
inline MPI_Datatype mpi_type( unsigned long  const &t ) { return MPI_UNSIGNED_LONG; }
inline MPI_Datatype mpi_type( long long      const &t ) { return MPI_LONG_LONG; }
inline MPI_Datatype mpi_type( float          const &t ) { return MPI_FLOAT; }
inline MPI_Datatype mpi_type( double         const &t ) { return MPI_DOUBLE; }
inline MPI_Datatype mpi_type( char           const &t ) { return MPI_CHAR; }
inline MPI_Datatype mpi_type( short          const &t ) { return MPI_SHORT; }
inline MPI_Datatype mpi_type( unsigned short const &t ) { return MPI_UNSIGNED_SHORT; }

template<typename T> inline std::vector<T> gather( T t, MPI_Comm comm ) {
	int size, rank;
	MPI_Comm_size( comm, &size );
	MPI_Comm_rank( comm, &rank );
	std::vector<T> v;
	if ( rank == 0 ) v.resize( size );
	MPI_Gather( &t, 1, mpi_type(t), v.data(), 1, mpi_type(t), 0, comm );
	return v;
}

}

}

#endif /* LIB_MPI_HEPLER_H_ */
