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
 * @file lib_mpi_helper.h
 * @author Y. H. Tang
 * @date 1 June 2015
 * @brief MPI data types used internally by MUI.
 */

#ifndef LIB_MPI_HELPER_H_
#define LIB_MPI_HELPER_H_

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

#endif /* LIB_MPI_HELPER_H_ */
