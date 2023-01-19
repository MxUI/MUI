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
 * @file lib_mpi_multidomain.h
 * @author Y. H. Tang
 * @date 1 June 2015
 * @brief Provides helper functions for creating and synchronising multiple
 * MUI interfaces for a single domain.
 */

#ifndef LIB_MPI_MULTIDOMAIN_H_
#define LIB_MPI_MULTIDOMAIN_H_

#include "lib_uri.h"
#include "comm_mpi.h"
#include "lib_mpi_helper.h"
#include "../uniface.h"

namespace mui
{

using namespace mpi;

template<class CONFIG>
inline std::vector<std::unique_ptr<uniface<CONFIG>>> create_uniface( std::string domain, std::vector<std::string> interfaces, MPI_Comm world = MPI_COMM_WORLD )
{
    // gather total number of interfaces (edges of the graph)
    int global_size, global_rank;
    MPI_Comm_size( world, &global_size );
    MPI_Comm_rank( world, &global_rank );

    std::map<int, std::string> map;

    std::vector<int> my_hashes;
    for( auto &i : interfaces ) {
        auto h = std::hash<std::string>()( i );
        my_hashes.push_back( h );
        map[h] = i;
    }

    // output for debugging
    if( !CONFIG::QUIET ) {
      for( auto &e : map ) {
        std::cout << "MUI [lib_mpi_multidomain]: Rank: " << global_rank << ", \"" << domain
        << "\" registered interface \"" << e.second
        << "\" as " << std::hex << e.first << std::dec << std::endl;
      }
    }
    else {
      if( global_rank == 0 ) {
        for( auto &e : map ) {
          std::cout << "MUI [lib_mpi_multidomain]: \"" << domain
          << "\" registered interface \"" << e.second
          << "\" as " << std::hex << e.first << std::dec << std::endl;
        }
      }
    }
    MPI_Barrier( world );


    int n_unique;
    std::set<int> unique_ifs;
    if( global_rank == 0 ) {
        std::vector<int> nifs = gather( static_cast<int>(interfaces.size()), world );
        std::vector<int> displs( global_size + 1, 0 );
        std::partial_sum( nifs.begin(), nifs.end(), displs.begin() + 1 );

        std::vector<int> all_hashes( displs.back() );
        MPI_Gatherv( my_hashes.data(), my_hashes.size(), MPI_INT, all_hashes.data(), nifs.data(), displs.data(), MPI_INT, 0, world );

        for( auto &i : all_hashes ) unique_ifs.insert( i );
        n_unique = unique_ifs.size();
        std::cout << "MUI Info [lib_mpi_multidomain]: " << n_unique << " distinct interface(s) found" << std::endl;
    } else {
        gather( static_cast<int>(interfaces.size()), world );
        MPI_Gatherv( my_hashes.data(), my_hashes.size(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, world );
    }

    MPI_Barrier( world );
    MPI_Bcast( &n_unique, 1, MPI_INT, 0, world );
    std::vector<int> uniq_hashes( n_unique );
    if( global_rank == 0 ) uniq_hashes.assign( unique_ifs.begin(), unique_ifs.end() );
    MPI_Bcast( uniq_hashes.data(), n_unique, MPI_INT, 0, world );

    std::vector<uniface<CONFIG>*> unifaces;
    for( auto &i : uniq_hashes ) {
        MPI_Comm comm_ifs;
        if( map.find( i ) != map.end() ) {
            MPI_Comm_split( world, 1, global_rank, &comm_ifs );
            int comm_rank;
            MPI_Comm_rank( comm_ifs, &comm_rank );
            if( comm_rank == 0 ) {
            	std::cout << "MUI [lib_mpi_multidomain]: Setting up interface " << map[i] << " [" << std::hex << i
            				<< std::dec << "] (rank ids are local to each interface)" << std::endl;
            }
            std::string full_uri( "mpi://" );
            full_uri = full_uri + domain + "/" + map[i];
            unifaces.push_back( new uniface<CONFIG>( new comm_mpi_smart( full_uri.c_str(), CONFIG::QUIET, comm_ifs ) ) );
        } else {
            MPI_Comm_split( world, 0, global_rank, &comm_ifs );
        }

        MPI_Barrier( world );
    }

    // Sort return vector into original interface order (order can be mangled due to use of hash values)
	std::vector<std::unique_ptr<uniface<CONFIG>>> unifaces_sorted;

	for (const auto &orig_inter : interfaces) {
		for (auto &uni : unifaces) {
			if( uni->uri_path().compare(orig_inter) == 0 ) {
				unifaces_sorted.push_back( std::unique_ptr<uniface<CONFIG>>(uni) );
				break;
			}
		}
	}

	return unifaces_sorted;
}

  // Ensure all interfaces in the range of the two iterators have hit sync_all
  template<typename IteratorT>
  void sync_all(IteratorT begin, IteratorT end,
		typename std::decay<decltype(**begin)>::type::time_type t) {
    for(auto it = begin; it != end; ++it) {
      auto& unif = *it;
      // If nothing pushed this just announces the current time 
      unif->commit(t);
    }

    for(auto it = begin; it != end; ++it) {
      auto& unif = *it;
      // This waits for all the other end of the interface to announce
      // a time greater than or equal to t
      unif->barrier(t);
    }
  }
}

#endif /* LIB_MPI_MULTIDOMAIN_H_ */
