/*
 * lib_mpi_multidomain.h
 *
 *  Created on: Jun 1, 2015
 *      Author: ytang
 */

#ifndef LIB_MPI_MULTIDOMAIN_H_
#define LIB_MPI_MULTIDOMAIN_H_

#include "lib_mpi_hepler.h"
#include "lib_uri.h"
#include "comm_mpi.h"
#include "uniface.h"

namespace mui
{

using namespace mpi;

template<class CONFIG>
inline std::vector<uniface<CONFIG> *> create_uniface( std::string domain, std::vector<std::string> interfaces, MPI_Comm world = MPI_COMM_WORLD )
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
    for( int i = 0; i < global_size; i++ ) {
        if( global_rank == i ) {
            for( auto &e : map ) {
                printf( "> rank %d \"%s\" registered interface %s as %08X\n", global_rank, domain.c_str(), e.second.c_str(), e.first );
            }
        }
        MPI_Barrier( world );
    }

    int n_unique;
    std::set<int> unique_ifs;
    if( global_rank == 0 ) {
        std::vector<int> nifs = gather( ( int )interfaces.size(), world );
        std::vector<int> displs( global_size + 1, 0 );
        std::partial_sum( nifs.begin(), nifs.end(), displs.begin() + 1 );

        std::vector<int> all_hashes( displs.back() );
        MPI_Gatherv( my_hashes.data(), my_hashes.size(), MPI_INT, all_hashes.data(), nifs.data(), displs.data(), MPI_INT, 0, world );

        for( auto &i : all_hashes ) unique_ifs.insert( i );
        n_unique = unique_ifs.size();
        printf( "> %d distinct interfaces found\n", n_unique );
    } else {
        gather( ( int )interfaces.size(), world );
        MPI_Gatherv( my_hashes.data(), my_hashes.size(), MPI_INT, NULL, NULL, NULL, MPI_INT, 0, world );
    }
    MPI_Bcast( &n_unique, 1, MPI_INT, 0, world );
    std::vector<int> uniq_hashes( n_unique );
    if( global_rank == 0 ) uniq_hashes.assign( unique_ifs.begin(), unique_ifs.end() );
    MPI_Bcast( uniq_hashes.data(), n_unique, MPI_INT, 0, world );

    std::vector<uniface<CONFIG> *> unifaces;
    for( auto &i : uniq_hashes ) {
        MPI_Comm comm_ifs;
        if( map.find( i ) != map.end() ) {
            MPI_Comm_split( world, 1, global_rank, &comm_ifs );
            int comm_rank;
            MPI_Comm_rank( comm_ifs, &comm_rank );
            if( comm_rank == 0 ) printf( "> setting up interface %s [%08X]\n", map[i].c_str(), i );
            std::string full_uri( "mpi://" );
            full_uri = full_uri + domain + "/" + map[i];
            unifaces.push_back( new uniface<CONFIG>( new comm_mpi_smart( full_uri.c_str(), comm_ifs ) ) );
        } else {
            MPI_Comm_split( world, 0, global_rank, &comm_ifs );
        }

        MPI_Barrier( world );
    }

    return unifaces;
}

}

#endif /* LIB_MPI_MULTIDOMAIN_H_ */
