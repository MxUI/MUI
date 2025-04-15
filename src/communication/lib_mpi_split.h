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
 * @file lib_mpi_split.h
 * @author Y. H. Tang
 * @date 14 March 2014
 * @brief Provides helper functions to generate (and finalize) a new MPI
 * comm world that can then be used by an already MPI enabled application
 * to ensure MPI_COMM_WORLD remains free for use by MUI.
 */

#ifndef LIB_MPI_SPLIT_H_
#define LIB_MPI_SPLIT_H_

#include <exception>
namespace mui {

  inline void mpi_finalize_after_split() {
    int flag;
    MPI_Finalized(&flag);
    if (!flag) MPI_Finalize();
  }

  inline MPI_Comm mpi_split_by_app( int argc=0, char **argv=NULL, int threadType=-1, int *thread_support=NULL, bool use_mpi_comm_split = true )
  {
    {
    int flag;
    MPI_Initialized(&flag);
    if( !flag ) {
        if(threadType != -1) {
          MPI_Init_thread( &argc, &argv, threadType, thread_support );
        }
        else {
          MPI_Init( &argc, &argv );
        }
        atexit( mpi_finalize_after_split );
      }
    }
    void* v;
    int flag;
    MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_APPNUM,&v,&flag);
    if (!flag) {
        std::cout << "MUI Info [lib_mpi_split.h]: Calling mpi_split_by_app() when run as a single application" << std::endl;
    }
    int appnum = *static_cast<int*>(v);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);
    MPI_Comm domain;
    if (use_mpi_comm_split) {
        MPI_Comm_split(MPI_COMM_WORLD,appnum,rank,&domain);
    } else {
        int size;
        MPI_Comm_size(MPI_COMM_WORLD,&size);
        std::vector<int> all_appnums(size);
        all_appnums[rank] = appnum;
        MPI_Allgather(MPI_IN_PLACE, 0, MPI_INT, all_appnums.data(), 1, MPI_INT, MPI_COMM_WORLD);
        std::vector<int> group_ranks;
        for (int i = 0; i < size; ++i) {
            if (all_appnums[i] == appnum) {
                group_ranks.push_back(i);
            }
        }
        MPI_Group world_group, new_group;
        MPI_Comm_group(MPI_COMM_WORLD, &world_group);
        MPI_Group_incl(world_group, group_ranks.size(), group_ranks.data(), &new_group);
        MPI_Comm_create_group(MPI_COMM_WORLD, new_group, 0, &domain);
        MPI_Group_free(&world_group);
        MPI_Group_free(&new_group);
    }

    return domain;
  }
}

#endif /* LIB_MPI_SPLIT_H_ */
