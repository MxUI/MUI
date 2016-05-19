/*
 * lib_mpi_split.h
 *
 *  Created on: Mar 14, 2014
 *      Author: ytang
 */

#ifndef LIB_MPI_SPLIT_H_
#define LIB_MPI_SPLIT_H_

namespace mui {

inline void mpi_finalize_after_split() {
	int flag;
	MPI_Finalized(&flag);
	if (!flag) MPI_Finalize();
}

inline MPI_Comm mpi_split_by_app( int argc=0, char **argv=NULL )
{
	{
		int flag;
		MPI_Initialized(&flag);
		if( !flag ) {
			MPI_Init( &argc, &argv );
			atexit( mpi_finalize_after_split );
		}
	}
	void* v;
	int flag;
	MPI_Comm_get_attr(MPI_COMM_WORLD,MPI_APPNUM,&v,&flag);
	int appnum = *(int*)v;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	MPI_Comm domain;
	MPI_Comm_split(MPI_COMM_WORLD,appnum,rank,&domain);
	return domain;
}

}

#endif /* LIB_MPI_SPLIT_H_ */
