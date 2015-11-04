#include <stdio.h>
#include "mpi.h"

int main( int argc, char **argv )
{
    int send[4], recv[4];
    int rank, size, k;

    MPI_Init( &argc, &argv );
    MPI_Comm_rank( MPI_COMM_WORLD, &rank );
    MPI_Comm_size( MPI_COMM_WORLD, &size );

    if (size != 4) {
       printf("Error!:# of processors must be equal to 4");
       printf("Programm aborting....");
       MPI_Abort(MPI_COMM_WORLD, 1);
    }
    for (k=0;k<size;k++) send[k] = (k+1) + rank*size;

    printf("%d : send = %d %d %d %d\n", rank, send[0], send[1], send[2], send[3]);

    MPI_Alltoall(send, 1, MPI_INT, recv, 1, MPI_INT, MPI_COMM_WORLD);

    printf("%d : recv = %d %d %d %d\n", rank, recv[0], recv[1], recv[2], recv[3]);

    MPI_Finalize();
    return 0;
}
