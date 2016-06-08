#include <stdio.h>
#include <mpi.h>
#include <omp.h>

int main (int argc, char * argv[])
{
  int rank, size, provided;
  MPI_Init_thread(&argc, &argv,
		  MPI_THREAD_FUNNELED, &provided);
  MPI_Comm_rank( MPI_COMM_WORLD, &rank );
  MPI_Comm_size( MPI_COMM_WORLD, &size );
#pragma omp parallel
  {
    printf( "Hello from processor %d"
	    " of %d, thread %d\n", rank, size,
	    omp_get_thread_num());
  }
  MPI_Finalize();
  return 0;
}
