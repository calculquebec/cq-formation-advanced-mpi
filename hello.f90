PROGRAM hello
USE omp_lib
implicit none
INCLUDE 'mpif.h' ! or USE mpi / USE mpi_f08 (MPI 3.0)

INTEGER rank, size, provided, ierr
CALL MPI_Init_Thread(MPI_THREAD_FUNNELED,&
                     provided, ierr)
CALL MPI_Comm_rank (MPI_COMM_WORLD, &
                     rank, ierr)
CALL MPI_Comm_size (MPI_COMM_WORLD, &
                     size, ierr)
!$OMP PARALLEL
   WRITE(*,*) 'Hello from processor ',&
     rank, ' of ', size, ' thread ',&
   omp_get_thread_num()
!$OMP END PARALLEL
CALL MPI_Finalize(ierr)

END PROGRAM hello

