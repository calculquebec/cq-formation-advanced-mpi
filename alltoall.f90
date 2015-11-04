program alltoall

  use mpi
  implicit none

  integer send(4),recv(3)

  integer rank,size,ierr,k
	
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, size, ierr)
  if(size/=4)then 
     write(*,*)'Error!:# of processors must be equal to 4'
     write(*,*)'Program aborting....'
     call MPI_ABORT(ierr,1)
  endif

  do k=1,size   
     send(k) = k + rank*size
  enddo

  write(*,*))rank, ': ', 'send = ',send

  call MPI_ALLTOALL(send, 2, MPI_REAL, recv, 1, &
       MPI_INTEGER, MPI_COMM_WORLD, ierr)

  write(*,*)rank,': ','recv = ',recv

  call MPI_FINALIZE(ierr)
end program alltoall
