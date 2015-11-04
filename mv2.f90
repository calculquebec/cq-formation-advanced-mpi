! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Data distribution of matrix: columnwise block striped
! based on example by
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6495
! 7 January 2009
!=================================================================


program MatVec_MPI
  !-----------------------------------------------------------------
  ! Description: 
  !    Computes the product of an m*n matrix and an n-vector.
  !-----------------------------------------------------------------
  !======= Inclusions ===========
  use mpi

  !======= Declarations =========
  implicit none

  integer :: m, n, i, is, ie, j, js, je, numprocs, myid, ierr, repeat
  integer, allocatable :: counts(:)
  double precision, allocatable :: A(:,:), x(:), myv(:), v(:)
  double precision :: stime, ftime, mynorm2, norm2

  !======= Internals ============
  
  ! initialize MPI
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

  ! input the size of the system
  if (myid == 0) then
     write(*,*) 'We will multiply an m*n matrix by an n-vector'
     read(*,*) m, n
     print *, 'm = ',m,', n = ',n
  endif

  ! root node sends m and n out to other processors
  call mpi_bcast(m, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
  call mpi_bcast(n, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

  if ((m < 1) .or. (n < 1)) then
     if (myid == 0) then
        write(*,*)  ' illegal input, m =',m,' and n =',n,&
             ' must both be positive'
     endif
     call mpi_abort(MPI_COMM_WORLD, 1, ierr)
     stop
  endif
  
  ! root node outputs parallelism information to screen
  if (myid == 0) then
     write(*,'(A,i5,A)') ' starting MPI with ', numprocs,' processes'
  endif

  ! determine this processor's interval
  is = myid*m/numprocs+1
  ie = (myid+1)*m/numprocs
  js = myid*n/numprocs+1
  je = (myid+1)*n/numprocs

  ! initialize the matrix and vectors 
  ! (only store matrix & vector parts that you need on this proc)
  allocate(A(m,js:je),x(js:je),myv(m),v(is:ie))
  do j=js,je
     do i=1,m
        A(i,j) = 1.d0/(1+(i-j)**2)
     enddo
     x(j) = 1.d0
  enddo

  allocate(counts(numprocs))
  do i=1,numprocs
     counts(i) = i*m/numprocs - (i-1)*m/numprocs
  enddo

  ! start timer
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  stime = mpi_wtime()

  do repeat=1,10000

     ! compute matrix-vector product
     do i=1,m
        myv(i) = 0.d0
     enddo

     do j=js,je
        do i=1,m
           myv(i) = myv(i) + A(i,j)*x(j)
        enddo
     enddo

     call mpi_reduce_scatter(...)
  enddo

  ! stop timer
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  ftime = mpi_wtime()

  mynorm2 = dot_product(v, v)
  call mpi_reduce(mynorm2, norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)

  ! output 2-norm of product and runtime to screen
  if (myid == 0) then
     write(*,*) ' 2-norm of product =',sqrt(norm2)
     write(*,*) ' runtime =',ftime-stime
  endif

  ! finalize MPI
  call mpi_finalize(ierr)  

  ! free vectors
  deallocate(A,x,v,counts)

end program MatVec_MPI
!=================================================================
