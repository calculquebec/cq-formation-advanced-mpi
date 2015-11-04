! -*- Mode: Fortran90; -*-
!-----------------------------------------------------------------
! Data distribution of matrix: checkerboard block
! based on example by
! Daniel R. Reynolds
! SMU, Mathematics
! Math 6495
! 7 January 2009
! and Michael J. Quinn, Parallel programming in C with MPI and OpenMP
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
  integer :: col_comm, row_comm, grid_comm, grid_id, status(MPI_STATUS_SIZE)
  integer :: grid_size(2), grid_coords(2)
  logical :: periodic(2)
  integer :: src, dest, coords(2)
  double precision, allocatable :: A(:,:), x(:), myv(:), v(:)
  double precision :: stime, ftime, mynorm2, norm2

  !======= Internals ============
  
  ! initialize MPI
  call mpi_init(ierr)
  call mpi_comm_size(MPI_COMM_WORLD, numprocs, ierr)
  call mpi_comm_rank(MPI_COMM_WORLD, myid, ierr)

  grid_size(:) = int(sqrt(dble(numprocs)))
  if (grid_size(1)*grid_size(2) /= numprocs) then
     write(*,*) 'Number of processors ', numprocs, ' not square'
     call mpi_abort(MPI_COMM_WORLD, 1, ierr)
     stop
  endif
    
  ! input the size of the system
  if (myid == 0) then
     write(*,*) 'We will multiply an m*n matrix by an n-vector'
     read(*,*) m, n
     print *, 'm = ',m,', n = ',n
  endif

  periodic(:) = .false.
  call mpi_cart_create(MPI_COMM_WORLD, 2, grid_size, periodic, .true., grid_comm, &
       ierr)
  call mpi_comm_rank(grid_comm, grid_id, ierr)
  call mpi_cart_coords(grid_comm, grid_id, 2, grid_coords, ierr)
  call mpi_comm_split(grid_comm, grid_coords(1), grid_coords(2), row_comm, ierr)
  call mpi_comm_split(grid_comm, grid_coords(2), grid_coords(1), col_comm, ierr)

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
  is = grid_coords(1)*m/grid_size(1)+1
  ie = (grid_coords(1)+1)*m/grid_size(1)
  js = grid_coords(2)*n/grid_size(2)+1
  je = (grid_coords(2)+1)*n/grid_size(2)

  ! initialize the matrix and vectors 
  ! (only store matrix & vector parts that you need on this proc)
  allocate(A(is:ie,js:je),x(js:je),myv(is:ie))
  do j=js,je
     do i=is,ie
        A(i,j) = 1.d0/(1+(i-j)**2)
     enddo
  enddo

  ! initialize x and v for the first column only
  if (grid_coords(2) == 0) then
     do j=js,je
        x(j) = 1.0
     enddo
     allocate(v(is:ie))
  endif

  ! start timer
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  stime = mpi_wtime()

  do repeat=1,10000

     ! Step 1: send elements from first column processes to first row processes
     if (grid_coords(2) == 0 .and. grid_coords(1) > 0) then
        coords(1) = 0
        coords(2) = grid_coords(1)
        call mpi_cart_rank(grid_comm, coords, dest, ierr)
        call mpi_send(x, je-js+1, MPI_DOUBLE_PRECISION, dest, 0, grid_comm, ierr);
     endif

     if (grid_coords(1) == 0 .and. grid_coords(2) > 0) then
        coords(1) = grid_coords(2)
        coords(2) = 0
        call mpi_cart_rank(grid_comm, coords, src, ierr)
        call mpi_recv(x, je-js+1, MPI_DOUBLE_PRECISION, src, MPI_ANY_TAG, grid_comm, status, ierr)
     endif

     ! Step 2: Broadcast along columns
     call mpi_bcast(...)

     ! compute matrix-vector product
     do i=is,ie
        myv(i) = 0.d0
     enddo

     do j=js,je
        do i=is,ie
           myv(i) = myv(i) + A(i,j)*x(j)
        enddo
     enddo

     call mpi_reduce(...)
  enddo

  ! stop timer
  call mpi_barrier(MPI_COMM_WORLD, ierr)
  ftime = mpi_wtime()

  ! output 2-norm of product and runtime to screen
  if (grid_coords(2) == 0) then
     mynorm2 = dot_product(v, v)
     call mpi_reduce(mynorm2, norm2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, col_comm, ierr)
     if (grid_coords(1) == 0) then
        write(*,'(A,I0,A,I0,A,I0,A,F14.6)') ' rank=',myid,' = (',grid_coords(1),',',grid_coords(2), &
             '), 2-norm of product =',sqrt(norm2)
     endif
  endif

  if (myid == 0) then
     write(*,*) ' runtime =',ftime-stime
  endif

  ! finalize MPI
  call mpi_finalize(ierr)  

  ! free vectors
  deallocate(A,x)

end program MatVec_MPI
!=================================================================
