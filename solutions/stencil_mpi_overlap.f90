program stencil

  use mpi_f08

  implicit none

  integer, parameter :: nsources = 3
  integer n, energy, niters, iter, sources(2,nsources)
  integer locnsources ! number of sources in my area
  integer locsources(2,nsources) ! sources local to my rank
  integer locx, locy
  character(len=32) :: arg
  double precision rheat, heat, t
  double precision, allocatable :: aold(:,:), anew(:,:)
  double precision, allocatable :: sbufeast(:), sbufwest(:), rbufeast(:), rbufwest(:)
  integer r, p, args(5), i
  integer px, py, rx, ry, north, south, east, west, bx, by, offx, offy
  type(MPI_Comm) comm

  call MPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, r)
  call MPI_Comm_size(comm, p)

  if (r==0) then
      ! argument checking
      if(command_argument_count() < 5) then
          write(*,'(A)')'usage: stencil_mpi <n> <energy> <niters> <px> <py>'
          call MPI_Finalize()
          stop
      endif
      
      call get_command_argument(1, arg) ! nxn grid
      read(arg,*)n
      call get_command_argument(2, arg) ! energy to be injected per iteration
      read(arg,*)energy
      call get_command_argument(3, arg) ! number of iterations
      read(arg,*)niters
      call get_command_argument(4, arg) ! 1st dim processes
      read(arg,*)px
      call get_command_argument(5, arg) ! 2nd dim processes
      read(arg,*)py
      if(px * py /= p) call MPI_Abort(comm, 1) ! abort if px or py are wrong
      if(mod(n,py) /= 0) call MPI_Abort(comm, 2) ! abort px needs to divide n
      if(mod(n,px) /= 0) call MPI_Abort(comm, 3) ! abort py needs to divide n
      
      ! distribute arguments
      args(:) = (/n, energy, niters, px,  py/)
      call MPI_Bcast(args, 5, MPI_INTEGER, 0, comm)
  else
      call MPI_Bcast(args, 5, MPI_INTEGER, 0, comm)
      n=args(1)
      energy=args(2)
      niters=args(3)
      px=args(4)
      py=args(5)
  endif

  ! determine my coordinates (x,y) -- r=x*a+y in the 2d processor array
  rx = mod(r, px)
  ry = r / px
  ! determine my four neighbors
  north = (ry-1)*px+rx
  if(ry-1 < 0)   north = MPI_PROC_NULL
  south = (ry+1)*px+rx
  if(ry+1 >= py) south = MPI_PROC_NULL
  west= ry*px+rx-1
  if(rx-1 < 0)   west = MPI_PROC_NULL
  east = ry*px+rx+1
  if(rx+1 >= px) east = MPI_PROC_NULL
  ! decompose the domain
  bx = n/px ! block size in x
  by = n/py ! block size in y
  offx = rx*bx ! offset in x
  offy = ry*by ! offset in y

  !printf("%i (%i,%i) - w: %i, e: %i, n: %i, s: %i\n", r, ry,rx,west,east,north,south);

  ! allocate two work arrays
  
  allocate(aold(0:bx+1,0:by+1), anew(0:bx+1,0:by+1)) ! 1-wide halo zones!
  aold(:,:) = 0d0
  anew(:,:) = 0d0

  ! initialize three heat sources
  sources = reshape([n/2,n/2,n/3,n/3,n*4/5,n*8/9],shape(sources))
 
  locnsources=0 ! number of sources in my area
  do i=1,nsources ! determine which sources are in my patch
    locx = sources(1,i) - offx
    locy = sources(2,i) - offy
    if(locx >= 0 .and. locx < bx .and. locy >= 0 .and. locy < by) then
      locnsources = locnsources + 1
      locsources(1,locnsources) = locx+1 ! offset by halo zone
      locsources(2,locnsources) = locy+1 ! offset by halo zone
    endif
  end do

  t=-MPI_Wtime() ! take time
  ! allocate communication buffers
  allocate(sbufeast(by),sbufwest(by)) ! send buffers
  allocate(rbufeast(by),rbufwest(by)) ! receive buffers

  do iter = 1, niters
    if (mod(iter,2)==1)then
       call update(aold,anew)
    else
       call update(anew,aold)
    endif

    ! optional - print image
    !if(iter == niters-1) printarr_par(iter, anew, n, px, py, rx, ry, bx, by, offx, offy, comm);
  enddo

  t=t+MPI_Wtime()
  ! get final heat in the system
  call MPI_Allreduce(heat, rheat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm);
  if(r==0) write(*,'(A,I0,A,F12.6,A,F9.6)')"[",r,"] last heat: ",rheat," time: ",t

  call MPI_Finalize()

contains

subroutine update(aold,anew)

  double precision, intent(inout) :: aold(0:bx+1,0:by+1)
  double precision, intent(inout) :: anew(0:bx+1,0:by+1)

  integer i,j
  type(MPI_Request) reqs(8)
  type(MPI_Status) status(8)

  ! refresh heat sources
  do i=1,locnsources
     aold(locsources(1,i),locsources(2,i)) = aold(locsources(1,i),locsources(2,i)) + energy ! heat source
  enddo

  ! exchange data with neighbors
  sbufeast(:) = aold(bx,1:by) ! pack
  sbufwest(:) = aold(1,1:by)
  call MPI_Isend(aold(1,1), bx, MPI_DOUBLE_PRECISION, north, 9, comm, reqs(1))
  call MPI_Isend(aold(1,by), bx, MPI_DOUBLE_PRECISION, south, 9, comm, reqs(2))
  call MPI_Isend(sbufeast, by, MPI_DOUBLE_PRECISION, east, 9, comm, reqs(3))
  call MPI_Isend(sbufwest, by, MPI_DOUBLE_PRECISION, west, 9, comm, reqs(4))
  call MPI_Irecv(aold(1,0), bx, MPI_DOUBLE_PRECISION, north, 9, comm, reqs(5))
  call MPI_Irecv(aold(1,by+1), bx, MPI_DOUBLE_PRECISION, south, 9, comm, reqs(6))
  call MPI_Irecv(rbufeast, by, MPI_DOUBLE_PRECISION, east, 9, comm, reqs(7))
  call MPI_Irecv(rbufwest, by, MPI_DOUBLE_PRECISION, west, 9, comm, reqs(8))

  ! update inner grid points
  heat=0d0 ! total heat in system
  do j=2,by-1
     do i=2,bx-1
        anew(i,j)=anew(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo

  call MPI_Waitall(8, reqs, status)

  ! update outer grid points (two elements less per loop to avoid double update)
  do j=1,by,by-1
     do i=2,bx-1
        anew(i,j)=anew(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo
  do j=1,by
     anew(1,j)=anew(1,j)/2 + (rbufwest(j)+aold(2,j)+aold(1,j-1)+aold(1,j+1))/8
     heat = heat + anew(1,j)
     anew(bx,j)=anew(bx,j)/2 + (aold(bx-1,j)+rbufeast(j)+aold(bx,j-1)+aold(bx,j+1))/8
     heat = heat + anew(bx,j)
  enddo

end subroutine update

end program stencil
