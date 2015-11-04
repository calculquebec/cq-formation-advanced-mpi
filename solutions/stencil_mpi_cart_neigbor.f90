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
  integer r, p, args(3), i
  integer px, py, rx, ry, north, south, bx, by, offx, offy
  type(MPI_Comm) comm, topocomm
  double precision, allocatable :: sbuf(:),rbuf(:)
  integer pdims(2), coords(2)
  logical periods(2)

  call MPI_Init()

  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, r)
  call MPI_Comm_size(comm, p)

  if (r==0) then
      ! argument checking
      if(command_argument_count() < 3) then
          write(*,'(A)')'usage: stencil_mpi <n> <energy> <niters>'
          call MPI_Finalize()
          stop
      endif
      
      call get_command_argument(1, arg) ! nxn grid
      read(arg,*)n
      call get_command_argument(2, arg) ! energy to be injected per iteration
      read(arg,*)energy
      call get_command_argument(3, arg) ! number of iterations
      read(arg,*)niters
      
      ! distribute arguments
      args(:) = (/n, energy, niters/)
      call MPI_Bcast(args, 3, MPI_INTEGER, 0, comm)
  else
      call MPI_Bcast(args, 3, MPI_INTEGER, 0, comm)
      n=args(1)
      energy=args(2)
      niters=args(3)
  endif

  pdims=0
  ! compute good (rectangular) domain decomposition
  call MPI_Dims_create(p, 2, pdims)
  px = pdims(1)
  py = pdims(2)

  ! create Cartesian topology
  periods = .false.
  call MPI_Cart_create(comm, 2, pdims, periods, .false., topocomm)

  ! get my local x,y coordinates
  call MPI_Cart_coords(topocomm, r, 2, coords)
  rx = coords(1)
  ry = coords(2)

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

  allocate(sbuf(2*bx+2*by),rbuf(2*bx+2*by)) ! send/recv buffer (west, east, north, south)

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
  integer counts(4), displs(4)
  
  ! refresh heat sources
  do i=1,locnsources
     aold(locsources(1,i),locsources(2,i)) = aold(locsources(1,i),locsources(2,i)) + energy ! heat source
  enddo

  ! exchange data with neighbors
  sbuf(1:by) = aold(1,1:by) ! pack west
  sbuf(by+1:by+by) = aold(bx,1:by) ! pack east
  sbuf(2*by+1:2*by+bx) = aold(1:bx,1) ! pack north
  sbuf(2*by+bx+1:2*by+2*bx) = aold(1:bx,by) ! pack south
  counts = [by, by, bx, bx]
  displs = [0, by, 2*by, 2*by+bx]

  call MPI_Neighbor_alltoallv(sbuf, counts, displs, MPI_DOUBLE_PRECISION, rbuf, counts, displs, MPI_DOUBLE_PRECISION, topocomm)

  aold(0,1:by) = rbuf(1:by)
  aold(bx+1,1:by) = rbuf(by+1:by+by)
  aold(1:bx,0) = rbuf(2*by+1:2*by+bx)
  aold(1:bx,by+1) = rbuf(2*by+bx+1:2*by+2*bx)

  ! update grid points
  heat=0d0 ! total heat in system
  do j=1,by
     do i=1,bx
        anew(i,j)=anew(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo
end subroutine update

end program stencil
