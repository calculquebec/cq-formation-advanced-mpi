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
  type(MPI_Datatype) north_south_type, east_west_type
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

  ! create north-south datatype
  call MPI_Type_contiguous(bx, MPI_DOUBLE_PRECISION, north_south_type)
  call MPI_Type_commit(north_south_type)
  ! create east-west type
  call MPI_Type_vector(by,1,bx+2,MPI_DOUBLE_PRECISION, east_west_type)
  call MPI_Type_commit(east_west_type)

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

  call MPI_Type_free(east_west_type)
  call MPI_Type_free(north_south_type)

  ! get final heat in the system
  call MPI_Allreduce(heat, rheat, 1, MPI_DOUBLE_PRECISION, MPI_SUM, comm);
  if(r==0) write(*,'(A,I0,A,F12.6,A,F9.6)')"[",r,"] last heat: ",rheat," time: ",t

  call MPI_Finalize()

contains

subroutine update(aold,anew)

  double precision, intent(inout) :: aold(0:bx+1,0:by+1)
  double precision, intent(inout) :: anew(0:bx+1,0:by+1)

  integer i,j
  integer counts(4), sz
  integer(MPI_Address_kind) displs(4), rdispls(4), a1, a2
  type(MPI_Datatype) types(4)
  
  ! refresh heat sources
  do i=1,locnsources
     aold(locsources(1,i),locsources(2,i)) = aold(locsources(1,i),locsources(2,i)) + energy ! heat source
  enddo

  ! exchange data with neighbors
  call MPI_Get_Address(aold(0,0),a1)
  call MPI_Get_Address(aold(1,0),a2)
  sz = a2 - a1
  counts = [1, 1, 1, 1]
  displs = [0, (bx-1)*sz, 0, (bx+2)*(by-1)*sz]
  types = [east_west_type, east_west_type, north_south_type, north_south_type]
  rdispls = [(bx+2)*sz, (bx+2)*sz+(bx+1)*sz, sz, (bx+2)*(by+1)*sz+sz]

  call MPI_Neighbor_alltoallw(aold(1,1), counts, displs, types, aold(0,0), counts, rdispls, types, topocomm)

  ! update grid points
  heat=0d0 ! total heat in system
  do j=1,by
     do i=1,bx
        anew(i,j)=aold(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo
end subroutine update

end program stencil
