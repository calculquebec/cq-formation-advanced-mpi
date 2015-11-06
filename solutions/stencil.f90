program stencil

  use mpi_f08

  implicit none

  integer, parameter :: nsources = 3
  integer n, energy, niters, iter, sources(2,nsources)
  double precision heat, t
  character(len=32) :: arg
  double precision, allocatable :: aold(:,:), anew(:,:)

  if(command_argument_count() < 3) then
    write(*,'(A)')'Usage: mpiexec -n 1 ./stencil 240 1 64000 ; display heat.svg &'
    stop
  endif

  call get_command_argument(1, arg) ! nxn grid
  read(arg,*)n
  call get_command_argument(2, arg) ! energy to be injected per iteration
  read(arg,*)energy
  call get_command_argument(3, arg) ! number of iterations
  read(arg,*)niters
  allocate(aold(0:n+1,0:n+1), anew(0:n+1,0:n+1)) ! 1-wide halo zones!
  aold(:,:) = 0d0
  anew(:,:) = 0d0

  call MPI_Init()

  sources = reshape([n/2,n/2,n/3,n/3,n*4/5,n*8/9],shape(sources))
  
  heat=0d0 ! total heat in system
  t=-MPI_Wtime()
  do iter=1,niters
     if (mod(iter,2)==1)then
        call update(aold,anew)
     else
        call update(anew,aold)
     endif
  enddo
  t=t+MPI_Wtime()
  if (mod(iter,2)==1)then
     call printarr(anew, n)
  else
     call printarr(aold, n)
  endif
  write(*,'(A,F12.6,A,F9.6)')"last heat: ",heat," time: ",t

  call MPI_Finalize()

contains

subroutine update(aold,anew)

  double precision, intent(in) :: aold(0:n+1,0:n+1)
  double precision, intent(inout) :: anew(0:n+1,0:n+1)

  integer i,j

  heat = 0d0
  do j=1,n
     do i=1,n
        anew(i,j)=aold(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo
  do i=1,nsources
     anew(sources(1,i),sources(2,i)) = anew(sources(1,i),sources(2,i)) + energy
  enddo
end subroutine update

subroutine printarr(a,n)
  ! does nothing right now, should record each "frame" as image
  integer,intent(in) :: n
  double precision, intent(in) :: a(0:n+1,0:n+1)
  integer, parameter :: size = 1
  integer rgb, i, j

  open(7,FILE='heat.svg',STATUS='unknown',ACCESS='sequential')
  write(7,'(A/A,I0,A,I0,A)')'<?xml version="1.0" encoding="UTF-8"?>','<svg xmlns="http://www.w3.org/2000/svg" version="1.1" width="',size*n,'" height="',size*n,'">'

  write(7,'(A,I0,A,I0,A)')'<rect x="0" y="0" width="',size*n,'" height="',size*n,'" style="stroke-width:1;fill:rgb(0,0,0);stroke:rgb(0,0,0)"/>'
  do i=1,n
     do j=1,n
        if (a(i,j) > 0) then
           rgb = nint(255.0*a(i,j))
        else
           rgb = 0
        endif
        if(rgb>255) rgb=255
        if(rgb/=0) then
           write(7,'(A,I0,A,I0,A,I0,A,I0,A,I0,A,I0,A)') '<rect x="',size*(i-1),'" y="',size*(j-1),'" width="',size,'" height="', &
                size,'" style="stroke-width:1;fill:rgb(',rgb,',0,0);stroke:rgb(',rgb,',0,0)"/>'
        endif
     enddo
  enddo
  write(7,'(A)')'</svg>'

  close(7)
end subroutine printarr

end program stencil
