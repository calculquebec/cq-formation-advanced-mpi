program stencil_omp

  use mpi_f08
  use omp_lib

  implicit none

  integer, parameter :: nsources = 3
  integer n, energy, niters, iter, sources(2,nsources), j
  double precision heat, t
  character(len=32) :: arg
  double precision, allocatable :: aold(:,:), anew(:,:)

  call get_command_argument(1, arg) ! nxn grid
  read(arg,*)n
  call get_command_argument(2, arg) ! energy to be injected per iteration
  read(arg,*)energy
  call get_command_argument(3, arg) ! number of iterations
  read(arg,*)niters
  allocate(aold(0:n+1,0:n+1), anew(0:n+1,0:n+1)) ! 1-wide halo zones!

  call MPI_Init()

  aold(:,0)=0d0
  anew(:,0)=0d0
!$omp parallel
!$omp do
  do j=1,n
     aold(:,j)=0d0
     anew(:,j)=0d0
  enddo
!$omp end do
  if (omp_get_thread_num()==omp_get_num_threads()-1)then
     aold(:,n+1)=0d0
     anew(:,n+1)=0d0
  endif
!$omp end parallel

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

!$omp parallel do reduction(+:heat)
  do j=1,n
     do i=1,n
        anew(i,j)=anew(i,j)/2 + (aold(i-1,j)+aold(i+1,j)+aold(i,j-1)+aold(i,j+1))/8
        heat = heat + anew(i,j)
     enddo
  enddo
!$omp end parallel do
  do i=1,nsources
     anew(sources(1,i),sources(2,i)) = anew(sources(1,i),sources(2,i)) + energy
  enddo
end subroutine update

subroutine printarr(a,n)
  ! does nothing right now, should record each "frame" as image
  integer,intent(in) :: n
  double precision, intent(in) :: a(0:n+1,0:n+1)
  integer, parameter :: size = 5
  integer rgb, i, j

  open(7,FILE='heat.svg',STATUS='unknown',ACCESS='sequential')
  write(7,'(A/A/A)')'<html>','<body>','<svg xmlns="http://www.w3.org/2000/svg" version="1.1">'

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
  write(7,'(A/A/A)')'</svg>','</body>','</html>'

  close(7)
end subroutine printarr

end program stencil_omp
