!Author: Dante Buhl
! CFD HW 4

program advect
  use LinAl
  use NumDE
  use CFD

  implicit none

  ! Global vars
  integer :: i, j, fn

  ! Simvars
  real(kind=kr) :: sf = 0.8
  real(kind=kr) :: xstart, xstop, tstart, tstop, dx, dt
  integer :: Nx, Nt, ic
  real(kind=kr), dimension(6) :: cfl_array, nx_array
  real(kind=kr), allocatable :: X(:), U(:,:), cell_bound(:), U_grid(:,:)

  ! Loop Vars
  integer :: counter = 0

  print *, " ----------------------------------------------------------- "
  print *, " "
  
  xstart = 0
  xstop = 1
  Nx = 32
  Nt = 100
  allocate(U(Nx+2,Nt), U_grid(Nx+1,Nt))

  ! this specifies which IC to use (1-11)
  ! 1 corresponds to IC 7.38
  ! 11 corresponds to IC 7.48
  ic = 1 
  
  call grid_init(xstart, xstop, nx, x, cell_bound, dx)
  call advect_init(U(:,1),U_grid(:,1),cell_bound,nx,dx,ic)
  call bc(U(:,1),nx,'p')
  call cfl(U(:,1),dx,sf,dt)

  do i = 1, Nt-1
    U(:,i+1) = advect_update(U(:,i),nx,dx,dt,'f','u')
    call bc(U(:,i),nx,'p')
    call cfl(U(:,i),dx,sf,dt)
  end do

  print *, " "
  print *, " ----------------------------------------------------------- "

  contains
    include 'advect_update.f90'  
end program advect
