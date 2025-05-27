!File: sim_time_est.f90
!Author: Dante Buhl
!Purpose: to estimate the number of timesteps needed to obtain simulations in a
!statistically stationary state. 

program driver

  use LinAl
  use NumDE
  use CFD

  implicit none

  ! Global vars
  !integer, parameter :: kr=kind(dble(1.0))
  integer :: i, j, fn

  ! Simvars
  real(kind=kr) :: cfl = 0.8
  real(kind=kr) :: xstart, xstop, tstart, tstop, dx, dt
  real(kind=kr) :: a = 1
  integer :: Nx, Nt
  real(kind=kr), dimension(6) :: cfl_array, nx_array
  real(kind=kr), allocatable :: x(:,:), u(:,:) 

  ! Loop Vars
  integer :: counter = 0

  !----------------------------------------------------------------------------
  print *, " ----------------------------------------------------------- "
  print *, " "
  print *, " Problem 6:"
  print *, " "
   
    fn = 20 
    open(fn, file='prob6.dat')
    cfl_array = (/0.8, 0.8, 1.0, 1.0, 1.2, 1.2/)
    nx_array = (/32, 128, 32, 128, 32, 128/)
    xstart = 0
    xstop = 1 
    tstart = 0
    tstop = 1 

    ! initializing data
    do j = 1, 6
      cfl = cfl_array(j)
      a = 1
      Nx = nx_array(j)
      dx = (xstop-xstart)/(Nx-1)
      dt = dx*cfl/a
      Nt = ceiling((tstop-tstart)/dt)
      allocate(x(Nx,Nt), u(Nx,Nt))
      do i = 1, Nx
        x(i,1) = xstart + (i-1)*dx
      end do  
      u = sin(2*pi*x) ! might have to average these together over the cells

      ! running timestepper
      do i = 1, Nt-1
        u(:,i+1) = LF_update_1D(u(:,i),u(:,i),dt,dx,Nx)
      end do 
      write(fn,*) Nx, Nt
      call writemat(u, Nx, Nt, fn)
      write(fn,*) ""
      deallocate(x,u)
    end do 
    close(fn)
    fn = fn + 1
  
  print *, " "
  print *, " ----------------------------------------------------------- "
  print *, " "
  print *, " Problem 7:"
  print *, " "

    open(fn, file='prob7.dat')
    cfl_array = (/0.8, 0.8, 1.0, 1.0, 1.2, 1.2/)
    nx_array = (/32, 128, 32, 128, 32, 128/)
    xstart = -1
    xstop = 1 
    tstart = 0
    tstop = 1 

    ! initializing data
    do j = 1, 6
      cfl = cfl_array(j)
      a = 1
      Nx = nx_array(j)
      dx = (xstop-xstart)/(Nx-1)
      dt = dx*cfl/a
      Nt = ceiling((tstop-tstart)/dt)
      allocate(x(Nx,Nt), u(Nx,Nt))
      do i = 1, Nx
        x(i,1) = xstart + (i-1)*dx
        if (abs(x(i,1)) .lt. 1./3) then
          u(i,1) = 1 ! might have to average these together over the cells
        else 
          u(i,1) = 0
        end if
      end do  

      ! running timestepper
      do i = 1, Nt-1
        u(:,i+1) = LF_update_1D(u(:,i),u(:,i),dt,dx,Nx)
      end do 
      write(fn,*) Nx, Nt
      call writemat(u, Nx, Nt, fn)
      write(fn,*) ""
      deallocate(x,u)
    end do 
    close(fn)
    fn = fn + 1

  print *, " "
  print *, " ----------------------------------------------------------- "
  print *, " "
  print *, " Problem 8:"
  print *, " "

    open(fn, file='prob8.dat')
    cfl_array = (/0.8, 0.8, 1.0, 1.0, 1.2, 1.2/)
    nx_array = (/32, 128, 32, 128, 32, 128/)
    xstart = 0
    xstop = 1 
    tstart = 0
    tstop = 1 

    ! initializing data
    do j = 1, 6
      cfl = cfl_array(j)
      a = 1
      Nx = nx_array(j)
      dx = (xstop-xstart)/(Nx-1)
      dt = dx*cfl/a
      Nt = ceiling((tstop-tstart)/dt)
      allocate(x(Nx,Nt), u(Nx,Nt))
      do i = 1, Nx
        x(i,1) = xstart + (i-1)*dx
      end do  
      u = sin(2*pi*x) ! might have to average these together over the cells
      
      ! running timestepper
      do i = 1, Nt-1
        u(:,i+1) = LW_update_1D(u(:,i),a,dt,dx,Nx)
      end do 
      write(fn,*) Nx, Nt
      call writemat(u, Nx, Nt, fn)
      write(fn,*) ""
      deallocate(x,u)
    end do 
    close(fn)
    fn = fn + 1

  print *, " "
  print *, " ----------------------------------------------------------- "
  print *, " "
  print *, " Problem 9:"
  print *, " "

    open(fn,file='prob9.dat')
    cfl_array = (/0.8, 0.8, 1.0, 1.0, 1.2, 1.2/)
    nx_array = (/32, 128, 32, 128, 32, 128/)
    xstart = -1
    xstop = 1 
    tstart = 0
    tstop = 1 

    ! initializing data
    do j = 1, 6
      cfl = cfl_array(j)
      a = 1
      Nx = nx_array(j)
      dx = (xstop-xstart)/(Nx-1)
      dt = dx*cfl/a
      Nt = ceiling((tstop-tstart)/dt)
      allocate(x(Nx,Nt), u(Nx,Nt))
      do i = 1, Nx
        x(i,1) = xstart + (i-1)*dx
        if (abs(x(i,1)) .lt. 1./3) then
          u(i,1) = 1 ! might have to average these together over the cells
        else 
          u(i,1) = 0
        end if
      end do  
      
      ! running timestepper
      do i = 1, Nt-1
        u(:,i+1) = LW_update_1D(u(:,i),a,dt,dx,Nx)
      end do 
      write(fn,*) Nx, Nt
      call writemat(u, Nx, Nt, fn)
      write(fn,*) ""
      deallocate(x,u)
    end do 
    close(fn)

  print *, " "
  print *, " ----------------------------------------------------------- "
end program driver
