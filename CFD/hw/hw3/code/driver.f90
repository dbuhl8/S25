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
  integer, parameter :: kr=kind(dble(1.0))
  integer :: i, j

  ! Simvars
  real(kind=kr), parameter :: invRo0 = 2., invRoF=8., Fr=0.1
  real(kind=kr) :: u0=2*invRo0, uF=2*invRoF, u ! max rms velocity
  real(kind=kr) :: dt0=.0005, dt ! timestep
  real(kind=kr) :: m=0.02 ! m = du/dt (averaged growth rate of u)
  real(kind=kr) :: cfl_sf = 0.8
  real(kind=kr) :: Re=600, Pe=60
  real(kind=kr) :: cfl_a=5.61,cfl_b=0.69

  ! Loop Vars
  integer :: counter = 0
 
  u = u0 
  dt = dt0
  
  print *, "Starting Loop"
  do while (u .lt. uF)
    u = u + m*dt
    dt = (dt0*u0)/u
    counter = counter + 1
  end do 
  print *, "Exitted Loop"
  print *, "Timesteps Needed:       ", counter
  print *, "Final velocity:         ", u
  print *, "Final Timestep Length:  ", dt

end program driver
