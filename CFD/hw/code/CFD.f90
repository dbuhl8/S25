!File: Numerical Differential Equations Module
!Author: Dante Buhl
!Dependencies: IBVP_1DCN uses a LAPACK routine for a linear solve
!Purpose: To store numerical methods for Differential Equations

module CFD

  use LinAl ! contains pi as a parameter
  use NumDE

  implicit none

  ! Global vars
  !integer, parameter :: kr=kind(dble(1.0))
  real(kind=kr) :: cfl_number

  contains 

  function LF_update_1D(U,a,dt,dx,nx) result(NU)
    implicit none
    ! here U and FU should be nx by 1 dimensional arrays, the output of the
    ! method should be another nx by 1 array 
    real(kind=kr) :: U(:), a, dt, dx
    integer :: nx

    real(kind=kr), dimension(nx) :: NU

    ! inside the boundary
    !NU(2:nx-1) = (U(1:nx-2) + U(3:nx))/2 - (a*dt/dx)*(U(3:nx) - U(1:nx-2))/2
    NU(2:nx-1) = (1 - a*dt/dx)*U(3:nx)/2 + (1 + a*dt/dx)*U(1:nx-2)/2

    ! on the boundary (periodic conditions)
    NU(1) = (U(nx) + U(2))/2 - (a*dt/dx)*(-U(nx) + U(2))/2
    NU(nx) = (U(nx-1) + U(1))/2 - (a*dt/dx)*(U(1) - U(nx-1))/2
  end function LF_update_1D

  function LW_update_1D(U, a, dt, dx, nx) result(NU)
    implicit none
    ! input
    real(kind=kr) :: U(:), a, dt, dx
    integer :: nx
    ! output 
    real(kind=kr), dimension(nx) :: NU

    ! in the boundary
    NU(2:nx-1) = U(2:nx-1) - (a*dt/dx)*(-U(1:nx-2) + U(3:nx))/2 +&
      (a*dt/dx)**2*(U(1:nx-2) -2*U(2:nx-1)+ U(3:nx))/2

    ! on the boundary (periodic conditions enforced)
    NU(1) = U(1) - (a*dt/dx)*(U(2) - U(nx))/2 + (a*dt/dx)**2*(U(2) -&
      2*U(1) + U(nx))/2
    NU(nx) = U(nx) -(a*dt/dx)*(-U(nx-1) + U(1))/2 +(a*dt/dx)**2*(U(1) -&
      2*U(nx) + U(nx-1))/2
  end function LW_update_1D

  function Upwind_update_1D(U,a,dt,dx,nx) result(NU)
    implicit none
    ! input
    real(kind=kr) :: U(:), a, dt, dx
    integer :: nx
    ! output 
    real(kind=kr), dimension(nx) :: NU

    NU(2:nx) = U(2:nx) - (a*dt/dx)*(U(2:nx) - U(1:nx-1))
    NU(1) = U(1) - (a*dt/dx)*(U(1) - U(nx))
  end function Upwind_update_1D

  include 'fog.f90'
  include 'plm.f90'
  include 'bc.f90'
  include 'cfl.f90'
  include 'grid_init.f90'
  include 'advect_init.f90' 
  include 'write_data.f90'

end module CFD
 
