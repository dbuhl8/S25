!File: Numerical Differential Equations Module
!Author: Dante Buhl
!Dependencies: IBVP_1DCN uses a LAPACK routine for a linear solve
!Purpose: To store numerical methods for Differential Equations

module CFD

  use LinAl ! contains pi as a parameter
  use NumDE

  implicit none

  ! Global vars
  integer, parameter :: kr=kind(dble(1.0))
  real(kind=kr) :: cfl_number

  contains 

  function LF_update_1D(U,FU,dt,dx,nx) result(NU)
    implicit none
    ! here U and FU should be nx by 1 dimensional arrays, the output of the
    ! method should be another nx by 1 array 
    real(kind=kr) :: U(:,:), FU(:,:), dt, dx
    integer :: nx

    real(kind=kr), dimension(nx, 1) :: NU

    ! inside the boundary
    NU(2:nx-1,1) = (U(1:nx-2,1) + U(3:nx,1))/2 - (dt/dx)*(FU(1:nx-2,1) - FU(3:nx,1))/2

    ! on the boundary (periodic conditions)
    NU(1,1) = (U(nx,1) + U(2,1))/2 - (dt/dx)*(FU(nx,1) - FU(2,1))/2
    NU(nx,1) = (U(1,1) + U(nx-1,1))/2 - (dt/dx)*(FU(1,1) - FU(nx-1,1))/2
  end function LF_update()

  subroutine LW_method()

  end subroutine LW_method

end module CFD
 
