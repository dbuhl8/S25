


function FOG(U, dt, dx, nx,ngc) result(NU)
  ! for the burgers equation 
  use NumDE
  implicit none
  ! return U(n+1) = U(n) + (del T)/(del X)(F_god_+ - F_god_-)
  ! need to compute the Godunov Fluxes (which at the cell boundaries are
  ! constant), and then return U(n+1)
  ! need a where construct in order to evalutate God_fluxes
  real(kind=kr) :: U(:), dt, dx, s
  integer :: nx, i,ngc, k
  real(kind=kr), dimension(nx+2*ngc) :: F
  real(kind=kr), dimension(nx) :: NU

  ! this is programmed to have ghost cells at each end of the array
  ! these cells will be routinely updated after each update

  ! computing Fp and Fm
  do i = 1, nx+2*ngc-1
    !k = i - ngc
    ! Fp
    if (U(i) .ge. U(i+1)) then
      s = (U(i) + U(i+1))/2
      if (s .gt. 0) then
        F(i) = (U(i)**2)/2
      else 
        F(i) = (U(i+1)**2)/2
      end if
    else 
      if (U(i) .ge. 0) then
        F(i) = (U(i)**2)/2
      else if (U(i+1) .le. 0) then
        F(i) = (U(i+1)**2)/2
      else
        F(i) = 0
      end if
    end if 
  end do 

  NU = U(ngc+1:nx+ngc) - (dt/dx)*(F(ngc+1:nx+ngc)-F(ngc:nx+ngc-1))
end function FOG
