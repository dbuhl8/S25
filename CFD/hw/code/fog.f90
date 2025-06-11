


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
  real(kind=kr), dimension(nx) :: NU, Fp, Fm

  ! this is programmed to have ghost cells at each end of the array
  ! these cells will be routinely updated after each update

  ! computing Fp and Fm
  do i = ngc+1, nx+ngc
    k = i - ngc
    ! Fp
    if (U(i) .ge. U(i+1)) then
      s = (U(i) + U(i+1))/2
      if (s .gt. 0) then
        Fp(k) = (U(i)**2)/2
      else 
        Fp(k) = (U(i+1)**2)/2
      end if
    else 
      if (U(i) .ge. 0) then 
        Fp(k) = (U(i)**2)/2
      else if ((U(i) .lt. 0) .and. (0 .lt. U(i+1))) then
        Fp(k) = 0
      else 
        Fp(k) = (U(i+1)**2)/2
      end if
    end if 
    ! Fm
    if (U(i-1) .ge. U(i)) then
      s = (U(i-1) + U(i))/2
      if (s .gt. 0) then
        Fm(k) = (U(i-1)**2)/2
      else 
        Fm(k) = (U(i)**2)/2
      end if
    else 
      if (U(i-1) .ge. 0) then 
        Fm(k) = (U(i-1)**2)/2
      else if ((U(i-1) .lt. 0) .and. (0 .lt. U(i))) then
        Fm(k) = 0
      else 
        Fm(k) = (U(i)**2)/2
      end if
    end if 
  end do 

  NU = U(ngc+1:nx+ngc) - (dt/dx)*(Fp-Fm)
end function FOG
