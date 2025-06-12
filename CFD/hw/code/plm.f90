
function PLM(U, dt, dx, nx,ngc, sl_type) result(NU)
  ! for the burgers equation 
  use NumDE
  implicit none
  ! return U(n+1) = U(n) + (del T)/(del X)(F_god_+ - F_god_-)
  ! need to compute the Godunov Fluxes (which at the cell boundaries are
  ! constant), and then return U(n+1)
  ! need a where construct in order to evalutate God_fluxes
  real(kind=kr) :: U(:), dt, dx, s, c
  integer :: nx, i,ngc, k, j
  character :: sl_type
  real(kind=kr), dimension(nx) :: NU
  real(kind=kr), dimension(nx+2*ngc) :: sl, F

  ! this is programmed to have ghost cells at each end of the array
  ! these cells will be routinely updated after each update

  ! computing Fp and Fm
  sl = sl_limiter(sl_type, U, dx, nx,ngc)
  do i = 1, nx+ngc+1
    !k = i - ngc
    !j = i - 1
    ! Fp
    s = (U(i) + U(i+1))/2
    c = s*dt/dx
    if (U(i) .ge. U(i+1)) then
      if (s .gt. 0) then
        F(i) = ((U(i) + (dx/2)*sl(i)*(1 - c))**2)/2
      else 
        F(i) = ((U(i+1) - (dx/2)*sl(i+1)*(1 + c))**2)/2
      end if
    else 
      if (U(i) .ge. 0) then 
        F(i) = ((U(i) + (dx/2)*sl(i)*(1 - c))**2)/2
      else if (U(i+1) .le. 0) then
        F(i) = ((U(i+1) - (dx/2)*sl(i+1)*(1 + c))**2)/2
      else 
        F(i) = 0
      end if
    end if 
    ! Fm
    !if (U(i-1) .ge. U(i)) then
      !s = (U(i-1) + U(i))/2
      !if (s .gt. 0) then
        !Fm(k) = ((U(i-1) + (dx/2)*sl(j-1)*(1 - s*dt/dx))**2)/2
      !else 
        !Fm(k) = ((U(i) - (dx/2)*sl(j)*(1 + s*dt/dx))**2)/2
      !end if
    !else 
      !if (U(i-1) .ge. 0) then 
        !Fm(k) = ((U(i-1) + (dx/2)*sl(j-1)*(1 - s*dt/dx))**2)/2
      !else if ((U(i-1) .lt. 0) .and. (0 .lt. U(i))) then
        !Fm(k) = 0
      !else 
        !Fm(k) = ((U(i) - (dx/2)*sl(j)*(1 + s*dt/dx))**2)/2
      !end if
    !end if 
  end do 

  NU = U(ngc+1:nx+ngc) - (dt/dx)*(F(ngc+1:nx+ngc)-F(ngc:ngc+nx-1))

  contains 
    include 'upwind.f90'
    include 'downwind.f90'
    include 'centered.f90'
    include 'minmod.f90'
    include 'mc.f90'
    include 'vanLeers.f90'

    function sl_limiter(sl_type, U, dx, nx,ngc) result(sl)
      implicit none
      character :: sl_type  
      real(kind=kr) :: U(:), dx
      integer :: nx,ngc
      real(kind=kr), dimension(nx+2*ngc) :: sl

      if (sl_type .eq. 'u') then
      ! upwind slope limiter
        sl = upwind(U, dx, nx,ngc)
      else if (sl_type .eq. 'd') then 
      ! downwind slope limiter
        sl = downwind(U, dx, nx,ngc)
      else if (sl_type .eq. 'c') then 
      ! centered slope limiter
        sl = centered(U, dx, nx,ngc)
      else if (sl_type .eq. 'm') then 
      ! minmod slope limiter
        sl = minmod(U, dx, nx,ngc)
      else if (sl_type .eq. 'o') then 
      ! mc slope limiter
        sl = mc(U, dx, nx,ngc)
      else if (sl_type .eq. 'v') then 
      ! van leer slope limiter
        sl = vanLeers(U, dx, nx,ngc)
      end if
    end function sl_limiter 
end function PLM
