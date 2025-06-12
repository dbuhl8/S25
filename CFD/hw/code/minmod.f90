

function minmod(U, dx, nx,ngc) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i, ngc
  real(kind=kr), dimension(nx+2*ngc) :: sl

  sl = 0
  do i = 2, nx+2*ngc-1
    sl(i) = mm_2arg((U(i)-U(i-1))/dx, (U(i+1)-U(i))/dx)
  end do 
end function minmod
