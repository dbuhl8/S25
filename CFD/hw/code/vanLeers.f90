

function vanLeers(U, dx, nx,ngc) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i, ngc
  real(kind=kr), dimension(nx+2*ngc) :: sl

  do i = 2, nx+2*ngc-1
    sl(i) = vl((U(i)-U(i-1))/dx, (U(i+1)-U(i-1))/dx)
  end do 
end function vanLeers
