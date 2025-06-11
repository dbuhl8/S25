

function vanLeers(U, dx, nx,ngc) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i, ngc
  real(kind=kr), dimension(nx+2*ngc-2) :: sl

  do i = ngc-1, nx+ngc
    sl(i) = vl((U(i+1)-U(i))/dx, (U(i+2)-U(i+1))/dx)
  end do 
end function vanLeers
