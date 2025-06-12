
function mc(U, dx, nx,ngc) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i, ngc
  real(kind=kr), dimension(nx+2*ngc) :: sl

  do i = 2, nx+2*ngc-1
    sl(i) = mm_3arg((U(i+1)-U(i-1))/(2*dx),2*(U(i+1)-U(i))/dx, 2*(U(i)-U(i-1))/dx)
  end do 
end function mc
