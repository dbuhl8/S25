

function vanLeers(U, dx, nx) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i
  real(kind=kr), dimension(nx) :: sl

  do i = 1, nx
    sl(i) = vl((U(i+1)-U(i))/dx, (U(i+2)-U(i+1))/dx)
  end do 

end function vanLeers
