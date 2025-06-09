

function minmod(U, dx, nx) result(sl)
  use NumDE
  implicit none

  real(kind=kr) :: U(:), dx
  integer :: nx, i
  real(kind=kr), dimension(nx) :: sl

  do i = 1, nx
    sl(i) = mm_2arg((U(i+1)-U(i))/dx, (U(i+2)-U(i+1))/dx)
  end do 
end function minmod
