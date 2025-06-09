

function upwind(U, dx,nx) result(sl)  
  use NumDE 
  implicit none
  real(kind=kr) :: U(:), dx
  integer :: nx
  real(kind=kr), dimension(nx) :: sl
  sl = (U(2:nx+1)-U(1:nx))/dx
end function upwind
