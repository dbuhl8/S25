

function centered(U, dx,nx) result(sl)  
  use NumDE 
  implicit none
  real(kind=kr) :: U(:), dx
  integer :: nx
  real(kind=kr), dimension(nx) :: sl
  sl = (U(3:nx+2)-U(1:nx))/(2*dx)
end function centered
