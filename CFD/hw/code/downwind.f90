

function downwind(U, dx,nx,ngc) result(sl)  
  use NumDE 
  implicit none
  real(kind=kr) :: U(:), dx
  integer :: nx,ngc
  real(kind=kr), dimension(nx+2*ngc) :: sl
  sl(2:nx+2*ngc-1) = (U(2:nx+2*ngc-1)-U(1:nx+2*ngc-2))/dx
end function downwind
