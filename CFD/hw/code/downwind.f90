

function downwind(U, dx,nx,ngc) result(sl)  
  use NumDE 
  implicit none
  real(kind=kr) :: U(:), dx
  integer :: nx,ngc
  real(kind=kr), dimension(nx+2*ngc-2) :: sl
  sl = (U(ngc+1:nx+ngc+2)-U(ngc:nx+ngc+1))/dx
end function downwind
