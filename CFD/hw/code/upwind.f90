

function upwind(U, dx,nx,ngc) result(sl)  
  use NumDE 
  implicit none
  real(kind=kr) :: U(:), dx
  integer :: nx, ngc
  real(kind=kr), dimension(nx+2*ngc-2) :: sl
  sl = (U(ngc:nx+ngc+1)-U(ngc-1:nx+ngc))/dx
end function upwind
