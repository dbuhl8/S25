

subroutine bc(U, nx,ngc, type)
  use NumDE
  implicit none
  real(kind=kr) :: U(:)
  integer :: nx, ngc
  character :: type

  if (type .eq. 'p') then
    ! periodic BC (fills ghost cells with opposite end)
    U(1:ngc) = U(nx+1:nx+ngc)
    U(nx+ngc+1:nx+2*ngc) = U(ngc+1:2*ngc)
  else if (type .eq. 'o') then
    ! outflow bc (specified in HW3)
    U(1:ngc) = U(ngc+1)
    U(nx+ngc+1:nx+2*ngc) = U(nx+ngc)
  else 
    print *, ' Boundary Condition Type not recognized '
  end if

end subroutine bc
