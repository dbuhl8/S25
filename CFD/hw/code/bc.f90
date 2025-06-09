

subroutine bc(U, nx, type)
  use NumDE
  implicit none
  real(kind=kr) :: U(:)
  integer :: nx
  character :: type

  if (type .eq. 'p') then
    ! periodic BC (fills ghost cells with opposite end)
    U(1) = U(nx+1)
    U(nx+2) = U(2)
  else if (type .eq. 'o') then
    ! outflow bc (not specified anywhere in notes)
  else 
    print *, ' Boundary Condition Type not recognized '
  end if

end subroutine bc
