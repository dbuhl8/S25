

function advect_update(U, nx, dx, dt, method, submethod) result(NU)
  use LinAL
  use NumDE
  use CFD
  implicit none
  real(kind=kr) :: U(:), dx, dt
  integer :: nx
  character :: method, submethod
  real(kind=kr), dimension(nx) :: NU

  if(method .eq. 'f') then
    NU = FOG(U,dt,dx,nx)
  else if (method .eq. 'p') then
    NU = PLM(U,dt,dx,nx,submethod)
  else 
    print *, 'Evolution Method not recognized'
  end if
end function advect_update
