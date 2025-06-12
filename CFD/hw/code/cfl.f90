

subroutine cfl(U, dx, sf, dt)
  use NumDE
  implicit none
  real(kind=kr) :: dx, sf, U(:), dt, eps=1e-4
  dt = sf*dx/max(maxval(abs(U)),eps)
end subroutine cfl
