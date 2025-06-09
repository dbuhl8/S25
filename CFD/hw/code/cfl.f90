

subroutine cfl(U, dx, sf, dt)
  use NumDE
  implicit none
  real(kind=kr) :: dx, sf, U(:), dt
  dt = sf*dx/maxval(abs(U))
end subroutine cfl
