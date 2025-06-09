

subroutine grid_init(xstart, xstop, nx, x_cell, x_bound, dx)
  use LinAL
  implicit none
  real(kind=kr), allocatable :: x_cell(:), x_bound(:)
  real(kind=kr) :: xstart, xstop, dx
  integer :: nx, i

  allocate(x_cell(nx), x_bound(nx+1))

  dx = (xstop-xstart)/(nx)
  x_bound(1) = xstart
  do i = 1, nx
    x_bound(i+1) = xstart + i*dx
    x_cell(i) = xstart + (i-0.5)*dx
  end do 
end subroutine grid_init
