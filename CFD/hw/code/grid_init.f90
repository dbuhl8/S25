

subroutine grid_init(xstart, xstop, nx,ngc, x_cell, x_bound, dx)
  use LinAL
  implicit none
  real(kind=kr), allocatable :: x_cell(:), x_bound(:)
  real(kind=kr) :: xstart, xstop, dx
  integer :: nx, i,ngc

  allocate(x_cell(nx+2*ngc), x_bound(nx+1+2*ngc))

  dx = (xstop-xstart)/(nx)
  x_bound(1) = xstart-(ngc*dx)
  do i = 1, nx+2*ngc
    x_bound(i+1) = xstart + (i-ngc)*dx
    x_cell(i) = xstart + (i-ngc-0.5)*dx
  end do 
end subroutine grid_init
