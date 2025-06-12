

subroutine advect_init(U_cell, U_grid, x_grid, nx,ngc, dx,ic)
  use NumDE
  implicit none

  real(kind=kr) :: U_grid(:), dx, U_cell(:), x_grid(:)
  integer :: nx, i, ic, ngc

  if (ic .eq. 1) then ! 7.38
    where(x_grid .le. 0.5)
      U_grid = 2
    elsewhere(x_grid .gt. 0.5)
      U_grid = -1
    end where
  else if (ic .eq. 2) then ! 7.39
    where(x_grid .le. 0.5)
      U_grid = 1
    elsewhere(x_grid .gt. 0.5) 
      U_grid = -1
    end where
  else if (ic .eq. 3) then  ! 7.40
    where(x_grid .le. 0.5) 
      U_grid = 1
    elsewhere(x_grid .gt. 0.5)
      U_grid = -2
    end where
  else if (ic .eq. 4) then  !7.41
    where(x_grid .le. 0.5)
      U_grid = -1
    elsewhere(x_grid .gt. 0.5)
      U_grid = 1
    end where
  else if (ic .eq. 5) then  !7.42
    U_grid = sin(2*pi*x_grid)
  else if (ic .eq. 6) then  !7.43
    where(x_grid .le. 0.3)
      U_grid = 2
    elsewhere((x_grid .le. 0.6) .and. (x_grid .gt. 0.3))
      U_grid = -1
    elsewhere(x_grid .gt. 0.6)
      U_grid = 3
    end where
  else if (ic .eq. 7) then  !7.44
    where(x_grid .le. 0.3)
      U_grid = -1
    elsewhere((x_grid .le. 0.6) .and. (x_grid .gt. 0.3))
      U_grid = 2
    elsewhere(x_grid .gt. 0.6)
      U_grid = -2
    end where
  else if (ic .eq. 8) then  !7.45
    where(x_grid .le. 0.3)
      U_grid = 4
    elsewhere((x_grid .le. 0.6) .and. (x_grid .gt. 0.3))
      U_grid = 2
    elsewhere(x_grid .gt. 0.6)
      U_grid = -1
    end where
  else if (ic .eq. 9) then  !7.46
    where(x_grid .le. 0.3)
      U_grid = 4
    elsewhere((x_grid .le. 0.6) .and. (x_grid .gt. 0.3))
      U_grid = 0
    elsewhere(x_grid .gt. 0.6)
      U_grid = -6
    end where
  else if (ic .eq. 10) then  !7.47
    where(x_grid .le. 0.3)
      U_grid = 4
    elsewhere((x_grid .le. 0.6) .and. (x_grid .gt. 0.3))
      U_grid = 0
    elsewhere(x_grid .gt. 0.6)
      U_grid = -4
    end where
  else if (ic .eq. 11) then  !7.48
    ! sin wave IC
    U_grid = sin(2*pi*x_grid)
  else 
    print *, 'Initial Condition not recognized'
  end if

  ! cell averaging
  do i = 1, nx
    U_cell(i+ngc) = cell_avg(U_grid(i+ngc:i+ngc+1), dx)
  end do 

  contains 
    function cell_avg(U,dx) result(avg)
      implicit none
      real(kind=kr) :: U(:), dx, avg
      avg = sum(U)/2
    end function cell_avg
end subroutine advect_init
