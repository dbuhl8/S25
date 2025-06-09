

subroutine write_data(A, ma, na, fn)
  use NumDE
  use LinAl
  !writes matrix A out to a file using a format
  implicit none

  real :: A(:,:)
  integer :: ma, na, fn, i

  do i = 1, na
    write(fn, "("//trim(str(ma))//"(F32.16, ' '))") A(:,i)
  end do 
end subroutine write_data
