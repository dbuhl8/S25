
function mm_2arg(a,b) result(s)
  use LinAl
  implicit none
  real(kind=kr) :: a, b, s
  if (a*b .gt. 0) then
    if(abs(a) .le. abs(b)) then
      s = a
    else 
      s= b
    end if
  else 
    s = 0
  end if
end function mm_2arg
