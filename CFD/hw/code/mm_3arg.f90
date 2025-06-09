
function mm_3arg(a,b,c) result(s)
  use LinAl
  implicit none
  real(kind=kr) :: a, b, c, s
  if (a*b*c .gt. 0) then
    if((abs(a) .le. abs(b)) .and. (abs(a) .le. abs(c))) then
      s = a
    else if (abs(b) .le. abs(c)) then
      s = b
    else 
      s = c
    end if
  else 
    s = 0
  end if
end function mm_3arg
