
function vl(a,b) result(s)
  use LinAl
  implicit none
  real(kind=kr) :: a, b, s
  if (a*b .le. 0) then
    s = 0
  else 
    s = 2*a*b/(a+b)
  end if
end function vl
