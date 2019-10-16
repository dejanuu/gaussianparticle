
subroutine vproteut(X,CA,SA)

! Transpose of vector rotation using Euler angles. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

use common
implicit none

real(kind=dp), dimension(3) ::  X(3),CA(3),SA(3)

call vrotz(X,CA(3),-SA(3))
call vroty(X,CA(2),-SA(2))
call vrotz(X,CA(1),-SA(1))

end subroutine vproteut




subroutine vroty(X,c,s)

! Vector rotation about the y-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

use common
implicit none

real(kind=dp) :: c,s,q
real(kind=dp), dimension(3) :: X

q   = c*X(3)+s*X(1)
X(1)=-s*X(3)+c*X(1)
X(3)=q

end subroutine vroty



subroutine vrotz(X,c,s)

! Vector rotation about the z-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

use common
implicit none

real(kind=dp) :: c,s,q
real(kind=dp), dimension(3) :: X

q   = c*X(1)+s*X(2)
X(2)=-s*X(1)+c*X(2)
X(1)=q

end subroutine vrotz








