module randev

! Random deviates:
!
! RNDU: uniform distribution within [0, 1)
!
! Required modules:
!

use parameters      ! Parameters
        
implicit none
  
contains


subroutine RNDG(r1)

! Returns a normally distributed random deviate with zero mean and 
! unit variance. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(inout) :: r1
integer :: flg,irnd
real(kind=dp) :: q1,q2,r2,x
save flg,r2
data flg/0/
common irnd

if (flg.eq.1) then
 r1=r2
 flg=0
 return
endif

flg=0

q1=1000.0d0
do while(q1>=1.0d0 .or. q1<=0.0d0)
    CALL RANDOM_NUMBER(x)
    r1=2.0d0*x-1.0d0
    CALL RANDOM_NUMBER(x)
    r2=2.0d0*x-1.0d0
    q1=r1**2+r2**2
enddo

q2=sqrt(-2.0d0*log(q1)/q1)
r1=r1*q2
r2=r2*q2
flg=1
end subroutine RNDG

end module randev

