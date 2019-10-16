
subroutine rang2(r1)

! Returns a normally distributed random deviate with zero mean and
! unit variance.

use common
implicit none

integer  :: flag
real (kind=dp) :: ran2,q1,q2,r1,r2
save flag,r2
data flag/0/


if (flag.eq.1) then
    r1=r2
    flag=0
    return
endif

flag=0
do
  call random_number(ran2)
  r1=2.0_dp*ran2-1.0_dp
  call random_number(ran2)
  r2=2.0_dp*ran2-1.0_dp
  q1=r1**2+r2**2
  if(q1 > 0.0_dp .and. q1 < 1.0_dp) exit
end do

q2=sqrt(-2.0_dp*log(q1)/q1)
r1=r1*q2
r2=r2*q2
flag=1

end subroutine rang2




