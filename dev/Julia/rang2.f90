
subroutine rang2(r1)

! Returns a normally distributed random deviate with zero mean and
! unit variance.

implicit none
integer, parameter ::                              &
sp = kind(1.0),                                  &
dp = selected_real_kind(2*precision(1.0_sp)),    &
qp = selected_real_kind(2*precision(1.0_dp))

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
10  call random_number(ran2)
r1=2.0_dp*ran2-1.0_dp
call random_number(ran2)
r2=2.0_dp*ran2-1.0_dp
q1=r1**2+r2**2
if (q1.ge.1.0_dp .or. q1.le.0.0_dp) goto 10

q2=sqrt(-2.0_dp*log(q1)/q1)
r1=r1*q2
r2=r2*q2
flag=1

end subroutine rang2




