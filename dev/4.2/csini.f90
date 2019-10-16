subroutine csini(CL,crlen,cs2d,cs4d,lmax)

! Computes the correlation length, and second and fourth derivatives
! for an autocorrelation function expressed in Legendre series with
! unnormalized coefficients. The coefficients are properly normalized
! in output.


use common
implicit none
integer :: l,lmax
real (kind=dp) :: crlen,cs2d,cs4d,norm
real (kind=dp), dimension(0:256) :: CL

!if (lmin.lt.2) stop 'Trouble in CSINI: lmin < 2.'

! Normalization:

norm=0.0_dp
do l=0,lmax
    norm=norm+CL(l)
end do
do l=0,lmax
    CL(l)=CL(l)/norm
end do

! Derivatives and correlation length:

cs2d=0.0_dp
cs4d=0.0_dp
do l=1,lmax
    cs2d=cs2d-CL(l)*l*(l+1.0_dp)/2.0_dp
    cs4d=cs4d+CL(l)*l*(l+1.0_dp)*(3.0_dp*l**2+3.0_dp*l-2.0_dp)/8.0_dp
end do
if (cs2d.eq.0.0_dp) then
    crlen=2.0_dp
else
    crlen=1.0_dp/sqrt(-cs2d)
endif
end subroutine csini


