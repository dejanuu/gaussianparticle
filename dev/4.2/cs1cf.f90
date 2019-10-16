subroutine cs1cf(CSCF,nu,lmin,lmax)

! Returns the Legendre coefficients for the correlation
! function with power-law Legendre coefficients.


use common
implicit none

integer :: l,lmin,lmax
real (kind=dp) :: nu,norm
real (kind=dp), dimension(0:256) :: CSCF

!WRITE(*,*) nu,lmin,lmax
do l=0,lmin-1
    CSCF(l)=0.0_dp
end do

norm=0.0_dp
do l=lmin,lmax
    CSCF(l)=1.0_dp/l**nu
    norm=norm+CSCF(l)
end do

do l=lmin,lmax
    CSCF(l)=CSCF(l)/norm
end do

end subroutine cs1cf






