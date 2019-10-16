
subroutine sgscf(ACF,BCF,SCFSTD,lmax)

! Generates the sample spherical harmonics coefficients for the
! logarithmic radial distance of the G-sphere. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

use common
implicit none

integer  :: l,lmax,m
real (kind=dp) :: rn
real(kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD


do l=0,lmax
    call rang2(rn)
    ACF(l,0)=rn*SCFSTD(l,0)
    BCF(l,0)=0.0_dp
end do

do m=1,lmax
    do l=m,lmax
        call rang2(rn)
        ACF(l,m)=rn*SCFSTD(l,m)
        call rang2(rn)
        BCF(l,m)=rn*SCFSTD(l,m)
    end do
end do


end subroutine sgscf


