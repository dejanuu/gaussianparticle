module lm

use specfunc
use voper
use parameters

contains

subroutine ARDAR(AR,DAR,ACF,BCF,AT,MUT,PHIT,ntri,nacf,lmin,lmax)

!use constants	 
implicit none
integer::ntri,nacf,l,m,lmin,lmax,j1,j2
real(kind=dp),dimension(260000)::AR,AT,MUT,PHIT
real(kind=dp)::DAR(7200,256),CPHI(0:256),SPHI(0:256),SGS,mu,phi,&
		 armn0,armn
real(kind=dp),dimension(0:256,0:256)::ACF,BCF,LEGP

armn0=4.0d0*pi/ntri
CPHI(0)=1.0d0
SPHI(0)=0.0d0

do j1=1,7200
    do j2=10,256
         DAR(j1,j2)=0.0d0
	end do
end do


armn=0.0d0
do j1=1,ntri
    	mu=MUT(j1)
	phi=PHIT(j1)

! Precomputation of sines, cosines, and associated Legendre functions:

    call LEGAN(LEGP,mu,lmax,0)
    do m=1,lmax
	call LEGAN(LEGP,mu,lmax,m)
	CPHI(m)=cos(m*phi)
	SPHI(m)=sin(m*phi)
    end do
    LEGP(0,0)=1.0d0/sqrt(2.0d0)

! Sum up:

    SGS=0.0d0
    do l=lmin,lmax
	SGS=SGS+LEGP(l,0)*ACF(l,0)
    end do
    do m=1,lmax
	do l=max(m,lmin),lmax
		SGS=SGS+LEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
	end do
    end do

AR(j1)=exp(SGS)*AT(j1)
armn=armn+AR(j1)
!	do j2=10,nacf
!		DAR(j1,j2)=LEGP(NCF(j2,1),NCF(j2,2))*CPHI(NCF(j2,2))*AR(j1)
!	end do
!	do j2=nacf+1,npar
!			DAR(j1,j2)=LEGP(NCF(j2,1),NCF(j2,2))*SPHI(NCF(j2,2))*AR(j1)
!	end do
end do
armn=armn/ntri

end


end module lm
