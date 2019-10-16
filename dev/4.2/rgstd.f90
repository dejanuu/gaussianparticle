! rgstd
! rgs
! sgs


subroutine rgstd(X,N,MU,PHI,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)

! Discrete triangle representation for a sample G-sphere.

use common
implicit none

integer  :: nnod,ntri,lmax,j1,j2
integer , dimension(260000,3) :: IT
real (kind=dp) :: RSPH1,rmax,beta,r,nu
real (kind=dp), dimension(130000,3) :: X
real (kind=dp), dimension(260000,3) :: N
real (kind=dp), dimension(130000) :: MU, PHI
real (kind=dp), dimension(0:256,0:256) :: ACF, BCF
real (kind=dp), dimension(3) :: X1, X2, X3




! Node coordinates:

rmax=0.0_dp
do j1=1,nnod
    r=RSPH1(ACF,BCF,MU(j1),PHI(j1),beta,lmax)
    nu=sqrt(1.0_dp-MU(j1)**2)
    X(j1,1)=r*nu*cos(PHI(j1))
    X(j1,2)=r*nu*sin(PHI(j1))
    X(j1,3)=r*MU(j1)
    if (r.ge.rmax) rmax=r
end do

! Outer unit triangle normals:

do j1=1,ntri
    do j2=1,3
        r=X(IT(j1,1),j2)
        X1(j2)=X(IT(j1,2),j2)-r
        X2(j2)=X(IT(j1,3),j2)-r
    end do
    call provec(X3,X1,X2)
    r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
    do j2=1,3
        N(j1,j2)=X3(j2)/r
    end do
end do
end subroutine rgstd


REAL(KIND=dp) FUNCTION RGS(ACF,BCF,mu,phi,beta,lmin,lmax)

! Radial distance in a given direction for a sample G-sphere.


use common
implicit none

integer  :: lmin,lmax
real (kind=dp) :: SGS,mu,phi,beta
real (kind=dp), dimension(0:256,0:256) :: ACF, BCF

RGS=exp(SGS(ACF,BCF,mu,phi,lmin,lmax)-0.5_dp*beta**2)
end

REAL(KIND=dp) FUNCTION SGS(ACF,BCF,mu,phi,lmin,lmax)

! Logarithmic radial distance in a given direction for a sample G-sphere.


use common
implicit none

integer  :: l,m,lmin,lmax
real (kind=dp) :: mu,phi
real (kind=dp), dimension(0:256,0:256) :: ACF, BCF, LEGP
real (kind=dp), dimension(256) :: CPHI, SPHI




if (lmax.eq.0) then
    SGS=ACF(0,0)
    return
endif

! Precomputation of sines, cosines, and associated Legendre functions:

call LEGA(LEGP,mu,lmax,0)
do m=1,lmax
    call LEGA(LEGP,mu,lmax,m)
    CPHI(m)=cos(m*phi)
    SPHI(m)=sin(m*phi)
end do
LEGP(0,0)=1.0_dp

! Sum up:

SGS=0.0_dp
do l=lmin,lmax
    SGS=SGS+LEGP(l,0)*ACF(l,0)
end do
do m=1,lmax
    do l=max(m,lmin),lmax
        SGS=SGS+LEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
    end do
end do
end


REAL(KIND=dp) FUNCTION RSPH1(ALM,BLM,mu,phi,beta,lmax)

! Computes the radius of the lognormal particle (spherical coordinate
! system).
!
! Copyright by Karri Muinonen. All rights are reserved.
! Version 3.1, 2003 September 12.

use common
implicit none


integer  :: lmax
real(kind=dp),dimension(0:256,0:256) :: ALM,BLM
real(kind=dp) :: SSPH1,mu,phi,beta

RSPH1=exp(SSPH1(ALM,BLM,mu,phi,lmax)-0.5_dp*beta**2)
end



REAL(KIND=dp) FUNCTION SSPH1(ALM,BLM,mu,phi,lmax)

! Computes the logradius in a given direction (spherical coordinate
! system).
!
! Copyright by Karri Muinonen. All rights are reserved.
! Version 3.1, 2003 September 12.

use common
implicit none


real(kind=dp),dimension(0:256,0:256) :: ALM,BLM,PLM
real(kind=dp) :: mu,phi
real(kind=dp),dimension(256) :: SPHI,CPHI
integer  :: l,m,lmax,temp


! Instant return if lmax=0:

if (lmax.eq.0) then
SSPH1=ALM(0,0)
return
endif

! Precomputation of sines, cosines, and associated Legendre functions:
temp=0
call PLMG(PLM,mu,lmax,temp)
do m=1,lmax
call PLMG(PLM,mu,lmax,m)
CPHI(m)=cos(m*phi)
SPHI(m)=sin(m*phi)
end do
PLM(0,0)=1.0_dp

! Sum for the logradius:

SSPH1=0.0_dp
do l=0,lmax
SSPH1=SSPH1+PLM(l,0)*ALM(l,0)
end do
do m=1,lmax
do l=m,lmax
SSPH1=SSPH1+PLM(l,m)*(ALM(l,m)*CPHI(m)+BLM(l,m)*SPHI(m))
end do
end do
end


subroutine PLMG(PLM,x,lmax,m)

! Computes associated Legendre functions from degree l=m
! up to l=lmax.
!
! Copyright by Karri Muinonen. All rights are reserved.
! Version 3.1, 2003 September 12.

use common
implicit none

integer  :: lmax,l,m,temp
real(kind=dp),dimension(0:256,0:256) :: PLM
real(kind=dp) :: FACTRL,x
complex(kind=dp),dimension(0:256,0:256,-2:2) :: PLMM
complex(kind=dp) :: i

i=cmplx(0.0_dp,1.0_dp,kind=dp)

! Check degree, orders, and argument:

if (lmax.lt.0) stop 'Trouble in PLMG: degree negative.'

if (m.gt.lmax .or. m.lt.0) stop 'Trouble in PLMG: order out of range.'

if (abs(x).gt.1.0_dp) stop 'Trouble in PLMG: argument out of range.'

! Compute associated Legendre functions with the help of
! the generalized spherical functions:

temp=0
call PLMMG(PLMM,x,lmax,m,temp)

do l=m,lmax
PLM(l,m)=real(i**m*sqrt(FACTRL(l+m)/FACTRL(l-m))*PLMM(l,m,0))
end do
end subroutine PLMG



subroutine PLMMG(PLMM,x,lmax,m1,m2)

! Computes generalized spherical functions from degree max(abs(m1),abs(m2))
! up to lmax.
!
! Copyright by Karri Muinonen. All rights are reserved.
! Version 3.1, 2003 September 12.

use common
implicit none


integer  :: lmax,l,m0,m1,m2,m12,p12
real(kind=dp) :: FACTRL,x
complex(kind=dp),dimension(0:256,0:256,-2:2) :: PLMM
complex(kind=dp) :: i

i=cmplx(0.0_dp,1.0_dp,kind=dp)

! Check degree, orders, and argument:

if (lmax.lt.0) stop 'Trouble in PLMMG: degree negative.'

if (abs(m1).gt.lmax .or. abs(m2).gt.min(2,lmax) .or. m1.lt.0) stop 'Trouble in PLMMG: order out of range.'

if (abs(x).gt.1.0_dp) stop 'Trouble in PLMMG: argument out of range.'

! Compute generalized spherical functions:

m0=max(abs(m1),abs(m2))
m12=abs(m1-m2)
p12=abs(m1+m2)

if (m0.gt.0) then

if (m12.ne.0 .and. p12.ne.0) then
PLMM(m0,m1,m2)=(-i)**m12/2.0_dp**m0*sqrt(FACTRL(2*m0)/(FACTRL(m12)*FACTRL(p12))*(1.0_dp-x)**m12*(1.0_dp+x)**p12)
elseif (m12.eq.0) then
PLMM(m0,m1,m2)=1.0_dp/2.0_dp**m0*sqrt(FACTRL(2*m0)/FACTRL(p12)*(1.0_dp+x)**p12)
else
PLMM(m0,m1,m2)=(-i)**m12/2.0_dp**m0*sqrt(FACTRL(2*m0)/FACTRL(m12)*(1.0_dp-x)**m12)
endif

if (m0.eq.lmax) return

PLMM(m0+1,m1,m2)=(2*m0+1)*(m0*(m0+1)*x-m1*m2)*PLMM(m0,m1,m2)/(m0*sqrt(dble((m0+1)**2-m2**2))*sqrt(dble((m0+1)**2-m1**2)))

if (m0+1.eq.lmax) return

do l=m0+1,lmax-1
PLMM(l+1,m1,m2)=((2*l+1)*(l*(l+1)*x-m1*m2)*PLMM(l,m1,m2)&
    -(l+1)*sqrt(dble((l**2-m1**2)*(l**2-m2**2)))*PLMM(l-1,m1,m2))/&
    (l*sqrt(dble(((l+1)**2-m1**2)*((l+1)**2-m2**2))))
end do

else

PLMM(0,0,0)=1.0_dp
if (lmax.eq.0) return
PLMM(1,0,0)=x
if (lmax.eq.1) return

do l=m0+1,lmax-1
PLMM(l+1,0,0)=((2*l+1)*x*PLMM(l,0,0)-l*PLMM(l-1,0,0))/(l+1)
end do

endif
end subroutine PLMMG



