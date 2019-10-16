module ellips
! Ellipsoid (defined octant-wise)
!
! RELLTD:   discrete triangle representation
! ARELLTD:  triangle areas
! RELL:     radial distance

use parameters
use voper

contains


subroutine RELLTD(X,N,MU,PHI,AX,rmax,IT,nnod,ntri)

! Discrete triangle representation for an ellipsoid.

! Author: Johanna Torppa (based on code by Karri Muinonen)
! Version: 2015-07-23

!use constants,only:rk
implicit none
integer::IT(260000,3),nnod,ntri,j1,j2
real(kind=dp)::X(130000,3),N(260000,3),MU(130000),PHI(130000),AX(6),&
		 X1(3),X2(3),X3(3),rmax,r,nu

! Node coordinates:

rmax=0.0d0
do j1=1,nnod
        r=RELL(AX,MU(j1),PHI(j1))
        nu=sqrt(1.0d0-MU(j1)**2)
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
    call VPRO(X3,X1,X2)
    r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
    do j2=1,3
         N(j1,j2)=X3(j2)/r
	end do
end do
end


subroutine ARELLTD(AR,MU,PHI,AX,rmax,IT,nnod,ntri)

! Triangle areas for an ellipsoid.

! Author: Johanna Torppa (based on code by Karri Muinonen)
! Version: 2015-07-23

!use constants,only:rk
implicit none
integer::IT(260000,3),nnod,ntri,j1,j2
real(kind=dp)::AR(260000),X(130000,3),MU(130000),PHI(130000),&
		 AX(6),X1(3),X2(3),X3(3),rmax,r,nu

! Node coordinates:

rmax=0.0d0
do j1=1,nnod
        r=RELL(AX,MU(j1),PHI(j1))
        nu=sqrt(1.0d0-MU(j1)**2)
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
    call VPRO(X3,X1,X2)
	AR(j1)=0.5d0*(X3(1)**2+X3(2)**2+X3(3)**2)
end do
end



function RELL(AX,mu,phi)

! Computes the radius of an ellipsoid (spherical coordinate system).

! Author: Johanna Torppa (based on code by Karri Muinonen)
! Version: 2015-07-23
!use constants,only:rk,pi
implicit none
real(kind=dp)::RELL,AX(6),mu,phi,a,b,c

if (phi.ge.0.0d0 .and. phi.lt.0.5d0*pi) then
        a=AX(1)
        b=AX(2)
elseif (phi.ge.0.5d0*pi .and. phi.lt.pi) then
        a=AX(3)
        b=AX(2)
elseif (phi.ge.pi .and. phi.lt.1.5d0*pi) then
        a=AX(3)
        b=AX(4)
else
        a=AX(1)
        b=AX(4)
endif
if (mu.ge.0.0d0) then
        c=AX(5)
else
        c=AX(6)
endif

RELL=a*b*c/sqrt(c**2*(1.0d0-mu**2)*(b**2*cos(phi)**2+&
     a**2*sin(phi)**2)+a**2*b**2*mu**2)
end



function GSDE(AX,mu,phi)

! Computes the Gaussian surface density for an ellipsoid (spherical 
! normal coordinate system). Note the discontinuity for octant 
! ellipsoids.
!
! Author: Johanna Torppa (based on code by Karri Muinonen)
! Version: 2015-07-23

!use constants,only:pi,rk
implicit none
real(kind=dp)::GSDE,AX(6),mu,phi,mu2,a,b,c

if (phi>=0.0d0 .and. phi<0.5d0*pi) then
        a=AX(1)
        b=AX(2)
elseif (phi>=0.5d0*pi .and. phi<pi) then
        a=AX(3)
        b=AX(2)
elseif (phi>=pi .and. phi<1.5d0*pi) then
        a=AX(3)
        b=AX(4)
else
        a=AX(1)
        b=AX(4)
endif
if (mu>=0.0d0) then
        c=AX(5)
else
        c=AX(6)
endif
       
mu2=mu**2
GSDE=(a*b*c)**2/((1.0d0-mu2)*((a*cos(phi))**2+(b*sin(phi))**2)+mu2*c**2)**2
end

end module ellips
