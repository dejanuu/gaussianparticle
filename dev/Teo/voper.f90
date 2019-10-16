module voper

! Vector calculus:
! VROTEU:  vector rotation using Euler angles
! VROTEUT: transpose of vector rotation using Euler angles
! VROTX:   vector rotation about the x-axis
! VROTY:   vector rotation about the y-axis
! VROTZ:   vector rotation about the z-axis
! VPRO:    vector product
! SPRO:    scalar product
! VDIF:    vector difference
! VDIFN:   normalized vector difference
! dot:     Dot product(duplicate, should be refactored)
! crossproduct: Vector product (same deal)
! modlo:   Modulo function
! vector:  Vector from 2 points (TODO: make it work)
! Required kind=dps:
!
use parameters      ! Parameters
implicit none


contains

subroutine VROTEU(X,CA,SA)

! Vector rotation using Euler angles. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(inout) :: X(3),CA(3),SA(3)

call VROTZ(X,CA(1),SA(1))
call VROTY(X,CA(2),SA(2))
call VROTZ(X,CA(3),SA(3))
end subroutine VROTEU



subroutine VROTEUT(X,CA,SA)

! Transpose of vector rotation using Euler angles. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

implicit none
real(kind=dp),intent(inout) :: X(3)
real(kind=dp),intent(in) :: CA(3),SA(3)

call VROTZ(X,CA(3),-SA(3))
call VROTY(X,CA(2),-SA(2))
call VROTZ(X,CA(1),-SA(1))
end subroutine VROTEUT



subroutine VROTX(X,c,s)

! Vector rotation about the x-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(in) :: c,s
real(kind=dp),intent(inout) :: X(3)
real(kind=dp) :: q

q   = c*X(2)+s*X(3)
X(3)=-s*X(2)+c*X(3)
X(2)=q
end subroutine VROTX



subroutine VROTY(X,c,s)

! Vector rotation about the y-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

 implicit none
 real(kind=dp),intent(in) :: c,s
 real(kind=dp),intent(inout) :: X(3)
 real(kind=dp) :: q

 q   = c*X(3)+s*X(1)
 X(1)=-s*X(3)+c*X(1)
 X(3)=q
end subroutine VROTY



subroutine VROTZ(X,c,s)

! Vector rotation about the z-axis. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(in) :: c,s
real(kind=dp),intent(inout) :: X(3)
real(kind=dp) :: q

q   = c*X(1)+s*X(2)
X(2)=-s*X(1)+c*X(2)
X(1)=q
end subroutine VROTZ



subroutine VPRO(XY,X,Y)

! Vector product. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(inout) :: XY(3),X(3),Y(3)

XY(1)=X(2)*Y(3)-X(3)*Y(2)
XY(2)=X(3)*Y(1)-X(1)*Y(3)
XY(3)=X(1)*Y(2)-X(2)*Y(1)    
end subroutine VPRO



subroutine SPRO(XY,X,Y)

! Scalar product. Version 2003-11-07.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
real(kind=dp),intent(inout) :: XY,X(3),Y(3)

XY=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)    
end subroutine SPRO



subroutine VDIF(XY,X,Y)

! Difference of two vectors. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

implicit none
real(kind=dp),intent(in) :: X(3),Y(3)
real(kind=dp),intent(inout) :: XY(3)
integer :: j1

do j1=1,3
 XY(j1)=X(j1)-Y(j1)
end do
end subroutine VDIF



subroutine VDIFN(XY,dxy,X,Y)

! Normalized difference of two vectors. Version 2003-11-07.
!
! Copyright (C) 2003 Karri Muinonen

implicit none
real(kind=dp),intent(in) :: X(3),Y(3)
real(kind=dp),intent(inout) :: dxy,XY(3)
integer :: j1

dxy=0.0d0
do j1=1,3
 XY(j1)=X(j1)-Y(j1)
 dxy=dxy+XY(j1)**2
end do
if (dxy.eq.0.0d0) &
stop 'Trouble in VDIFN: zero vector.'
dxy=sqrt(dxy)
do j1=1,3
 XY(j1)=XY(j1)/dxy
end do
end subroutine VDIFN

! From here till end: Torppa's code for Minkowski.

! The cross product of two vectors

subroutine crossproduct(avect,bvect,crossvect)

! use constants,only:rk
implicit none
real(kind=dp)::avect(3),bvect(3),crossvect(3)

crossvect(1)=avect(2)*bvect(3)-avect(3)*bvect(2)
crossvect(2)=avect(3)*bvect(1)-avect(1)*bvect(3)
crossvect(3)=avect(1)*bvect(2)-avect(2)*bvect(1)
end



! The dot product of two vectors
function dot(avect,bvect)

! use constants,only:rk
implicit none
real(kind=dp),intent(in)::avect(3),bvect(3)
real(kind=dp)::dot
integer n

dot=0.
do n=1,3
    dot=dot+avect(n)*bvect(n)
end do
end function dot



! The modulo function
function modlo(index,modsize)
! use constants,only:rk
implicit none
integer index,modsize,modlo

if (index.gt.modsize) then
        modlo=index-modsize
else
        if (index.lt.1) then
             modlo=index+modsize
        else
              modlo=index
        end if
end if
end function modlo



!A procedure for constructing a vector from two points
subroutine mvector(start,fin,vect,x,y,z)

implicit none
integer::fin,start
real(kind=dp)::vect(3),x(0:nmax),y(0:nmax),z(0:nmax)

vect(1)=x(fin)-x(start)
vect(2)=y(fin)-y(start)
vect(3)=z(fin)-z(start)
end subroutine

end module voper
