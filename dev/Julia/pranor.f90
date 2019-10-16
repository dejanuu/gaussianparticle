
! Particle orientation
!
! PRANOR:  random orientation using Euler angles




subroutine pranor(XN,NT,XN0,NT0,nnod,ntri)

! Random particle orientation using Euler angles.
!
! Version: 2015-02-05
! Author: Karri Muinonen

implicit none

integer, parameter ::                              &
sp = kind(1.0),                                  &
dp = selected_real_kind(2*precision(1.0_sp)),    &
qp = selected_real_kind(2*precision(1.0_dp))
integer  :: nnod,ntri,j1,j2
real(kind=dp) :: gamma,alpha,pi,ran2
real(kind=dp), dimension(130000,3) :: XN, XN0
real(kind=dp), dimension(260000,3) :: NT, NT0
real(kind=dp), dimension(3) :: X, CEU, SEU



pi=4.0_dp*atan(1.0_dp)

! Euler angles in the laboratory reference frame (K) for
! expressing a vector given in K in the particle reference
! frame (K'):

call random_number(ran2)
gamma=2.0_dp*pi*ran2
alpha=2.0_dp*pi*ran2
CEU(1)=cos(gamma)
CEU(2)=1.0_dp-2.0_dp*ran2
CEU(3)=cos(alpha)
SEU(1)=sin(gamma)
SEU(2)=sqrt(1.0_dp-CEU(2)**2)
SEU(3)=sin(alpha)

! Nodes:

do j1 = 1, nnod
    do j2 = 1, 3
        X(j2)=XN0(j1,j2)
    end do
    call vproteut(X,CEU,SEU)
    do j2 = 1, 3
        XN(j1,j2)=X(j2)
   end do
end do

! Normals:

do j1 = 1, ntri
    do j2 = 1, 3
        X(j2)=NT0(j1,j2)
    end do
    call vproteut(X,CEU,SEU)
    do j2 = 1, 3
        NT(j1,j2)=X(j2)
    end do
end do
end subroutine pranor










