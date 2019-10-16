MODULE SIRISRAY

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants

  public
  
  interface raypart0
    module procedure raypart0r, raypart0c
  end interface
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation from the particle to the ray coordinate system.
subroutine partray0(X,mu,nu,cphi,sphi)
       
  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), intent(in) :: mu, nu, cphi, sphi 
  real(kind=dp):: q1, q2, q3 

  q1 = X(1)
  q2 = X(2)
  q3 = X(3)
  X(1) = q1*mu*cphi+q2*mu*sphi-q3*nu
  X(2) = -q1*sphi+q2*cphi
  X(3) = q1*nu*cphi+q2*nu*sphi+q3*mu

end subroutine partray0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns the Mueller matrix, propagation direction, parallel and 
! perpendicular coordinate axes, ray location, and the chord number.
subroutine raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,NISBOX, &
           F,KE,KF,HL,HR,X,MA,pout,nis,jbox)

  real(kind=dp), dimension(4,4,50), intent(in) :: FBOX
  real(kind=dp), dimension(3,50), intent(in) :: KEBOX, KFBOX
  complex(kind=dp), dimension(3,50), intent(in) :: HLBOX, HRBOX 
  real(kind=dp), dimension(3,50), intent(in) :: XBOX
  real(kind=dp), dimension(2,50), intent(in) :: MABOX
  integer, dimension(50), intent(in) :: PBOX, NISBOX 
  real(kind=dp), dimension(4,4), intent(out) :: F
  real(kind=dp), dimension(3), intent(out) :: KE, KF
  complex(kind=dp), dimension(3), intent(out) :: HL, HR
  real (kind=dp), dimension(3), intent(out) :: X
  real (kind=dp), dimension(2), intent(out) :: MA
  integer, intent(out) :: pout, nis
  integer, intent(in) :: jbox
  integer :: j1, j2

  KE = KEBOX(1:3,jbox)
  KF = KFBOX(1:3,jbox) ! onko KFBOX vai KEBOX (Karrin koodissa); testasin, ei nayta vaikuttavan tuloksiin...
  HL = HLBOX(1:3,jbox)
  HR = HRBOX(1:3,jbox)
  X = XBOX(1:3,jbox)    
  MA = MABOX(1:2,jbox)

  do j1 = 1, 4
    do j2 = 1, 4
      F(j1,j2) = FBOX(j1,j2,jbox)
    end do
  end do

  pout = PBOX(jbox)
  nis = NISBOX(jbox)

end subroutine raygetd1c  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Initial unit Mueller matrix, propagation direction, and parallel 
! and perpendicular coordinate axes.
subroutine rayinic(F,KE,KF,HL,HR)

  real(kind=dp), dimension(4,4), intent(out) :: F
  real(kind=dp), dimension(3), intent(out) :: KE, KF
  complex(kind=dp), dimension(3), intent(out) :: HL, HR

  KE = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
  KF = (/ 0.0_dp, 0.0_dp, 1.0_dp /)
  HL = cmplx((/ 1.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp,dp)
  HR = cmplx((/ 0.0_dp, 1.0_dp, 0.0_dp /), 0.0_dp, dp)
  F = reshape((/ 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
    0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/4, 4/))

end subroutine rayinic


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation from the ray to the particle coordinate system.
! Implementation for interface raypart0
subroutine raypart0r(X,mu,nu,cphi,sphi)

  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), intent(in) :: mu, nu, cphi, sphi
  real(kind=dp) :: q1, q2, q3

  q1 = X(1)
  q2 = X(2)
  q3 = X(3)
  X(1) = q1*mu*cphi-q2*sphi+q3*nu*cphi
  X(2) = q1*mu*sphi+q2*cphi+q3*nu*sphi
  X(3) = -q1*nu+q3*mu

end subroutine raypart0r


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation from the ray to the particle coordinate system (complex).
! Implementation for interface raypart0
subroutine raypart0c(X,mu,nu,cphi,sphi)

  complex(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), intent(in) :: mu, nu, cphi, sphi
  complex(kind=dp) :: q1, q2, q3

  q1 = X(1)
  q2 = X(2)
  q3 = X(3)
  X(1) = q1*mu*cphi-q2*sphi+q3*nu*cphi
  X(2) = q1*mu*sphi+q2*cphi+q3*nu*sphi
  X(3) = -q1*nu+q3*mu
       
end subroutine raypart0c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stores the Mueller matrix, propagation direction, parallel and 
! perpendicular coordinate axes, ray location, chord number,
! and the triangle number.
subroutine rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX, &
                XBOX,MABOX,PBOX,NISBOX,F,KE,KF,HL,HR,X,MA, &
                pout,nis,jbox)

  real(kind=dp), dimension(4,4,50), intent(inout) :: FBOX
  real(kind=dp), dimension(3,50), intent(inout) :: KEBOX, KFBOX
  complex(kind=dp), dimension(3,50), intent(inout) :: HLBOX, HRBOX  
  real(kind=dp), dimension(3,50), intent(inout) :: XBOX
  real(kind=dp), dimension(2,50), intent(inout) :: MABOX
  integer, dimension(50), intent(inout) :: PBOX, NISBOX
  real(kind=dp), dimension(4,4), intent(in) :: F
  real(kind=dp), dimension(3), intent(in) :: KE, KF
  complex(kind=dp), dimension(3), intent(in) :: HL, HR
  real(kind=dp), dimension(3), intent(in) :: X
  real(kind=dp), dimension(2), intent(in) :: MA
  integer, intent(in) :: pout, nis, jbox
  integer :: j1, j2

  do j1 = 1, 3
    KEBOX(j1,jbox) = KE(j1)
    KFBOX(j1,jbox) = KF(j1)
    HLBOX(j1,jbox) = HL(j1)
    HRBOX(j1,jbox) = HR(j1)
    XBOX(j1,jbox) = X(j1)
  end do

  MABOX(1,jbox) = MA(1)
  MABOX(2,jbox) = MA(2)

  do j1 = 1, 4
    do j2 = 1, 4
      FBOX(j1,j2,jbox)=F(j1,j2)
    end do
  end do

  PBOX(jbox) = pout
  NISBOX(jbox) = nis

end subroutine rayputd1c


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stores the ray into the scattering phase matrix.
subroutine raysca(S,F,k)

  real(kind=dp), dimension(361,4,4), intent(inout) :: S   
  real (kind=dp), dimension(4,4), intent(in) :: F
  integer, intent(in) :: k
  integer :: j1, j2

  do j1 = 1, 4
    do j2 = 1,4
      S(k,j1,j2) = S(k,j1,j2)+F(j1,j2)
    end do
  end do

end subroutine raysca


END MODULE SIRISRAY
