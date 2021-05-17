MODULE SIRISRADTRANS

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  use sirismath
  use sirismaterial
  use sirisnumint
  use sirisray

  public
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Attenuates the Mueller matrix, and computes the absorption 
! cross section. 
subroutine absorb(F,qabs,len,abscf)
      
  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), intent(inout) :: qabs
  real(kind=dp), intent(in) :: len, abscf
  real(kind=dp) :: attcf

  if (abs(abscf*len) < 20.0_dp) then 
    attcf = exp(abscf*len)
  else
    attcf=0.0_dp
  end if

  qabs = qabs+(1.0_dp-attcf)*F(1,1)

  F = attcf*F

end subroutine absorb


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Attenuates the Mueller matrix, and computes the absorption
! cross section.
subroutine absorbrt(F,qabs,ot)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), intent(inout) :: qabs
  real(kind=dp), intent(in) :: ot
  integer :: j1, j2

  qabs=qabs+(1.0_dp-ot)*F(1,1)
  do j1 = 1, 4
    do j2 = 1, 4
      F(j1,j2)=ot*F(j1,j2)
    end do
  end do

end subroutine absorbrt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Maps the cosines of scattering angles to random numbers for
! the spline scattering phase function of the spline scattering
! phase matrix. Version 2008-04-10.
subroutine cspspliv(CSRN,XP,YP,YP2,nrn,np,ncm)

  real(kind=dp), dimension(0:), intent(out) :: CSRN
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  integer, intent(in) :: nrn, np, ncm
  integer :: j1
  real(kind=dp) :: rn
  real(kind=dp), parameter :: cstol=1.0e-6

  CSRN(0)=-1.0_dp
  do j1 = 1, nrn-1
    rn=dble(j1)/nrn
    CSRN(j1)=CSVRTBIS(XP,YP,YP2,rn,-1.0_dp,1.0_dp,cstol,np,ncm)
  end do
  CSRN(nrn)=1.0_dp

end subroutine cspspliv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real(kind=dp) function CSVRTBIS(XP,YP,YP2,rn,x1,x2,xacc,np,ncm)

  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  real(kind=dp), intent(in) ::rn, x1, x2, xacc
  integer, intent(in) :: np, ncm
  integer :: j
  integer, parameter :: JMAX=40
  real(kind=dp) :: xmid, dx, fmid, f

  fmid=PCDFSPLIV(x2,XP,YP,YP2,np,ncm)-rn
  f=PCDFSPLIV(x1,XP,YP,YP2,np,ncm)-rn

  if(f*fmid >= 0.0_dp) stop 'Trouble in CSVRTBIS: root must be bracketed.'

  if(f < 0.0_dp) then
    CSVRTBIS=x1
    dx=x2-x1
  else
    CSVRTBIS=x2
    dx=x1-x2
  endif
  do j = 1, JMAX
    dx=dx*0.5_dp
    xmid=CSVRTBIS+dx
    fmid=PCDFSPLIV(xmid,XP,YP,YP2,np,ncm)-rn
    if(fmid <= 0.0_dp) CSVRTBIS=xmid
    if(abs(dx) < xacc .or. fmid == 0.0_dp) return
  end do
  stop 'Trouble in CSVRTBIS: too many bisections.'

end function CSVRTBIS


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation of the Mueller matrix from the left (F <- KF).
subroutine frotl(F,c2psi,s2psi)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), intent(in) :: c2psi, s2psi
  real(kind=dp) :: q
  integer :: j1

  do j1 = 1, 4
    q = c2psi*F(2,j1)+s2psi*F(3,j1)
    F(3,j1) = -s2psi*F(2,j1)+c2psi*F(3,j1)
    F(2,j1) = q
  end do

end subroutine frotl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation of the Mueller matrix from the left (F2 =K*F1).
subroutine frotlc(F2,F1,H1,H2,G1,G2)

  real(kind=dp), dimension(4,4), intent(out) :: F2
  real(kind=dp), dimension(4,4), intent(in) :: F1
  complex(kind=dp), dimension(3), intent(in) :: H1, H2, G1, G2
  real(kind=dp), dimension(2,4) :: A      
  real(kind=dp), dimension(4,4) :: K
  complex(kind=dp) :: img
  complex(kind=dp), dimension(4) :: B
  complex(kind=dp), dimension(2,2) :: Q

  img = cmplx(0.0_dp,1.0_dp,dp)

  call proscac(Q(1,1),H1,G1)
  call proscac(Q(1,2),H1,G2)
  call proscac(Q(2,1),H2,G1)
  call proscac(Q(2,2),H2,G2)

  A(1,1) = 0.5_dp*(abs(Q(1,2))**2+abs(Q(1,1))**2)
  A(1,2) = 0.5_dp*(abs(Q(1,2))**2-abs(Q(1,1))**2)
  A(1,3) = real(Q(1,2)*conjg(Q(1,1)),dp)
  A(1,4) = aimag(Q(1,2)*conjg(Q(1,1)))

  A(2,1) = 0.5_dp*(abs(Q(2,2))**2+abs(Q(2,1))**2)
  A(2,2) = 0.5_dp*(abs(Q(2,2))**2-abs(Q(2,1))**2)
  A(2,3) = real(Q(2,2)*conjg(Q(2,1)),dp)
  A(2,4) = aimag(Q(2,2)*conjg(Q(2,1)))

  B(1) = 0.5_dp*(Q(2,2)*conjg(Q(1,2))+Q(2,1)*conjg(Q(1,1)))
  B(2) = 0.5_dp*(Q(2,2)*conjg(Q(1,2))-Q(2,1)*conjg(Q(1,1)))
  B(3) = 0.5_dp*(Q(2,2)*conjg(Q(1,1))+Q(2,1)*conjg(Q(1,2)))
  B(4) = -0.5_dp*img*(Q(2,2)*conjg(Q(1,1))-Q(2,1)*conjg(Q(1,2)))

  K(1,1) = A(1,1)+A(2,1)
  K(1,2) = A(1,2)+A(2,2)
  K(1,3) = A(1,3)+A(2,3)
  K(1,4) = A(1,4)+A(2,4)

  K(2,1) = -A(1,1)+A(2,1)
  K(2,2) = -A(1,2)+A(2,2)
  K(2,3) = -A(1,3)+A(2,3)
  K(2,4) = -A(1,4)+A(2,4)

  K(3,1) = real(2.0_dp*B(1),dp)
  K(3,2) = real(2.0_dp*B(2),dp)
  K(3,3) = real(2.0_dp*B(3),dp)
  K(3,4) = real(2.0_dp*B(4),dp)

  K(4,1) = -aimag(2.0_dp*B(1))
  K(4,2) = -aimag(2.0_dp*B(2))
  K(4,3) = -aimag(2.0_dp*B(3))
  K(4,4) = -aimag(2.0_dp*B(4))

  call MMMM(F2,K,F1)

end subroutine frotlc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation of the Mueller matrix from the left (F2 =K*F1).
subroutine frotlc2l(F2,F1,H1,H2,G1,G2)
    
  real(kind=dp), dimension(4,4), intent(out) :: F2
  real(kind=dp), dimension(4,4), intent(in) :: F1
  complex(kind=dp), dimension(3), intent(in) :: H1, H2, G1, G2
  integer :: j1, j2
  real(kind=dp) :: norm
  real(kind=dp), dimension(4,4) :: K
  complex(kind=dp) :: D
  complex(kind=dp), dimension(2) :: U, V
  complex(kind=dp), dimension(2,2) :: Q

  call PROSCAC(Q(1,1),H1,G1)
  call PROSCAC(Q(1,2),H1,G2)
  call PROSCAC(Q(2,1),H2,G1)
  call PROSCAC(Q(2,2),H2,G2)

  D=Q(1,1)*Q(2,2)-Q(1,2)*Q(2,1)
  U(1)=Q(2,2)/D
  U(2)=-Q(1,2)/D
  V(1)=-Q(2,1)/D
  V(2)=Q(1,1)/D

  K(1,1)=.5_dp*(abs(U(1))**2+abs(V(1))**2+abs(U(2))**2+abs(V(2))**2)
  K(1,2)=.5_dp*(abs(U(1))**2-abs(V(1))**2+abs(U(2))**2-abs(V(2))**2)
  K(1,3)=real(U(1)*conjg(V(1)),dp)+real(U(2)*conjg(V(2)),dp)
  K(1,4)=aimag(U(1)*conjg(V(1)))+aimag(U(2)*conjg(V(2)))

  K(2,1)=.5_dp*(abs(U(1))**2+abs(V(1))**2-abs(U(2))**2-abs(V(2))**2)
  K(2,2)=.5_dp*(abs(U(1))**2-abs(V(1))**2-abs(U(2))**2+abs(V(2))**2)
  K(2,3)=real(U(1)*conjg(V(1)),dp)-real(U(2)*conjg(V(2)),dp)
  K(2,4)=aimag(U(1)*conjg(V(1)))-aimag(U(2)*conjg(V(2)))

  K(3,1)=real(U(1)*conjg(U(2)),dp)+real(V(1)*conjg(V(2)),dp)
  K(3,2)=real(U(1)*conjg(U(2)),dp)-real(V(1)*conjg(V(2)),dp)
  K(3,3)=real(U(1)*conjg(V(2)),dp)+real(U(2)*conjg(V(1)),dp)
  K(3,4)=aimag(U(1)*conjg(V(2)))+aimag(U(2)*conjg(V(1)))

  K(4,1)=-aimag(U(1)*conjg(U(2)))-aimag(V(1)*conjg(V(2)))
  K(4,2)=-aimag(U(1)*conjg(U(2)))+aimag(V(1)*conjg(V(2)))
  K(4,3)=-aimag(U(1)*conjg(V(2)))+aimag(U(2)*conjg(V(1)))
  K(4,4)=real(U(1)*conjg(V(2)),dp)-real(U(2)*conjg(V(1)),dp)

  call mmmm(F2,K,F1)

  norm=F1(1,1)/F2(1,1)
  do j1 = 1, 4
    do j2 = 1, 4
      F2(j2,j1) = norm*F2(j2,j1)
    end do
  end do

end subroutine frotlc2l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rotation of the Mueller matrix from the right (F <- FK).
subroutine frotr(F,c2psi,s2psi)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), intent(in) :: c2psi, s2psi
  real(kind=dp) :: q
  integer :: j1

  do j1 = 1, 4
    q = c2psi*F(j1,2)-s2psi*F(j1,3)
    F(j1,3) = s2psi*F(j1,2)+c2psi*F(j1,3)
    F(j1,2) = q
  end do

end subroutine frotr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For an incident Mueller matrix and ray coordinate system (subscript 1), 
! computes the reflected (subscript 1) and refracted (subscript 2) 
! Mueller matrices and ray coordinate systems.
subroutine incide(F1,KE1,KF1,HL1,HR1,MA1,F2,KE2,KF2,HL2,HR2,MA2, &
                        N,m1,m2,totref)

  real(kind=dp), dimension(4,4), intent(inout) :: F1
  real(kind=dp), dimension(3), intent(inout) :: KE1, KF1
  complex(kind=dp), dimension(3), intent(inout) :: HL1, HR1
  real(kind=dp), dimension(2), intent(inout) :: MA1
  real(kind=dp), dimension(4,4), intent(out) :: F2
  real(kind=dp), dimension(3), intent(inout) :: KE2, KF2
  complex(kind=dp), dimension(3), intent(out) :: HL2, HR2
  real(kind=dp), dimension(2), intent(inout) :: MA2  
  real(kind=dp), dimension(3), intent(in) :: N
  complex(kind=dp), intent(in) :: m1, m2
  integer, intent(out) :: totref
  integer :: j1
  real(kind=dp):: ti, ci, ct, st, t11, t12, t33, t34, r11, r12, r33, r34, &
    mre, mim, ctl, ctr, nke, nkf, kekf1, kekf2, kekf
  real(kind=dp), dimension(2) :: MAI  
  real(kind=dp), dimension(3) :: KEI, KFI, T1, T2
  real(kind=dp), dimension(4,4) :: FI
  complex(kind=dp) :: m, frl, frr, ftl, ftr, fki, fkt
  complex(kind=dp), dimension(3) :: HLI, HRI, HLI0, HRI0, &
    NC, KC1, KC2, KCI




  ! Relative refractive index:
  m = m2/m1
  mre = real(m,dp)
  mim = aimag(m)

  ! Incident ray:
  MAI = MA1
  KEI = KE1
  KFI = KF1
  HLI0 = HL1
  HRI0 = HR1
  FI = F1

  call prosca(nke,N,KEI)
  call prosca(nkf,N,KFI)
  call prosca(kekf,KEI,KFI)

  do j1 = 1, 3
    NC(j1) = cmplx(N(j1),0.0_dp,dp)
    KCI(j1) = cmplx(MAI(1)*KEI(j1),MAI(2)*KFI(j1),dp)
  end do
  call refindapp(MA1,m1,kekf)

  ! OBLIQUE INCIDENCE:
  if (nke > -0.9999999_dp) then

    ! Two unit tangent vectors at the surface:
    ti = sqrt(1.0_dp - nke**2)
    T2 = (KEI-nke*N)/ti
    call provec(T1,T2,N)

    ! Complex unit vectors in the plane of incidence and
    ! rotation of the input Mueller matrix
    call provecn(HRI,KCI,NC)
    call provecn(HLI,HRI,KCI)
    call frotlc(F1,FI,HRI,HLI,HRI0,HLI0)

    ! Unit direction vectors for the reflected and refracted ray:
    call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T2, &
      fki,fkt,nke,nkf,m1,m2,st)

    ! Complex unit vectors of the reflected rays:
    KC1 = cmplx(MA1(1)*KE1, MA1(2)*KF1,dp)

    call proveccn(HR1,KC1,NC)
    call proveccn(HL1,HR1,KC1)

    ! Total Fresnel reflection, new reflected Mueller matrix:
    if(st >= 1.0_dp) then
      ci = -nke
      ct = sqrt(st**2-1.0_dp)
      frl = cmplx(mre*ci,-ct,dp)/cmplx(mre*ci,ct,dp)
      frr = cmplx(ci,-mre*ct,dp)/cmplx(ci,mre*ct,dp)
      r11 = 1.0_dp
      r12 = 0.0_dp
      r33 = real(frl*conjg(frr),dp)
      r34 = aimag(frl*conjg(frr))
      call reflect(F1,r11,r12,r33,r34)
      totref = 1
      return
    end if

    ! Ordinary Fresnel refraction and reflection, 
    ! new refracted and reflected Mueller matrices:
    frr = RR(fki,fkt)
    frl = RL(m1,m2,fki,fkt)
    ftr = TR(fki,fkt)
    ftl = TL(m1,m2,fki,fkt)
    ctl = abs(real(conjg(m)*fkt/(m*fki)))   
    ctr = abs(real(fkt/fki))

    r11 = 0.5_dp*(abs(frl)**2+abs(frr)**2)
    r12 = 0.5_dp*(abs(frl)**2-abs(frr)**2)
    r33 = real(frl*conjg(frr),dp)
    r34 = aimag(frl*conjg(frr))

    t11 = 1.0_dp-r11
    t12 = 0.5_dp*(ctl*abs(ftl)**2-ctr*abs(ftr)**2)
    t33 = sqrt(ctl*ctr)*real(ftl*conjg(ftr),dp)
    t34 = sqrt(ctl*ctr)*aimag(ftl*conjg(ftr))

    call refract(F1,F2,t11,t12,t33,t34)
    call reflect(F1,r11,r12,r33,r34)

    ! Complex unit vectors of the refracted rays:
    KC2 = cmplx(MA2(1)*KE2, MA2(2)*KF2,dp)

    call proveccn(HR2,KC2,NC)
    call proveccn(HL2,HR2,KC2)
    totref = -1

  ! NORMAL INCIDENCE:
  else

    ! Two unit tangent vectors at the surface:
    T1 = real(HRI0,dp)
    call provecn(T2,T1,KEI)
    call provec(T1,N,T2)

    ! Complex unit vectors in the plane of incidence and
    ! rotation of the input Mueller matrix:
    call proveccn(HRI,KCI,NC)   
    call proveccn(HLI,HRI,KCI)
    call frotlc(F1,FI,HRI,HLI,HRI0,HLI0)

    ! Unit direction vectors for the reflected and refracted ray:
    call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T2, &
       fki,fkt,nke,nkf,m1,m2,st)

    ! Complex unit vectors of the reflected rays:
    call prosca(kekf1,KE1,KF1)
    if (abs(kekf1) < 0.9999999_dp) then
      KC1 = cmplx(MA1(1)*KE1, MA1(2)*KF1,dp)
      call proveccn(HR1,KC1,NC)
      call proveccn(HL1,HR1,KC1)
    else
      HL1 = -HLI0
      HR1 = HRI0
    end if

    ! Complex unit vectors of the refracted rays:
    call prosca(kekf2,KE2,KF2)
    if (abs(kekf2) < 0.9999999_dp) then
      KC2 = cmplx(MA2(1)*KE2, MA2(2)*KF2,dp)
      call proveccn(HR2,KC2,NC)
      call proveccn(HL2,HR2,KC2)
    else
      HL2 = HLI0
      HR2 = HRI0
    end if

    ! Ordinary Fresnel refraction and reflection, 
    ! new refracted and reflected Mueller matrices:
    r11 = abs(RL(m1,m2,fki,fkt))**2
    r12 = 0.0_dp
    r33 = -r11
    r34 = 0.0_dp
    t11 = 1.0_dp-r11
    t12 = 0.0_dp
    t33 = t11
    t34 = 0.0_dp

    call refract(F1,F2,t11,t12,t33,t34)
    call reflect(F1,r11,r12,r33,r34)
    totref = -1
  end if 

end subroutine incide


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For an incident Mueller matrix and ray coordinate system (subscript 1), 
! computes the reflected (subscript 1) and refracted (subscript 2) 
! Mueller matrices and ray coordinate systems.
subroutine incide2l(F1,KE1,KF1,HL1,HR1,MA1,F2,KE2,KF2,HL2,HR2,MA2, &
                        N,m1,m2,totref,preventTR_in)

  real(kind=dp), dimension(4,4), intent(inout) :: F1
  real(kind=dp), dimension(3), intent(inout) :: KE1, KF1
  complex(kind=dp), dimension(3), intent(inout) :: HL1, HR1
  real(kind=dp), dimension(2), intent(inout) :: MA1
  real(kind=dp), dimension(4,4), intent(out) :: F2
  real(kind=dp), dimension(3), intent(inout) :: KE2, KF2
  complex(kind=dp), dimension(3), intent(out) :: HL2, HR2
  real(kind=dp), dimension(2), intent(inout) :: MA2  
  real(kind=dp), dimension(3), intent(in) :: N
  complex(kind=dp), intent(in) :: m1, m2
  integer, intent(out) :: totref
  integer :: j1
  real(kind=dp):: ti, ci, ct, st, t11, t12, t33, t34, r11, r12, r33, r34, &
    mre, mim, ctl, ctr, nke, nkf, kekf1, kekf2, dm, norm1,kekf
  real(kind=dp), parameter :: dmtol=1e-10_dp, ktol=1e-10_dp
  real(kind=dp), dimension(2) :: MAI  
  real(kind=dp), dimension(3) :: KEI, KFI, T1, T2
  real(kind=dp), dimension(4,4) :: FI
  complex(kind=dp) :: m, frl, frr, ftl, ftr, fki, fkt
  complex(kind=dp), dimension(3) :: HLI, HRI, HLI0, HRI0, &
    NC, KC1, KC2, KCI
  logical :: totref0
   integer :: preventTR
   integer,optional :: preventTR_in
    if(present(preventTR_in))then
        preventTR=preventTR_in
    else
        preventTR=0
    endif
totref=0
    dm=abs(m1-m2)

! Total refraction:

    if (dm.lt.dmtol) then
        MA2(1)=MA1(1)
        MA2(2)=MA1(2)

        KE2=KE1
        KF2=KF1
        HL2=HL1
        HR2=HR1

        F2=F1
        totref=2
        return
    end if


! Relative refractive index:

    m = m2/m1
    mre = real(m,dp)
    mim = aimag(m)

! Incident ray:

    MAI = MA1
    KEI = KE1
    KFI = KF1
    HLI0 = HL1
    HRI0 = HR1
    FI = F1

    call prosca(nke,N,KEI)
    call prosca(nkf,N,KFI)
    call prosca(kekf,KEI,KFI)

    if(abs(MAI(2)) .ge. ktol) then
        do j1 = 1, 3
            NC(j1) = cmplx(N(j1),0.0_dp,kind=dp)
            KCI(j1) = cmplx(MAI(1)*KEI(j1),MAI(2)*KFI(j1),kind=dp)
        end do
    else
        do j1 = 1, 3
            NC(j1)=cmplx(N(j1),0.0_dp,kind=dp)
            KCI(j1)=cmplx(MAI(1)*KEI(j1),0.0_dp,kind=dp)
        end do
    end if

! OBLIQUE INCIDENCE:

    if (nke .gt. -0.9999999_dp) then

! Two unit tangent vectors at the surface:

        ti = sqrt(1.0_dp - nke**2)
        norm1 = 0.0_dp
        do j1 = 1, 3
            T2(j1) = (KEI(j1)-nke*N(j1))/ti
            norm1 = norm1 + T2(j1)**2
        end do
        norm1 = sqrt(norm1)
        do j1 = 1, 3
            T2(j1) = T2(j1)/norm1
        end do
        call provecn(T1,T2,N)

! Complex unit vectors in the plane of incidence and
! rotation of the input Mueller matrix:

        call proveccn(HRI,KCI,NC)
        call proveccn(HLI,HRI,KCI)
        call frotlc2l(F1,FI,HLI,HRI,HLI0,HRI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel2l(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                fki,fkt,nke,nkf,m1,m2,st,totref0,preventTR)

! Complex unit vectors of the reflected rays:

        if (abs(MA1(2)) .ge. ktol) then
            KC1 = cmplx(MA1(1)*KE1, MA1(2)*KF1,kind=dp)
        else
            KC1 = cmplx(MA1(1)*KE1, 0.0_dp,kind=dp)
        end if

        call proveccn(HR1,KC1,NC)
        call proveccn(HL1,HR1,KC1)

! Total Fresnel reflection, new reflected Mueller matrix:

        if (abs(st) >= 1.0_dp) then
           ci = -nke
           ct = sqrt(st**2-1.0_dp)
           frl = cmplx(mre*ci,-ct,kind=dp)/cmplx(mre*ci,ct,kind=dp)
           frr = cmplx(ci,-mre*ct,kind=dp)/cmplx(ci,mre*ct,kind=dp)
           r11 = 1.0_dp
           r12 = 0.0_dp
           r33 = real(frl*conjg(frr),dp)
           r34 = aimag(frl*conjg(frr))
           call reflect(F1,r11,r12,r33,r34)
           totref = 1
           return
        end if

! Ordinary Fresnel refraction and reflection, 
! new refracted and reflected Mueller matrices:


        
        
        frr = RR(fki,fkt)
        frl = RL(m1,m2,fki,fkt)
        ftr = TR(fki,fkt)
        ftl = TL(m1,m2,fki,fkt)
        

        r11 = 0.5_dp*(abs(frl)**2+abs(frr)**2)
        r12 = 0.5_dp*(abs(frl)**2-abs(frr)**2)
        r33 = real(frl*conjg(frr),dp)
        r34 = aimag(frl*conjg(frr))


        ctr=(1.0_dp-abs(frr)**2)/abs(ftr)**2
        ctl=(1.0_dp-abs(frl)**2)/abs(ftl)**2

        if(ctr>=0.0_dp .and. ctl>=0.0_dp) then
            t11 = 1.0_dp-r11
            t12 = 0.5_dp*(ctl*abs(ftl)**2-ctr*abs(ftr)**2)
            t33 = sqrt(ctl*ctr)*real(ftl*conjg(ftr),dp)
            t34 = sqrt(ctl*ctr)*aimag(ftl*conjg(ftr))

            call refract(F1,F2,t11,t12,t33,t34)

            if (abs(MA2(2)) .ge. ktol) then
                KC2 = cmplx(MA2(1)*KE2, MA2(2)*KF2,kind=dp)
            else
                KC2 = cmplx(MA2(1)*KE2, 0.0_dp,kind=dp)
            end if

            call proveccn(HR2,KC2,NC)
            call proveccn(HL2,HR2,KC2)
        else
            totref0 = .true.
        endif

        call reflect(F1,r11,r12,r33,r34)


        if(totref0) then
            F2=FI-F1
            totref=3
        endif
    
! NORMAL INCIDENCE:

  else

! Two unit tangent vectors at the surface:

        T1 = real(HRI0,dp)
        call provecn(T2,T1,KEI)
        call provecn(T1,N,T2)

! Complex unit vectors in the plane of incidence and
! rotation of the input Mueller matrix:

        call proveccn(HRI,KCI,NC)   
        call proveccn(HLI,HRI,KCI)
        call frotlc2l(F1,FI,HLI,HRI,HLI0,HRI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel2l(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                 fki,fkt,nke,nkf,m1,m2,st,totref0,preventTR)

! Complex unit vectors of the reflected rays:

        call prosca(kekf1,KE1,KF1)
        if (abs(kekf1) .lt. 0.9999999_dp) then
            KC1 = cmplx(MA1(1)*KE1, MA1(2)*KF1,kind=dp)
            call proveccn(HR1,KC1,NC)
            call proveccn(HL1,HR1,KC1)
        else
            HL1 = -HLI0
            HR1 = HRI0
        end if

! Complex unit vectors of the refracted rays:

        call prosca(kekf2,KE2,KF2)
        if (abs(kekf2) .lt. 0.9999999_dp) then
            KC2 = cmplx(MA2(1)*KE2, MA2(2)*KF2,kind=dp)
            call proveccn(HR2,KC2,NC)
            call proveccn(HL2,HR2,KC2)
        else
          HL2 = HLI0       
          HR2 = HRI0
        end if

! Ordinary Fresnel refraction and reflection, 
! new refracted and reflected Mueller matrices:

        r11 = abs(RL(m1,m2,fki,fkt))**2
        r12 = 0.0_dp
        r33 = -r11
        r34 = 0.0_dp
        t11 = 1.0_dp-r11
        t12 = 0.0_dp
        t33 = t11
        t34 = 0.0_dp

        call refract(F1,F2,t11,t12,t33,t34)
        call reflect(F1,r11,r12,r33,r34)
  end if 

end subroutine incide2l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Approximate Rayleigh scattering algorithm for the
! generation of a new propagation direction and Mueller
! matrix.
subroutine incrt(F,K,EL,ER,CSRN,XP,YP,YP2,np,nrn)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), dimension(3), intent(inout) :: K, EL, ER
  real(kind=dp), dimension(0:), intent(in) :: CSRN
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  integer, intent(in) :: np, nrn
  integer :: j1, j2, mrn
  real(kind=dp) :: nt, c2psi, s2psi, rn, cthe, sthe, phi, nn, ran2
  real(kind=dp), dimension(3) :: K1, EL1, ER1, T1, T2, N
  real(kind=dp), dimension(4,4) :: F1, P

  ! Temporary storage:
  do j1 = 1, 4
    do j2 = 1, 4
      F1(j1,j2)=F(j1,j2)
    end do
  end do
  do j1 = 1, 3
    K1(j1)=K(j1)
    EL1(j1)=EL(j1)
    ER1(j1)=ER(j1)
  end do

  ! New scattering direction (approximately):
  call random_number(ran2)

  rn = nrn*ran2
  mrn = int(rn)
  cthe = (real(mrn+1,dp)-rn)*CSRN(mrn)+(rn-real(mrn,dp))*CSRN(mrn+1)

  sthe = sqrt(1.0_dp-cthe**2)

  call random_number(ran2)
  phi = 2.0_dp*pi*ran2
  do j1 = 1, 3
    K(j1)=sthe*cos(phi)*EL1(j1)+sthe*sin(phi)*ER1(j1)+cthe*K1(j1)
  end do

  ! Auxiliary coordinate system:
  do j1 = 1, 3
    N(j1)=K(j1)-K1(j1)
    T2(j1)=K(j1)+K1(j1)
  end do

  call provec(T1,T2,N)

  nt=0.0_dp
  do j1 = 1, 3
    nt=nt+T1(j1)**2
  end do
  nt=sqrt(nt)
  do j1 = 1, 3
    T1(j1)=T1(j1)/nt
  end do

  ! Rotation of the input Mueller matrix:
  c2psi=(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))**2- (EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))**2
  s2psi=-2.0_dp*(EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))*(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))
  call frotl(F1,c2psi,s2psi)

  ! New direction vectors for the reflected ray:
  do j1 = 1, 3
    ER(j1)=T1(j1)
  end do
  call provec(EL,ER,K)

  call pspliv(P,XP,YP,YP2,acos(cthe),np,0)
  call mmmm(F,P,F1)

  nn=F1(1,1)/F(1,1)
  do j1 = 1, 4
    do j2 = 1, 4
      F(j1,j2)=F(j1,j2)*nn
    end do
  end do

end subroutine incrt


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mueller matrix operation on another Mueller matrix.
subroutine mmmm(M3,M2,M1)
 
  real(kind=dp), dimension(4,4), intent(out) :: M3
  real (kind=dp), dimension(4,4), intent(in) :: M2, M1
  integer  :: j1, j2, j3

  do j3 = 1, 4
    do j2 = 1, 4
      M3(j2,j3) = 0.0_dp
      do j1 = 1, 4
        M3(j2,j3) = M3(j2,j3)+M2(j2,j1)*M1(j1,j3)
      end do
    end do
  end do                

end subroutine mmmm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the cumulative distribution function for the spline
! scattering phase function of the spline scattering phase matrix.
! Version 2008-04-10.
real(kind=dp) function PCDFSPLIV(xi,XP,YP,YP2,np,ncm)

  real(kind=dp), intent(in) :: xi
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  integer, intent(in) :: np, ncm
  integer :: j1, temp
  real(kind=dp), dimension(4,4) :: P
  real(kind=dp), dimension(512) :: X, WX

  call gauleg(-1.0_dp,xi,X,WX,ncm)
  PCDFSPLIV=0.0_dp
  do j1 = 1, ncm
    temp=1
    call PSPLIV(P,XP,YP,YP2,acos(X(j1)),np,temp)
    PCDFSPLIV=PCDFSPLIV+WX(j1)*0.5_dp*P(1,1)
  end do

end function PCDFSPLIV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input of discretized scattering phase matrix.
! Version from multiparticle code
subroutine pmatrix1_multi(P,np)

  real(kind=dp), dimension(0:360,4,4), intent(out) :: P
  integer, intent(in) :: np
  integer :: j1, fu
  real(kind=dp) :: p11, p12, p22, p33, p34, p44, the
  real(kind=dp) :: tmp
  character(len=file_name_length) :: infile

  ! Reset:
  P(:,:,:)=0.0_dp

  ! Scattering phase matrix:
  j1 = command_argument_count()
  if(j1 <= 1) then
    infile = 'pmatrix1.in'
  else
    call get_command_argument(2,infile)
  endif
  write(output_unit, '(A,A,A)') "Reading diffuse particle scattering matrix input from file '", trim(infile), "'"
  open(newunit=fu, file=trim(infile), status='old', action='read')

  do j1 = 0, np
    read(fu,*) the,p11,p12,tmp,tmp, &
      tmp,p22,tmp,tmp, &
      tmp,tmp,p33,p34, &
      tmp,tmp,tmp,p44
    P(j1,1,1)=p11
    P(j1,1,2)=p12
    P(j1,2,1)=p12
    P(j1,2,2)=p22
    P(j1,3,3)=p33
    P(j1,3,4)=p34
    P(j1,4,3)=-p34
    P(j1,4,4)=p44
  end do
  close(fu)

end subroutine pmatrix1_multi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Input of discretized scattering phase matrix.
! Version from single-particle code
subroutine pmatrix1_single(P,np,argnum)

  real(kind=dp), dimension(0:360,4,4), intent(out) :: P
  integer, intent(in) :: np
  integer, intent(in), optional :: argnum
  integer :: j1, fu, argn
  real(kind=dp) :: p11, p12, p22, p33, p34, p44, the
  character(len=file_name_length) :: infile

  ! Reset:
  P(:,:,:)=0.0_dp

  ! Scattering phase matrix:
  j1 = command_argument_count()
  if(present(argnum)) then
    argn = argnum
  else
    argn = 1
  end if
  if(j1 <= argn) then
    infile = 'pmatrix1.in'
  else
    call get_command_argument(argn+1,infile)
  endif
  write(output_unit, '(A,A,A)') "Reading diffuse particle scattering matrix input from file '", trim(infile), "'"
  open(newunit=fu, file=trim(infile), status='old', action='read')

  do j1 = 0, np
    read(fu,*) the,p11,p12,p22,p33,p34,p44
    P(j1,1,1)=p11
    P(j1,1,2)=p12
    P(j1,2,1)=p12
    P(j1,2,2)=p22
    P(j1,3,3)=p33
    P(j1,3,4)=p34
    P(j1,4,3)=-p34
    P(j1,4,4)=p44
  end do
  close(fu)

end subroutine pmatrix1_single


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Rayleigh scattering phase matrix.
! Version 3.1, 2003 September 12.
subroutine praylv(P,cthe)

  real(kind=dp), dimension(4,4), intent(out) :: P
  real(kind=dp), intent(in) :: cthe
  integer :: j1, j2

  ! Reset:
  do j1 = 1, 4
    do j2 = 1, 4
      P(j1,j2)=0.0_dp
    end do
  end do

  ! Phase matrix:
  P(1,1)=0.75_dp*(1.0_dp+cthe**2)
  P(2,2)=P(1,1)
  P(1,2)=-0.75_dp*(1.0_dp-cthe**2)
  P(2,1)=P(1,2)
  P(3,3)=1.5_dp*cthe
  P(4,4)=P(3,3)

end subroutine praylv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spline scattering phase matrix. Version 2008-04-09.
subroutine pspliv(P,XP,YP,YP2,the,np,pfflg)

  real(kind=dp), dimension(4,4), intent(out) :: P
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  real(kind=dp), intent(in) :: the
  integer, intent(in) :: np, pfflg
  integer :: j1
  real(kind=dp), dimension(361) :: YPA, YPA2

  ! Reset:
  P(:,:)=0.0_dp

  ! Phase matrix:
  do j1 = 1, np+1
    YPA(j1) =YP(j1,1,1)
    YPA2(j1)=YP2(j1,1,1)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(1,1))
  if (pfflg == 1) return

  do j1 = 1, np+1
    YPA(j1) =YP(j1,1,2)
    YPA2(j1)=YP2(j1,1,2)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(1,2))

  do j1 = 1, np+1
    YPA(j1) =YP(j1,2,2)
    YPA2(j1)=YP2(j1,2,2)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(2,2))

  do j1 = 1, np+1
    YPA(j1) =YP(j1,3,3)
    YPA2(j1)=YP2(j1,3,3)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(3,3))

  do j1 = 1, np+1
    YPA(j1) =YP(j1,3,4)
    YPA2(j1)=YP2(j1,3,4)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(3,4))

  do j1 = 1, np+1
    YPA(j1) =YP(j1,4,4)
    YPA2(j1)=YP2(j1,4,4)
  end do
  call splint(XP,YPA,YPA2,np+1,the,P(4,4))

  P(2,1)= P(1,2)
  P(4,3)=-P(3,4)
    
end subroutine pspliv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Spline scattering phase matrix. Version 2008-04-09.
subroutine psplivi(P,CSRN,XP,YP,YP2,pnorm,dthe,nrn,np,ncm)

  real(kind=dp), dimension(0:360,4,4), intent(in) :: P
  real(kind=dp), dimension(0:), intent(out) :: CSRN
  real(kind=dp), dimension(361), intent(out) :: XP
  real(kind=dp), dimension(361,4,4), intent(out) :: YP,YP2
  real(kind=dp), intent(out) :: pnorm
  real(kind=dp), intent(in) :: dthe
  integer, intent(in) :: nrn, np, ncm
  integer :: j1, j2, j3, temp
  real(kind=dp), dimension(4,4) :: R
  real(kind=dp), dimension(512) :: XI,WXI
  real(kind=dp), dimension(361) :: YPA, YPA2


  ! Compute the spline representation of the scattering phase matrix:
  do j1 = 1, np+1
    XP(j1)=(j1-1)*dthe
  end do
  do j3 = 1, 4
    do j2 = 1, 4
      do j1 = 1, np+1
        YP(j1,j2,j3)=P(j1-1,j2,j3)
        YPA(j1)     =P(j1-1,j2,j3)
      end do
      call spline(XP,YPA,np+1,0.0_dp,0.0_dp,YPA2)
      do j1 = 1, np+1
        YP2(j1,j2,j3)=YPA2(j1)
      end do
    end do
  end do

  ! Normalize:
  call gauleg(-1.0_dp,1.0_dp,XI,WXI,ncm)
  pnorm=0.0_dp
  do j1 = 1, ncm
    temp = 1
    call pspliv(R,XP,YP,YP2,acos(XI(j1)),np,temp)
    pnorm=pnorm+WXI(j1)*R(1,1)
  end do
  do j1 = 1, np+1
    do j2 = 1, 4
      do j3 = 1, 4
        YP(j1,j2,j3) =2.0_dp*YP(j1,j2,j3)/pnorm
        YP2(j1,j2,j3)=2.0_dp*YP2(j1,j2,j3)/pnorm
      end do
    end do
  end do

  ! Map scattering angle cosines with random numbers:
  call cspspliv(CSRN,XP,YP,YP2,nrn,np,ncm)

end subroutine psplivi


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mueller matrix multiplication by Fresnel reflection matrix (F <- RF).
subroutine reflect(F,r11,r12,r33,r34)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), intent(in) :: r11, r12, r33, r34
  integer :: j1
  real(kind=dp) :: q

  do j1 = 1, 4
    q = r11*F(1,j1)+r12*F(2,j1)
    F(2,j1) = r12*F(1,j1)+r11*F(2,j1)
    F(1,j1) = q
    q = r33*F(3,j1)+r34*F(4,j1)
    F(4,j1) = -r34*F(3,j1)+r33*F(4,j1)
    F(3,j1) = q
  end do

end subroutine reflect


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Mueller matrix multiplication by Fresnel refraction matrix (G <- TF).
subroutine refract(F,G,t11,t12,t33,t34)

  real(kind=dp), dimension(4,4), intent(in) :: F
  real(kind=dp), dimension(4,4), intent(out) :: G
  real(kind=dp), intent(in) :: t11, t12, t33, t34
  integer :: j1

  do j1 = 1, 4
    G(1,j1) = t11*F(1,j1)+t12*F(2,j1)
    G(2,j1) = t12*F(1,j1)+t11*F(2,j1)
    G(3,j1) = t33*F(3,j1)+t34*F(4,j1)
    G(4,j1) = -t34*F(3,j1)+t33*F(4,j1)
  end do
    
end subroutine refract


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fresnel reflection coefficient for the parallel polarization.
function RL(m1,m2,ki,kt) result(res)

  complex(kind=dp), intent(in) :: m1, m2, ki, kt
  complex(kind=dp) :: res
  complex(kind=dp) :: q1, q2

  q1 = m1**2
  q2 = m2**2
  res = (q2*ki-q1*kt)/(q2*ki+q1*kt)
     
end function RL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fresnel reflection coefficient for the perpendicular polarization.
function RR(ki,kt) result(res)

  complex(kind=dp), intent(in) :: ki, kt
  complex (kind=dp) :: res
  
  res = (ki-kt)/(ki+kt)

end function RR


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Monte Carlo radiative transfer.
! Version: 2018 September 13
subroutine rtmc(F,KE,KF,HL,HR,X,XP,YP,YP2,CSRN,qabs,len,abscf,omg,np,nrn)

  real(kind=dp), dimension(4,4), intent(inout) :: F
  real(kind=dp), dimension(3), intent(inout) :: KE, KF
  complex(kind=dp), dimension(3), intent(inout) :: HL
  complex(kind=dp), dimension(3), intent(out) :: HR
  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  real(kind=dp), dimension(0:), intent(in) :: CSRN
  real(kind=dp), intent(inout) :: qabs
  real(kind=dp), intent(in) :: len
  real(kind=dp), intent(in) :: abscf, omg
  integer, intent(in) :: np, nrn
  integer :: j1
  real(kind=dp) :: kekf, lenabs
  real(kind=dp), dimension(3) :: EL, ER

  do j1 = 1, 3
    X(j1) = X(j1)+len*KE(j1)
  end do

  call prosca(kekf,KE,KF)
  lenabs=len*abs(kekf)

  call absorb(F,qabs,lenabs,abscf)
  call absorbrt(F,qabs,omg)

  call ehk(EL,ER,HL,KE)
  call incrt(F,KE,EL,ER,CSRN,XP,YP,YP2,np,nrn)

  do j1 = 1, 3
    KF(j1)=KE(j1)
    HL(j1)=cmplx(EL(j1),0.0_dp,dp)
    HR(j1)=cmplx(ER(j1),0.0_dp,dp)
  end do

end subroutine rtmc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generation of a new propagation direction and Mueller matrix.
subroutine scart(I,K,ER,EL,CSRN,XP,YP,YP2,np,nrn,pflg)

  real(kind=dp), dimension(4,4), intent(inout) :: I
  real(kind=dp), dimension(3), intent(inout) :: K, ER, EL
  real(kind=dp), dimension(0:), intent(in) :: CSRN
  real(kind=dp), dimension(361), intent(in) :: XP
  real(kind=dp), dimension(361,4,4), intent(in) :: YP, YP2
  integer, intent(in) :: np, nrn, pflg
  integer :: j1, j2, mrn, temp
  real(kind=dp) :: ran2, nt1, c2psi, s2psi, rn, cthe, sthe, phi, nn
  real(kind=dp), parameter :: gtol=1.0e-6
  real(kind=dp), dimension(4,4) :: I1
  real(kind=dp), dimension(3) :: K1, ER1, EL1, T1, T2, N1
  real(kind=dp), dimension(4,4) :: P

  ! Temporary storage:
  do j1=1,4
    do j2=1,4
      I1(j1,j2) =I(j1,j2)
    end do
  end do

  do j1 = 1, 3
    K1(j1)=K(j1)
    EL1(j1)=EL(j1)
    ER1(j1)=ER(j1)
  end do

  ! Generate polar scattering angle, and compute scattering phase matrix:
  if(pflg == 1) then
    call random_number(ran2)
    rn=nrn*ran2
    mrn=int(rn)
    cthe=(real(mrn+1,dp)-rn)*CSRN(mrn)+(rn-real(mrn,dp))*CSRN(mrn+1)
    temp=0
    call pspliv(P,XP,YP,YP2,acos(cthe),np,temp)

  elseif(pflg == 2) then
    call random_number(ran2)
    rn=2.0_dp*(1.0_dp-2.0_dp*ran2)
    cthe=(sqrt(1.0_dp+rn**2)+rn)**(1.0_dp/3.0_dp)-(sqrt(1.0_dp+rn**2)-rn)**(1.0_dp/3.0_dp)
    call praylv(P,cthe)

  endif

  ! Generate azimuthal scattering angle:
  call random_number(ran2)
  phi=2.0_dp*pi*ran2

  sthe=sqrt(1.0_dp-cthe**2)

  do j1 = 1, 3
    K(j1)=sthe*cos(phi)*EL1(j1)+sthe*sin(phi)*ER1(j1)+cthe*K1(j1)
  end do

  ! Auxiliary coordinate system:
  do j1 = 1, 3
    N1(j1)=K(j1)-K1(j1)
    T2(j1)=K(j1)+K1(j1)
  end do

  call provec(T1,T2,N1)
  nt1=0.0_dp
  do j1 = 1, 3
    nt1=nt1+T1(j1)**2
  end do
  nt1=sqrt(nt1)
  do j1 = 1, 3
    T1(j1)=T1(j1)/nt1
  end do

  ! Rotation of the input Mueller matrix:
  c2psi=(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))**2-(EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))**2
  s2psi=-2.0d0*(EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))*(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))
  call frotl(I1,c2psi,s2psi)

  ! New direction vectors for the reflected ray:
  do j1=1, 3
    ER(j1)=T1(j1)
  end do
  call provec(EL,ER,K)
  call MMMM(I,P,I1)
  nn=I1(1,1)/I(1,1)
  do j1=1,4
    do j2=1,4
      I(j1,j2)=I(j1,j2)*nn
    end do
  end do

end subroutine scart


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Stores the scattered Mueller matrix into the scattering phase matrix.
subroutine scatter(S,FOUT,KOUT,ELOUT,ER,qsca,bin,bin0,nbin)
 
  real(kind=dp), dimension(361,4,4), intent(inout) :: S
  real(kind=dp), dimension(4,4), intent(inout) :: FOUT
  real(kind=dp) , dimension(3), intent(in) :: KOUT, ELOUT, ER
  real(kind=dp), intent(inout) :: qsca
  real(kind=dp), intent(in) :: bin, bin0
  integer, intent(in) :: nbin
  integer :: kbin
  real(kind=dp) :: c2psi, s2psi, psi, ran2

  call random_number(ran2)

  ! Scattering efficiency:
  qsca = qsca + FOUT(1,1)

  ! First and last bins, i.e., forward and backward directions,
  ! accounting for random orientation:
  if(KOUT(3) > bin0) then ! forward scattering

    c2psi = ER(2)**2-ELOUT(2)**2
    s2psi = 2.0_dp*ELOUT(2)*ER(2)
    call frotr(FOUT,c2psi,-s2psi)

    call random_number(ran2)
    psi = 2.0_dp*pi*ran2
    c2psi = cos(psi)**2-sin(psi)**2
    s2psi = 2.0_dp*sin(psi)*cos(psi)
    call frotl(FOUT,c2psi, s2psi)
    call frotr(FOUT,c2psi,-s2psi)

    kbin = 1
    
    call raysca(S,FOUT,kbin)

  else if(KOUT(3) < -bin0) then ! backscattering

    c2psi = ER(2)**2-ELOUT(2)**2
    s2psi = -2.0_dp*ELOUT(2)*ER(2)
    call frotr(FOUT,c2psi,-s2psi)
    call random_number(ran2)        
    psi = 2.0_dp*pi*ran2
    c2psi = cos(psi)**2-sin(psi)**2
    s2psi = 2.0_dp*sin(psi)*cos(psi)
    call frotl(FOUT,c2psi,-s2psi)
    call frotr(FOUT,c2psi,-s2psi)

    kbin = nbin+1        
    call raysca(S,FOUT,kbin)

  ! Other directions:
  else
            
    c2psi = ((-ER(1)*KOUT(2)+ER(2)*KOUT(1))**2 - &
      (ELOUT(1)*KOUT(2)-ELOUT(2)*KOUT(1))**2) / &
      (1.0-KOUT(3)**2)
    s2psi = (2.0_dp*(ELOUT(1)*KOUT(2)-ELOUT(2)*KOUT(1))* &
      (-ER(1)*KOUT(2)+ER(2)*KOUT(1))) / &
      (1.0_dp-KOUT(3)**2)

    call frotl(FOUT,c2psi,s2psi)

    c2psi = (KOUT(1)**2-KOUT(2)**2)/(1.0_dp-KOUT(3)**2)
    s2psi = 2.0_dp*KOUT(1)*KOUT(2)/(1.0_dp-KOUT(3)**2)
    call frotr(FOUT,c2psi,-s2psi)
    kbin = int(0.5_dp+acos(KOUT(3))/bin)+1        
    call raysca(S,FOUT,kbin)

  endif

end subroutine scatter   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For given directions of constant phase (KEI) and amplitude of an
! incident inhomogeneous plane wave (KFI), SNEL computes the corresponding
! directions for the reflected (KE1, KF1) and refracted inhomogeneous plane
! waves (KE2, KF2).
subroutine snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T2, &
                fki,fkt,nke,nkf,m1,m2,sthe)

  real(kind=dp), dimension(3), intent(in) :: KEI, KFI 
  real(kind=dp), dimension(2), intent(in) :: MAI
  real(kind=dp), dimension(3), intent(inout) :: KE1, KF1
  real(kind=dp), dimension(2), intent(out) :: MA1
  real(kind=dp), dimension(3), intent(inout) :: KE2, KF2
  real(kind=dp), dimension(2), intent(out) :: MA2
  real(kind=dp), dimension(3), intent(in) :: N, T2
  complex(kind=dp), intent(out) :: fki, fkt
  real(kind=dp), intent(in) :: nke, nkf
  complex(kind=dp), intent(in) :: m1, m2
  real(kind=dp), intent(out) :: sthe
  real(kind=dp) :: n1, k1, n2, k2, ns, ks, d1, d2, &
    nke2, nkf2, tkf, cthe, cpsi, spsi, cphi, q1, q2
  real(kind=dp), dimension(3) :: T3

  ! Refractive indices:
  n1 = real(m1,dp)
  k1 = aimag(m1)
  d1 = n1**2-k1**2

  n2 = real(m2,dp)
  k2 = aimag(m2)
  d2 = n2**2-k2**2

  ! Apparent refractive index in medium 1; auxiliary
  ! products of refractive index and sines:
   MA1 = MAI

  if(abs(nke) < 1.0_dp) then
    ns = MAI(1)*sqrt(1.0_dp-nke**2)
  else
    ns = 0.0_dp
  end if

  if(abs(nkf) < 1.0_dp) then
    ks = MAI(2)*sqrt(1.0_dp-nkf**2)
  else
    ks = 0.0_dp
  end if

  ! Tangent vector in the plane perpendicular
  ! to constant amplitude; azimuthal angle between
  ! the planes of constant amplitude and phase:
  if (nkf > -0.9999999_dp) then
    tkf = sqrt(abs(1.0_dp-nkf**2))
    T3 = (KFI - nkf*N)/tkf
    call prosca(cphi,T2,T3)
  else
    T3 = 0.0_dp
    cphi = 1.0_dp
  end if

  ! Apparent refractive index in medium 2:
  q1 = (ns**2+ks**2+d2)**2 - &
    4.0_dp*((ns*ks)**2+d2*ns**2-(n2*k2-ns*ks*cphi)**2)
  q2 = ns**2+ks**2+d2

  if(q1 > 0.0_dp) then
    MA2(1) = sqrt(0.5_dp*(q2+sqrt(q1)))
  else
    MA2(1) = sqrt(0.5_dp*q2)
  end if

  if(MA2(1)**2 .gt. d2) then  
    MA2(2) = sqrt(MA2(1)**2-d2)
  else
    MA2(2) = 0.0_dp
  end if

  ! Unit direction vectors of for the reflected ray:
  KE1 = KEI - 2.0_dp*nke*N
  KF1 = KFI - 2.0_dp*nkf*N

  ! Unit direction vectors of for the refracted ray,
  ! total reflection:
  sthe = ns/MA2(1)
  if(sthe < 1.0_dp) then
    cthe = sqrt(1.0_dp-sthe**2)
    KE2 = sthe*T2 - cthe*N
  else 
    return
  end if

  if(MA2(2) /= 0.0_dp) then
    spsi = ks/MA2(2)
    if(spsi < 1.0_dp) then
      cpsi = sqrt(1.0_dp-spsi**2)
    else
      spsi = 1.0_dp
      cpsi = 0.0_dp
    end if
    KF2 = spsi*T3-cpsi*N
  else
    KF2 = KE2
  end if

  ! Preparation for Fresnel coefficients:
  call prosca(nke2,N,KE2)
  call prosca(nkf2,N,KF2)
  fki = cmplx(MA1(1)*abs(nke),MA1(2)*abs(nkf),dp)
  fkt = cmplx(MA2(1)*abs(nke2),MA2(2)*abs(nkf2),dp)

end subroutine snel


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! For given directions of constant phase (KEI) and amplitude of an
! incident inhomogeneous plane wave (KFI), SNEL computes the corresponding
! directions for the reflected (KE1, KF1) and refracted inhomogeneous plane
! waves (KE2, KF2).
! Version for 2-layer code
subroutine snel2l(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                     fki,fkt,nke,nkf,m1,m2,sthe,totref,preventTR_in)
  real(kind=dp), dimension(3), intent(in) :: KEI, KFI 
  real(kind=dp), dimension(2), intent(in) :: MAI
  real(kind=dp), dimension(3), intent(inout) :: KE1, KF1
  real(kind=dp), dimension(2), intent(out) :: MA1
  real(kind=dp), dimension(3), intent(inout) :: KE2, KF2
  real(kind=dp), dimension(2), intent(out) :: MA2
  real(kind=dp), dimension(3), intent(in) :: N, T1,T2
  logical, intent(out) :: totref
  complex(kind=dp), intent(out) :: fki, fkt
  real(kind=dp), intent(in) :: nke, nkf
  complex(kind=dp), intent(in) :: m1, m2
  real(kind=dp), intent(out) :: sthe
  real(kind=dp) :: n1, k1, n2, k2, ns, ks, d1, d2, &
    nke2, nkf2, tkf, cthe, cpsi, spsi, cphi, q1, q2, dd
  integer :: j1
  real(kind=dp) :: Na, Nb, norm1, norm2
  real(kind=dp), dimension(3) :: T3
  integer :: preventTR
  integer,optional :: preventTR_in
  if(present(preventTR_in))then
     preventTR=preventTR_in
  else
     preventTR=0
  endif

  totref = .false.
! Refractive indices:

  n1 = real(m1)
  k1 = aimag(m1)
  d1 = n1**2-k1**2

  n2 = real(m2)
  k2 = aimag(m2)
  d2 = n2**2-k2**2

! Apparent refractive index in medium 1; auxiliary
! products of refractive index and sines:

  MA1 = MAI

  if (abs(nke).lt.1.0_dp) then
    ns = MAI(1)*sqrt(1.0_dp-nke**2)
  else
    ns = 0.0_dp
  end if

  if (abs(nkf) .lt. 1.0_dp) then
    ks = MAI(2)*sqrt(1.0_dp-nkf**2)
  else
    ks = 0.0_dp
  end if

! Tangent vector in the plane perpendicular
! to constant amplitude; azimuthal angle between
! the planes of constant amplitude and phase:

    if (nkf .gt. -0.9999999_dp) then
        tkf = sqrt(1.0_dp-nkf**2)
        norm1 = 0.0_dp
        do j1 = 1, 3
            T3(j1) = (KFI(j1) - nkf*N(j1))/tkf
            norm1 = norm1 + T3(j1)**2
        end do
        norm1 = sqrt(norm1)
        do j1 = 1, 3
            T3(j1) = T3(j1)/norm1
        end do
        call prosca(cphi,T2,T3)
    else
        T3 = T2
        cphi = 1.0_dp
    end if

! Apparent refractive index in medium 2:


    
    if(abs(nke)>=1.0_dp) then
        Na=0.0_dp
        dd=0.0_dp
    else
        Na=MAI(1)**2*(1.0_dp-nke**2)
        dd = sqrt(1.0_dp-nke**2)
    endif
    
    if(abs(nkf)>=1.0_dp) then
        Nb=0.0_dp
    else
        Nb=MAI(2)**2*(1.0_dp-nkf**2)
    endif


    q1 = (ns**2+ks**2+d2)**2 - &
        4.0_dp*((ns*ks)**2+d2*ns**2-(n2*k2-ns*ks*cphi)**2)
    q2 = ns**2+ks**2+d2
   

    !tmp = sqrt(0.5_dp*(q2+sqrt(q1)))
    !write(6,*) tmp

    !tmp = sqrt(0.5_dp*(q2+sqrt(q1)))
    !if(tmp**2>d2) then
    !    tmp2 = sqrt(tmp**2-d2)
    !    tmp3 = sqrt(0.5_dp*(q2+sqrt(q1))-d2)
    !else
    !    tmp2=0.0_dp
    !    tmp3=0.0_dp
    !endif

    q1 = Na**2+Nb**2+d2**2-2.0_dp*Na*d2+2.0_dp*Nb*d2 &
         -2.0_dp*Na*Nb+4.0_dp*(n2*k2)**2- &
        8.0_dp*(ns*ks*cphi*n2*k2)+4.0_dp*Na*Nb*cphi**2
    q2 = Na+Nb+d2
    
    

    MA2(1) = sqrt(0.5_dp*(q2+sqrt(q1)))

    !write(6,*) tmp-MA2(1),tmp,MA2(1)


    if((MA1(2)<epsilon(k2) .and. k2<epsilon(k2)) .or. 0.5_dp*(q2+sqrt(q1))-d2<epsilon(q1)) then
        MA2(1) = n2
        MA2(2) = 0.0_dp
    else
        MA2(2) = sqrt(0.5_dp*(q2+sqrt(q1))-d2)
    endif


    !write(6,*) MA1
    !write(6,*) MA2,sqrt(0.5_dp*(q2-sqrt(q1))),0.5_dp*(q2+sqrt(q1))-d2

 ! Unit direction vectors of for the reflected ray:
    norm1 = 0.0_dp
    norm2 = 0.0_dp
    do j1 = 1, 3
        KE1(j1)=KEI(j1)-2.0_dp*nke*N(j1)
        KF1(j1)=KFI(j1)-2.0_dp*nkf*N(j1)
        norm1=norm1+KE1(j1)**2
        norm2=norm2+KF1(j1)**2
    end do
    norm1=sqrt(norm1)
    norm2=sqrt(norm2)

    do j1 = 1, 3
        KE1(j1)=KE1(j1)/norm1
        KF1(j1)=KF1(j1)/norm2
    end do

    !if (abs(k1) .lt. ktol) then
    !    KF1=KE1
    !end if


! Unit direction vectors of for the refracted ray,
! total reflection:

    sthe = ns/MA2(1)
    if (sthe >= 1.0_dp) then
        !totref=.true.
        return
    end if


    cthe = sqrt(1.0_dp-sthe**2)
    norm1 = 0.0_dp
    do j1 = 1, 3
        KE2(j1) = sthe*T2(j1) - cthe*N(j1)
        norm1 = norm1+KE2(j1)**2
    end do
    norm1 = sqrt(norm1)
    do j1 = 1, 3
        KE2(j1) = KE2(j1)/norm1
    end do
   

    if (abs(MA2(2)) > 0.0_dp) then
        spsi = ks/MA2(2)
        if (spsi .lt. 1.0_dp) then
          cpsi = sqrt(1.0_dp-spsi**2)
        else
          spsi = 1.0_dp
          cpsi = 0.0_dp
        end if
        norm1 = 0.0_dp
        do j1 = 1, 3
            KF2(j1) = spsi*T3(j1)-cpsi*N(j1)
            norm1 = norm1+KF2(j1)**2
        end do
        norm1 = sqrt(norm1)
        do j1 = 1, 3
            KF2(j1) = KF2(j1)/norm1
        end do
    else
        KF2 = KE2
    end if

    !write(6,*) "MA2",MA2,sqrt(0.5_dp*(q2-sqrt(q1)))
    !write(6,*) abs(n2-sqrt(0.5_dp*(q2-sqrt(q1)))),abs(n2-sqrt(0.5_dp*(q2+sqrt(q1))))
    !if(abs(n2-sqrt(0.5_dp*(q2-sqrt(q1))))<abs(n2-sqrt(0.5_dp*(q2+sqrt(q1))))) then
    !    totref =.true.
    !endif

    !write(6,*) dot_product(N,KE2)
    !write(6,*) nke

    !if(ns/real(m2,kind=dp)>1.0_dp) then
    if(real(m1,kind=dp)*dd/real(m2,kind=dp)>1.0_dp .and. preventTR==1) then
        !write(6,*) "totref"
        totref =.true.
    endif
    
    !if (abs(k2) .lt. ktol) then
    !    KF2 = KE2
    !end if

! Preparation for Fresnel coefficients:
    call prosca(nke2,N,KE2)
    call prosca(nkf2,N,KF2)
    fki = cmplx(MA1(1)*abs(nke),MA1(2)*abs(nkf),kind=dp)
    fkt = cmplx(MA2(1)*abs(nke2),MA2(2)*abs(nkf2),kind=dp)

end subroutine snel2l


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fresnel refraction coefficient for the parallel polarization.
function TL(m1,m2,ki,kt) result(res)

  complex(kind=dp), intent(in) :: m1, m2, ki, kt
  complex(kind=dp) :: res
  complex(kind=dp) :: q1, q2

  q1 = m1**2
  q2 = m2**2
  res = 2.0_dp*m1*m2*ki/(q2*ki+q1*kt)

end function TL


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fresnel refraction coefficient for the perpendicular polarization.
function TR(ki,kt) result(res)

  complex(kind=dp), intent(in) :: ki, kt
  complex(kind=dp) :: res

  res = 2.0d0*ki/(ki+kt)

end function TR


END MODULE SIRISRADTRANS
