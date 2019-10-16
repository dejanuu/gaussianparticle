subroutine incide(F1,KE1,KF1,HL1,HR1,MA1,F2,KE2,KF2,HL2,HR2,MA2, &
                        N,m1,m2,totref)

! For an incident Mueller matrix and ray coordinate system (subscript 1), 
! computes the reflected (subscript 1) and refracted (subscript 2) 
! Mueller matrices and ray coordinate systems.

    implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
    integer :: totref, j1
  real (kind=dp):: ti,ci,ct,st,t11,t12,t33,t34,r11,r12,r33,r34, &
      mre,mim,ctl,ctr,nke,nkf,kekf1,kekf2
  real(kind=dp) ::  kekf
  real (kind=dp), dimension(3) :: KE1, KF1, KE2, KF2, KEI, KFI, T1, T2, N
  real (kind=dp), dimension(2) :: MA1, MA2, MAI  
  real (kind=dp), dimension(4,4) :: F1, F2, FI
  complex (kind=dp) :: RR,RL,TL,TR,m,frl,frr,ftl,ftr,m1,m2,fki,fkt
  complex (kind=dp), dimension(3) :: HL1, HR1, HL2, HR2, HLI, HRI, HLI0, HRI0, &
      NC, KC1, KC2, KCI

! Relative refractive index:

    m = m2/m1
    mre = dreal(m)
    mim = dimag(m)

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
      NC(j1) = dcmplx(N(j1),0.0_dp)
        KCI(j1) = dcmplx(MAI(1)*KEI(j1),MAI(2)*KFI(j1))
  end do
    call refindapp(MA1,m1,kekf)

! OBLIQUE INCIDENCE:

    if (nke .gt. -0.9999999_dp) then

! Two unit tangent vectors at the surface:

        ti = dsqrt(1.0_dp - nke**2)
        T2 = (KEI-nke*N)/ti
        call provec(T1,T2,N)

! Complex unit vectors in the plane of incidence and
! rotation of the input Mueller matrix:

        call proveccn(HRI,KCI,NC)
        call proveccn(HLI,HRI,KCI)
        call frotlc(F1,FI,HRI,HLI,HRI0,HLI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                fki,fkt,nke,nkf,kekf,m1,m2,st)

! Complex unit vectors of the reflected rays:

    KC1 = dcmplx(MA1(1)*KE1, MA1(2)*KF1)

        call proveccn(HR1,KC1,NC)
        call proveccn(HL1,HR1,KC1)

! Total Fresnel reflection, new reflected Mueller matrix:

        if (st .ge. 1.0_dp) then
          ci = -nke
           ct = sqrt(st**2-1.0_dp)
           frl = dcmplx(mre*ci,-ct)/dcmplx(mre*ci,ct)
           frr = dcmplx(ci,-mre*ct)/dcmplx(ci,mre*ct)
           r11 = 1.0_dp
           r12 = 0.0_dp
           r33 = dreal(frl*dconjg(frr))
           r34 = dimag(frl*dconjg(frr))
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
        ctl = abs(dreal(dconjg(m)*fkt/(m*fki)))   
        ctr = abs(dreal(fkt/fki))

        r11 = 0.5_dp*(cdabs(frl)**2+cdabs(frr)**2)
        r12 = 0.5_dp*(cdabs(frl)**2-cdabs(frr)**2)
        r33 = dreal(frl*dconjg(frr))
        r34 = dimag(frl*dconjg(frr))

        t11 = 1.0_dp-r11
        t12 = 0.5_dp*(ctl*cdabs(ftl)**2-ctr*cdabs(ftr)**2)
        t33 = dsqrt(ctl*ctr)*dreal(ftl*dconjg(ftr))
        t34 = dsqrt(ctl*ctr)*dimag(ftl*dconjg(ftr))

        call refract(F1,F2,t11,t12,t33,t34)
        call reflect(F1,r11,r12,r33,r34)

! Complex unit vectors of the refracted rays:

        KC2 = dcmplx(MA2(1)*KE2, MA2(2)*KF2)

        call proveccn(HR2,KC2,NC)
        call proveccn(HL2,HR2,KC2)
        totref = -1

! NORMAL INCIDENCE:

  else

! Two unit tangent vectors at the surface:

      T1 = dreal(HRI0)
        call provecn(T2,T1,KEI)
        call provec(T1,N,T2)

! Complex unit vectors in the plane of incidence and
! rotation of the input Mueller matrix:

        call proveccn(HRI,KCI,NC)   
        call proveccn(HLI,HRI,KCI)
        call frotlc(F1,FI,HRI,HLI,HRI0,HLI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                 fki,fkt,nke,nkf,kekf,m1,m2,st)

! Complex unit vectors of the reflected rays:

        call prosca(kekf1,KE1,KF1)
        if (abs(kekf1) .lt. 0.9999999_dp) then
      KC1 = dcmplx(MA1(1)*KE1, MA1(2)*KF1)
           call proveccn(HR1,KC1,NC)
           call proveccn(HL1,HR1,KC1)
        else
            HL1 = -HLI0
            HR1 = HRI0
        end if

! Complex unit vectors of the refracted rays:

        call prosca(kekf2,KE2,KF2)
        if (abs(kekf2) .lt. 0.9999999_dp) then
      KC2 = dcmplx(MA2(1)*KE2, MA2(2)*KF2)
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

! *******************************************************************

function RL(m1,m2,ki,kt) result(res)

! Fresnel reflection coefficient for the parallel polarization.

  implicit none

    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
    complex (kind=dp) :: m1, m2, ki, kt, q1, q2, res

    q1 = m1**2
  q2 = m2**2
  res = (q2*ki-q1*kt)/(q2*ki+q1*kt)
       
end function RL


function RR(ki,kt) result(res)

! Fresnel reflection coefficient for the perpendicular polarization.

    implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
    complex (kind=dp) :: ki, kt, res    

       res = (ki-kt)/(ki+kt)
       
end function RR


function TL(m1,m2,ki,kt) result(res)

! Fresnel refraction coefficient for the parallel polarization.

  implicit none

    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
    complex (kind=dp) :: m1, m2, ki, kt, q1, q2, res       

       q1 = m1**2
       q2 = m2**2
       res = 2.0_dp*m1*m2*ki/(q2*ki+q1*kt)
       
end function TL


function TR(ki,kt) result(res)

! Fresnel refraction coefficient for the perpendicular polarization.

  implicit none
       
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
    complex (kind=dp) :: ki, kt, res       

       res = 2.0d0*ki/(ki+kt)

end function TR
       
       
