subroutine incide(F1,KE1,KF1,HL1,HR1,MA1,F2,KE2,KF2,HL2,HR2,MA2, &
                        N,m1,m2,totref)

! For an incident Mueller matrix and ray coordinate system (subscript 1), 
! computes the reflected (subscript 1) and refracted (subscript 2) 
! Mueller matrices and ray coordinate systems.

    use common
    implicit none
    
    integer :: totref, j1
    real (kind=dp):: ti,ci,ct,st,t11,t12,t33,t34,r11,r12,r33,r34, &
      mre,mim,ctl,ctr,nke,nkf,kekf1,kekf2
    real(kind=dp) ::  kekf,ni,c2psi,s2psi,nke2,nkf2,tol,dm, &
        dmtol,norm1,norm2,ktol
    real (kind=dp), dimension(3) :: KE1, KF1, KE2, KF2, KEI, KFI, T1, T2, N
    real (kind=dp), dimension(2) :: MA1, MA2, MAI
    real (kind=dp), dimension(4,4) :: F1, F2, FI
    complex (kind=dp) :: RR,RL,TL,TR,m,frl,frr,ftl,ftr,m1,m2,fki,fkt,normc
    complex (kind=dp), dimension(3) :: HL1, HR1, HL2, HR2, HLI, HRI, HLI0, HRI0, &
      NC, KC1, KC2, KCI
    parameter (dmtol=1.0e-10_dp,ktol=1.0e-10_dp)
    logical :: totref0

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
        call frotlc(F1,FI,HLI,HRI,HLI0,HRI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                fki,fkt,nke,nkf,kekf,m1,m2,st,totref0)

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
        
        ctr=(1.0_dp-abs(frr)**2)/abs(ftr)**2
        ctl=(1.0_dp-abs(frl)**2)/abs(ftl)**2

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

!  Complex unit vectors of the refracted rays:

        if (abs(MA2(2)) .ge. ktol) then
            KC2 = cmplx(MA2(1)*KE2, MA2(2)*KF2,kind=dp)
        else
            KC2 = cmplx(MA2(1)*KE2, 0.0_dp,kind=dp)
        end if

        call proveccn(HR2,KC2,NC)
        call proveccn(HL2,HR2,KC2)

        if(totref0) then

            !!MUISTA KERÄTÄ ENERGIA TALTEEN JOSSAIN
            !!TAI ANNA TAKAISIN F1:lle
            totref=1
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
        call frotlc(F1,FI,HLI,HRI,HLI0,HRI0)

! Unit direction vectors for the reflected and refracted ray:

        call snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                 fki,fkt,nke,nkf,kekf,m1,m2,st,totref0)

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

end subroutine incide


! *******************************************************************

function RL(m1,m2,ki,kt) result(res)

    ! Fresnel reflection coefficient for the parallel polarization.
    
    use common
    implicit none
    
    complex (kind=dp) :: m1, m2, ki, kt, q1, q2, res
    
    q1 = m1**2
    q2 = m2**2
    res = (q2*ki-q1*kt)/(q2*ki+q1*kt)
    
    end function RL
    
    
    function RR(ki,kt) result(res)
    
    ! Fresnel reflection coefficient for the perpendicular polarization.
    
    use common
    implicit none
    
    complex (kind=dp) :: ki, kt, res    
    
    res = (ki-kt)/(ki+kt)
    
    end function RR
    
    
    function TL(m1,m2,ki,kt) result(res)
    
    ! Fresnel refraction coefficient for the parallel polarization.
    
    use common
    implicit none
    
    complex (kind=dp) :: m1, m2, ki, kt, q1, q2, res       
    
    q1 = m1**2
    q2 = m2**2
    res = 2.0_dp*m1*m2*ki/(q2*ki+q1*kt)
    
    end function TL
    
    
    function TR(ki,kt) result(res)
    
    ! Fresnel refraction coefficient for the perpendicular polarization.
    
    use common
    implicit none
    
    complex (kind=dp) :: ki, kt, res       
    
    res = 2.0_dp*ki/(ki+kt)
    
    end function TR
    
    
