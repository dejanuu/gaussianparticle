subroutine snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                     fki,fkt,nke,nkf,kekf,m1,m2,sthe,totref)

! For given directions of constant phase (KEI) and amplitude of an
! incident inhomogeneous plane wave (KFI), SNEL computes the corresponding
! directions for the reflected (KE1, KF1) and refracted inhomogeneous plane
! waves (KE2, KF2).

  use common
  implicit none
  
  integer :: j1
  real (kind=dp) :: n1,k1,n2,k2,ns,ks,d1,d2,nke,nkf, &
    nke2,nkf2,tkf,kekf,cthe,sthe,cpsi,spsi,cphi,q1,q2, &
    ni,ki,norm1,norm2,ktol,Na,Nb,tmp,tmp2,tmp3
  real (kind=dp), dimension(2) :: MA1, MA2, MAI
  real (kind=dp), dimension(3) :: KE1, KF1, KE2, KF2, &
    KEI, KFI, N, T1, T2, T3
  complex (kind=dp) :: m1,m2,fki,fkt
  logical :: totref
  parameter (ktol=1.0e-10_dp)


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
    else
        Na=MAI(1)**2*(1.0_dp-nke**2)
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


    if((MA1(2)<epsilon(k2) .and. k2<epsilon(k2)) .or. 0.5_dp*(q2+sqrt(q1))<epsilon(q1)) then
        MA2(1) = n2
        MA2(2) = 0.0_dp
    else
        MA2(2) = sqrt(0.5_dp*(q2+sqrt(q1))-d2)
    endif

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
    if(abs(n2-sqrt(0.5_dp*(q2-sqrt(q1))))<abs(n2-sqrt(0.5_dp*(q2+sqrt(q1))))) then

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
       
end subroutine snel
