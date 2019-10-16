subroutine snel(KEI,KFI,MAI,KE1,KF1,MA1,KE2,KF2,MA2,N,T1,T2, &
                     fki,fkt,nke,nkf,kekf,m1,m2,sthe)

! For given directions of constant phase (KEI) and amplitude of an
! incident inhomogeneous plane wave (KFI), SNEL computes the corresponding
! directions for the reflected (KE1, KF1) and refracted inhomogeneous plane
! waves (KE2, KF2).

  implicit none
  
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
  real (kind=dp) :: n1,k1,n2,k2,ns,ks,d1,d2,nke,nkf, &
    nke2,nkf2,tkf,kekf,cthe,sthe,cpsi,spsi,cphi,q1,q2
  real (kind=dp), dimension(2) :: MA1, MA2, MAI 
  real (kind=dp), dimension(3) :: KE1, KF1, KE2, KF2, &
    KEI, KFI, N, T1, T2, T3 
  complex (kind=dp) :: m1,m2,fki,fkt

! Refractive indices:

    n1 = dreal(m1)
    k1 = dimag(m1)
    d1 = n1**2-k1**2

    n2 = dreal(m2)
    k2 = dimag(m2)
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

    if (q1 .gt. 0.0_dp) then
      MA2(1) = sqrt(0.5_dp*(q2+sqrt(q1)))
    else
      MA2(1) = sqrt(0.5_dp*q2)
    end if
    
    if (MA2(1)**2 .gt. d2) then  
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
    if (sthe .lt. 1.0_dp) then
    cthe = sqrt(1.0_dp-sthe**2)
        KE2 = sthe*T2 - cthe*N
    else 
        return
    end if

    if (MA2(2) .ne. 0.0_dp) then
      spsi = ks/MA2(2)
        if (spsi .lt. 1.0_dp) then
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
    fki = dcmplx(MA1(1)*abs(nke),MA1(2)*abs(nkf))
    fkt = dcmplx(MA2(1)*abs(nke2),MA2(2)*abs(nkf2))
       
end subroutine snel
