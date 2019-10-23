MODULE SIRISINTERFACE

  use sirisradtrans

  public
  
contains


subroutine scatter_surface(KE1, KF1, HL1, HR1, MA1, F1, preserve1, N,   m1re,m1im, &
                        &  m2re, m2im,  KE2, KF2, HL2, HR2, MA2, F2, preserve2)
    use iso_c_binding, only : c_double,c_int
    implicit none

    integer, parameter :: dp = c_double
    integer(c_int), intent(inout) :: preserve1, preserve2
    real(c_double), intent(inout) :: KE1(3),KF1(3),F1(4,4),HL1(6),HR1(6),MA1(2)
    real(c_double), intent(inout) :: KE2(3),KF2(3),F2(4,4),HL2(6),HR2(6),MA2(2)
    real(c_double), intent(inout) :: N(3)
    real(c_double), intent(inout) :: m2re,m2im,m1re,m1im
   
    complex(c_double) :: HLF1(3),HRF1(3),HLF2(3),HRF2(3)
    complex(c_double) :: m1,m2
    real(c_double) :: nk,coll
    integer :: totref

    call convertToSiris(HL1,HR1,HLF1,HRF1,m1re,m1im,m2re,m2im,m1,m2)

    N = N/sqrt(dot_product(N,N))
    nk = dot_product(N,KE1)
    if(nk>0.0_dp) N=-N
    preserve1=0
    preserve2=0
    
    call incide(F1,KE1,KF1,HLF1,HRF1,MA1, &
                F2,KE2,KF2,HLF2,HRF2,MA2, &
                N,m1,m2,totref)

    !Two rays
    if(totref==0) then
        call convertToTracer(HLF2,HRF2,HL2,HR2)
        call convertToTracer(HLF1,HRF1,HL1,HR1)
        preserve1=1
        preserve2=1
    !Everything went through
    else if(totref==2) then
        call convertToTracer(HLF2,HRF2,HL2,HR2)
        preserve2=1
    !total reflection
    else if(totref==1) then
        call convertToTracer(HLF1,HRF1,HL1,HR1)
        preserve1=1
    else if(totref==3) then
        call convertToTracer(HLF1,HRF1,HL1,HR1)
        preserve2=3
        preserve1=1
    endif
end subroutine


!Diffuse scattering
subroutine scatter_diffuse(F,K,KF,HL,HR,XP,CSRN,coeffs,np,nrn,rands)
    use iso_c_binding, only : c_double,c_int
    use mathroutines, only : crossProduct
    implicit none
    integer, parameter :: dp = c_double 
    integer(c_int), intent(inout) :: np,nrn 
    real(c_double), intent(inout) :: F(4,4)
    real(c_double), intent(inout) :: XP(np)
    real(c_double), intent(inout) :: CSRN(nrn+1),coeffs(4,6,np)
    real(c_double), intent(inout) :: HL(6),HR(6)
    real(c_double), intent(inout)  :: rands(2)
    real(c_double), intent(inout) :: K(3),KF(3)
    real(kind=dp) :: EL(3), ER(3)

    EL(:)=(/HL(1),HL(3),HL(5)/)    
    ER = crossProduct(K,EL)
    EL = crossProduct(ER,K)

    ER = ER/sqrt(dot_product(ER,ER))
    EL = EL/sqrt(dot_product(EL,EL))

    call incrt_mp(F,K,EL,ER,CSRN,XP,coeffs,np,nrn,rands)
    KF(:)=K(:)
    HL(:)=(/EL(1),0.0_dp,EL(2),0.0_dp,EL(3),0.0_dp/)
    HR(:)=(/ER(1),0.0_dp,ER(2),0.0_dp,ER(3),0.0_dp/)
    
end subroutine 


subroutine incrt_mp(I,K,EL,ER,CSRN,XP,coeffs,np,nrn,rands)
    use iso_c_binding, only : c_double,c_int
    use sirisconstants, only : pi
    use mathRoutines, only : crossProduct
    use splinetools, only : interpolate
    implicit none
    integer, parameter :: dp = c_double 
    integer(kind=c_int), intent(in) :: np,nrn
    real(kind=c_double), intent(inout) :: I(4,4),K(3),ER(3),EL(3)
    real(kind=c_double), intent(in) :: XP(np),coeffs(4,6,np),CSRN(nrn+1),rands(2)
        
    real(kind=c_double) :: I1(4,4),K1(3),EL1(3),ER1(3),nt,c2psi,s2psi,N(3)
    real(kind=c_double) :: rn,cthe,tmpP(6),P(4,4),phi,sthe,N1(3),T1(3),T2(3),nn,norm
    integer :: temp,mrn,j1

    ! Temporary storage:
    I1(:,:) = I(:,:)
    
    K1(:)=K(:)
    EL1(:)=EL(:)
    ER1(:)=ER(:)
    
    ! Generate polar scattering angle, and compute scattering phase matrix:
    rn=nrn*rands(1)
    mrn=int(rn)
    cthe=(mrn+1.0_dp-rn)*CSRN(mrn+1)+(rn-mrn)*CSRN(mrn+2)
    sthe=sqrt(1.0_dp-cthe**2)

    ! Generate azimuthal scattering angle:
    
    phi=2.0_dp*pi*rands(2)
    K(:)=sthe*cos(phi)*EL1(:)+sthe*sin(phi)*ER1(:)+cthe*K1(:)
    
    !write(6,*) sqrt(dot_product(EL1,EL1))
    !write(6,*) sqrt(dot_product(ER1,ER1))
    norm = sqrt(dot_product(K(:),K(:)))
    K(:) = K(:)/norm
    
    ! Auxiliary coordinate system:
    N(:)=K(:)-K1(:)
    T2(:)=K(:)+K1(:)
    
    T1(:) =  crossProduct(T2,N)

    nt=0.0_dp
    do j1 = 1, 3
        nt=nt+T1(j1)**2
    end do
    nt=sqrt(nt)
    T1(:)=T1(:)/nt
    
    ! Rotation of the input Mueller matrix:
    c2psi=(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))**2-(EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))**2
    s2psi=-2.0_dp*(EL1(1)*T1(1)+EL1(2)*T1(2)+EL1(3)*T1(3))*(ER1(1)*T1(1)+ER1(2)*T1(2)+ER1(3)*T1(3))
    call frotl(I1,c2psi,s2psi)
    
    ! New direction vectors for the reflected ray:    
    ER(:)=T1(:)

    EL = crossProduct(ER,K)

    call interpolate(tmpP,xp,coeffs,acos(cthe))
    P=0.0_dp
    P(1,1)=tmpP(1)
    P(1,2)=tmpP(2)
    P(2,1)=tmpP(2)
    P(2,2)=tmpP(3)
    P(3,3)=tmpP(4)
    P(3,4)=tmpP(5)
    P(4,3)=-tmpP(5)
    P(4,4)=tmpP(6)
        
    I = matmul(P,I1)
    nn=I1(1,1)/I(1,1)
    I(:,:)=I(:,:)*nn
end subroutine 
    
    
subroutine convertToSIRIS(HL,HR,HLF,HRF,m1re,m1im,m2re,m2im,m1,m2)
    use iso_c_binding, only : c_double, c_int
    implicit none
    integer, parameter :: dp = c_double
    real(kind=c_double), intent(in) :: HL(6),HR(6),m1re,m1im,m2re,m2im
    complex(kind=c_double), intent(out) :: HLF(3),HRF(3),m1,m2

    m1 = cmplx(m1re,m1im,kind=dp)
    m2 = cmplx(m2re,m2im,kind=dp)

    HLF(1)=cmplx(HL(1),HL(2),kind=dp)
    HLF(2)=cmplx(HL(3),HL(4),kind=dp)
    HLF(3)=cmplx(HL(5),HL(6),kind=dp)

    HRF(1)=cmplx(HR(1),HR(2),kind=dp)
    HRF(2)=cmplx(HR(3),HR(4),kind=dp)
    HRF(3)=cmplx(HR(5),HR(6),kind=dp)
end subroutine


subroutine convertToTracer(HLF,HRF,HL,HR)
    use iso_c_binding, only : c_double,c_int
    implicit none
    integer, parameter :: dp = c_double
    real(kind=c_double), intent(out) :: HL(6),HR(6)
    complex(kind=dp), intent(in) :: HLF(3),HRF(3)

    HL(:)=(/real(HLF(1)),aimag(HLF(1)),real(HLF(2)),aimag(HLF(2)),real(HLF(3)),aimag(HLF(3))/)
    HR(:)=(/real(HRF(1)),aimag(HRF(1)),real(HRF(2)),aimag(HRF(2)),real(HRF(3)),aimag(HRF(3))/)
end subroutine


subroutine create_diffuse_scatterer(P,CSRN,XP,coeffs,np,nrn)
    use splinetools, only : initPresentation,bisectionMod, interpolate2
    use mathRoutines, only : gaussLQuads
    use iso_c_binding, only : c_double,c_int
    implicit none
    integer(kind=c_int), intent(inout) :: np,nrn
    real(kind=c_double), intent(inout) :: P(7,np)
    real(kind=c_double), intent(inout) :: CSRN(nrn+1)
    real(kind=c_double), intent(inout) :: XP(np)
    real(kind=c_double), intent(inout) :: coeffs(4,6,np)
    integer, parameter :: dp = c_double
    integer, parameter :: ncm = 32
    real(kind=c_double) :: XI(ncm),WXI(ncm)
    real(kind=c_double) :: pnorm
    real(kind=c_double) :: R,rn
    integer :: j1
    real(kind=dp), parameter :: tol = (10.0_dp**(-6.0_dp))

    call initPresentation(xp,coeffs,P,np,7)
    
    call gaussLQuads(-1.0_dp,1.0_dp,XI,WXI)
    pnorm = 0.0

    do j1=1,ncm
        call interpolate2(R,xp,coeffs,acos(XI(j1)))
        pnorm = pnorm+WXI(j1)*R
    end do

    do j1=1,np
        coeffs(:,:,j1)=2.0_dp*coeffs(:,:,j1)/pnorm
    enddo
    
    CSRN(1)=-1.0_dp
    do j1=2,nrn
        rn=(j1-1)/(nrn*1.0_dp)
        CSRN(j1)=bisectionMod(rn,-1.0_dp,1.0_dp,tol,ncm,xp,coeffs,np)
    enddo
    CSRN(nrn+1)=1.0_dp

end subroutine


subroutine finalize_ray(F,K,HL,theta,phi,nbins,rands)
    use iso_c_binding, only : c_double,c_int
    use mathRoutines, only : crossProduct
    implicit none
    integer, parameter :: dp = c_double

    real(c_double), intent(inout) :: F(4,4),K(3),HL(6),theta,phi,rands
    integer(c_int), intent(inout) :: nbins
    real(c_double) ::  mui,nui,cphi,sphi,EL(3),ER(3)

    mui = cos(theta)
    nui = sqrt(1.0_dp-mui**2)

    cphi = cos(phi)
    sphi = sin(phi)

    EL(:)=(/HL(1),HL(3),HL(5)/)    
    ER(:) = crossProduct(K,EL)
    EL(:) = crossProduct(ER,K)

    ER(:) = ER(:)/sqrt(dot_product(ER,ER))
    EL(:) = EL(:)/sqrt(dot_product(EL,EL))

    call partray0(K,mui,nui,cphi,sphi)
    call partray0(EL,mui,nui,cphi,sphi)
    call partray0(ER,mui,nui,cphi,sphi)
    call scatter_mod(F,K,EL,ER,nbins,rands)

end subroutine


subroutine scatter_mod(FOUT,KOUT,ELOUT,ER,nbin,rands)
    ! Stores the scattered Mueller matrix into the scattering phase matrix.
    use iso_c_binding, only : c_double,c_int
    use sirisconstants, only : pi
    implicit none
    integer, parameter :: dp=c_double
    integer, intent(in) :: nbin
    real (kind=dp), intent(in) :: KOUT(3), ELOUT(3), ER(3),rands
    real (kind=dp), intent(inout):: FOUT(4,4) 
    real (kind=dp) :: bin, bin0, c2psi, s2psi, psi
    
    bin = pi/nbin
    bin0 = cos(0.5_dp*bin)
    ! First and last bins, i.e., forward and backward directions,
    ! accounting for random orientation:

    if (KOUT(3).gt.bin0) then ! forward scattering
        c2psi = ER(2)**2-ELOUT(2)**2
        s2psi = 2.0_dp*ELOUT(2)*ER(2)
        call frotr(FOUT,c2psi,-s2psi)
        psi = 2.0_dp*pi*rands
        c2psi = cos(psi)**2-sin(psi)**2
        s2psi = 2.0_dp*sin(psi)*cos(psi)
        call frotl(FOUT,c2psi, s2psi)
        call frotr(FOUT,c2psi,-s2psi)
    else if (KOUT(3).lt.-bin0) then ! backscattering
        c2psi = ER(2)**2-ELOUT(2)**2
        s2psi = -2.0_dp*ELOUT(2)*ER(2)
        call frotr(FOUT,c2psi,-s2psi)
        psi = 2.0_dp*pi*rands
        c2psi = cos(psi)**2-sin(psi)**2
        s2psi = 2.0_dp*sin(psi)*cos(psi)
        call frotl(FOUT,c2psi,-s2psi)
        call frotr(FOUT,c2psi,-s2psi)
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
    endif

end subroutine       


END MODULE SIRISINTERFACE
