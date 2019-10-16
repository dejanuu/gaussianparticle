! Radiative transfer

! Includes subroutines and functions:

! scart:        scattering process
! pmatrix1:     input scattering phase matrix
! psplivi:      initialization of spline interpolation for phase matrix
! pspliv:       spline interpolation for phase atrix
! cspspliv:     cosine-random-number map for spline phase function
! CSVRTBIS:     bisectioning algorithm for cosine-random-number mapping
! PCDFSPLIV:    c.d.f. of spline phase function
! kepnm:        Kepler's equation solved by Newton's method



subroutine pmatrix1(P,np)

! Input of discretized scattering phase matrix.



    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: np,j1,j2,j3
    real(kind=dp) :: p11,p12,p22,p33,p34,p44,the
    real(kind=dp), dimension(0:360,4,4) :: P
    real(kind=dp) :: tmp
    character (len=128) :: infile

    ! Reset:

    do j1 = 0, np
        do j2 = 1, 4
            do j3 = 1, 4
                P(j1,j2,j3)=0.0_dp
            end do
        end do
    end do

    ! Scattering phase matrix:
    IF(IARGC() .EQ. 0 .OR. IARGC() .EQ. 1) THEN
        infile = 'pmatrix1.in'
    ELSE
        CALL GETARG(2,infile)
    ENDIF
    open(unit=1, file=infile, status='old')

    do j1 = 0, np
        read(1,*) the,p11,p12,tmp,tmp, &
                  &   tmp,p22,tmp,tmp, &
                  &   tmp,tmp,p33,p34, &
                  &   tmp,tmp,tmp,p44
        P(j1,1,1)=p11
        P(j1,1,2)=p12
        P(j1,2,1)=p12
        P(j1,2,2)=p22
        P(j1,3,3)=p33
        P(j1,3,4)=p34
        P(j1,4,3)=-p34
        P(j1,4,4)=p44
    end do
    close(1)
end subroutine pmatrix1




subroutine scart(I,K,ER,EL,CSRN,XP,YP,YP2,pmx,np,nrn,pflg)

! Generation of a new propagation direction and Mueller matrix.

implicit none
integer, parameter ::                              &
sp = kind(1.0),                                     &
dp = selected_real_kind(2*precision(1.0_sp)),       &
qp = selected_real_kind(2*precision(1.0_dp))
integer  :: j1,j2,np,nrn,mrn,pflg, temp
real(kind=dp) :: ran2,nt1,c2psi,s2psi,rn,cthe,sthe,phi,nn,pmx,gtol,pi
real(kind=dp), dimension(4,4) :: I,I1
real(kind=dp), dimension(3) :: K,K1,ER,ER1,EL,EL1,T1,T2,N1
real(kind=dp), dimension(4,4) :: P
real(kind=dp), dimension(0:1000) :: CSRN
real(kind=dp), dimension(361) :: XP
real(kind=dp), dimension(361,4,4) :: YP, YP2


parameter (gtol=1.0d-6)

! New scattering direction
!call random_seed (size=koko)
!    allocate(seed(koko))
!    seed = seed1
!    call random_seed (put=seed)

pi=4.0_dp*atan(1.0_dp)

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

if (pflg.eq.1) then
call random_number(ran2)
rn=nrn*ran2
mrn=int(rn)
cthe=(dble(mrn+1)-rn)*CSRN(mrn)+(rn-dble(mrn))*CSRN(mrn+1)
temp=0
call pspliv(P,XP,YP,YP2,acos(cthe),np,temp)

elseif (pflg.eq.2) then
call random_number(ran2)
rn=2.0_dp*(1.0_dp-2.0_dp*ran2)
cthe=(sqrt(1.0_dp+rn**2)+rn)**(1.0_dp/3.0_dp)-(sqrt(1.0_dp+rn**2)-rn)**(1.0_dp/3.0_dp)
call praylv(P,cthe)

!    else
!        call random_number(ran2)
!        if (ran2.lt.whg) then
!           if (abs(ghg1).gt.gtol) then
!                call random_number(ran2)
!                cthe=(1.0_dp+ghg1**2-(1.0_dp-ghg1**2)**2/1.0_dp+(1.0_dp-2.0_dp*ran2)*ghg1)**2)/(2.0_dp*ghg1)
!            else
!                cthe=1.0_dp-2.0_dp*ran2
!            endif
!        else
!            if (abs(ghg2).gt.gtol) then
!                call random_number(ran2)
!                cthe=(1.0_dp+ghg2**2-(1.0_dp-ghg2**2)**2/(1.0_dp+(1.0_dp-2.0_dp*ran2)*ghg2)**2)/(2.0_dp*ghg2)
!            else
!                call random_number(ran2)
!                cthe=1.0_dp-2.0_dp*ran2
!            endif
!        endif
!        call pdhgv(P,whg,ghg1,ghg2,pmx,cthe)
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
call PROVEC(EL,ER,K)
call MMMM(I,P,I1)
nn=I1(1,1)/I(1,1)
do j1=1,4
    do j2=1,4
        I(j1,j2)=I(j1,j2)*nn
    end do
end do


end subroutine scart














subroutine psplivi(P,CSRN,XP,YP,YP2,pnorm,dthe,nrn,np,ncm)

! Spline scattering phase matrix. Version 2008-04-09.
!
! Copyright (C) 2008 Karri Muinonen


    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: nrn,np,ncm,j1,j2,j3, temp
    real(kind=dp) :: dthe,pnorm
    real(kind=dp), dimension(4,4) :: R
    real(kind=dp), dimension(0:360,4,4) :: P
    real(kind=dp), dimension(0:1000) :: CSRN
    real(kind=dp), dimension(361) :: XP,YPA,YPA2
    real(kind=dp), dimension(361,4,4) :: YP,YP2
    real(kind=dp), dimension(512) :: XI,WXI


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



subroutine pspliv(P,XP,YP,YP2,the,np,pfflg)

! Spline scattering phase matrix. Version 2008-04-09.
!
! Copyright (C) 2008 Karri Muinonen

    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: np,pfflg,j1,j2
    real(kind=dp) :: the
    real(kind=dp), dimension(4,4) :: P

    real(kind=dp), dimension(361) :: XP,YPA,YPA2
    real(kind=dp), dimension(361,4,4) :: YP,YP2



    ! Reset:

    do j1 = 1, 4
        do j2 = 1, 4
            P(j1,j2)=0.0_dp
        end do
    end do

    ! Phase matrix:

    do j1 = 1, np+1
        YPA(j1) =YP(j1,1,1)
        YPA2(j1)=YP2(j1,1,1)
    end do
    call splint(XP,YPA,YPA2,np+1,the,P(1,1))
    if (pfflg.eq.1) return

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



subroutine cspspliv(CSRN,XP,YP,YP2,nrn,np,ncm)

! Maps the cosines of scattering angles to random numbers for
! the spline scattering phase function of the spline scattering
! phase matrix. Version 2008-04-10.
!
! Copyright (C) 2008 Karri Muinonen.


    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: nrn,np,ncm,j1
    real(kind=dp) :: CSVRTBIS,rn,cstol
    real(kind=dp), dimension(0:1000) :: CSRN
    real(kind=dp), dimension(361) :: XP
    real(kind=dp), dimension(361,4,4) :: YP,YP2

    parameter (cstol=1.0d-6)

    CSRN(0)=-1.0_dp
    do j1 = 1, nrn-1
        rn=dble(j1)/nrn
        CSRN(j1)=CSVRTBIS(XP,YP,YP2,rn,-1.0_dp,1.0_dp,cstol,np,ncm)
    end do
    CSRN(nrn)=1.0_dp
end subroutine cspspliv



double precision function CSVRTBIS(XP,YP,YP2,rn,x1,x2,xacc,np,ncm)


    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: JMAX,j,np,ncm
    real(kind=dp) ::PCDFSPLIV,x1,x2,xmid,dx,xacc,fmid,f,rn
    real(kind=dp), dimension(361) :: XP
    real(kind=dp), dimension(361,4,4) :: YP,YP2

    parameter (JMAX=40)

    fmid=PCDFSPLIV(x2,XP,YP,YP2,np,ncm)-rn
    f=PCDFSPLIV(x1,XP,YP,YP2,np,ncm)-rn

    if(f*fmid.ge.0.0_dp) stop 'Trouble in CSVRTBIS: root must be bracketed.'

        if(f.lt.0.0_dp)then
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
        if(fmid.le.0.0_dp)CSVRTBIS=xmid
        if(abs(dx).lt.xacc .or. fmid.eq.0.0_dp) return
    end do
    stop 'Trouble in CSVRTBIS: too many bisections.'
end function CSVRTBIS


double precision function PCDFSPLIV(xi,XP,YP,YP2,np,ncm)

! Computes the cumulative distribution function for the spline
! scattering phase function of the spline scattering phase matrix.
! Version 2008-04-10.
!
! Copyright (C) 2008 Karri Muinonen.

    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: np,ncm,j1, temp
    real(kind=dp) :: xi
    real(kind=dp), dimension(361) :: XP
    real(kind=dp), dimension(361,4,4) :: YP,YP2
    real(kind=dp), dimension(4,4) :: P
    real(kind=dp), dimension(512) :: X,WX



    call gauleg(-1.0_dp,xi,X,WX,ncm)
    PCDFSPLIV=0.0_dp
    do j1 = 1, ncm
        temp=1
        call PSPLIV(P,XP,YP,YP2,acos(X(j1)),np,temp)
        PCDFSPLIV=PCDFSPLIV+WX(j1)*0.5_dp*P(1,1)
    end do
end function PCDFSPLIV



subroutine kepnm(ea,e,ma)

! Kepler's equation solved by Newton's method. Version 2008-04-09.
!
! Copyright (C) 2008 Karri Muinonen

    implicit none
    integer, parameter ::                              &
    sp = kind(1.0),                                     &
    dp = selected_real_kind(2*precision(1.0_sp)),       &
    qp = selected_real_kind(2*precision(1.0_dp))
    real(kind=dp) :: ea,e,ma,f,f1,f2,f3,dea,pi,tol
    parameter(pi=3.1415926535898d0)
    parameter(tol=1.0d-12)

    ! Initialize:

    ea=ma+0.85_dp*e*sign(1.0_dp,sin(ma-int(ma/(2.0_dp*pi))*2.0_dp*pi))
    dea=1.0d12

    ! Iterate:

    10 if (abs(dea).lt.tol) return

    f3=e*cos(ea)
    f2=e*sin(ea)
    f1=1.0_dp-f3
    f=ea-f2-ma

    dea=-f/f1
    dea=-f/(f1+f2*dea/2.0_dp)
    dea=-f/(f1+f2*dea/2.0_dp+f3*dea**2/6.0_dp)

    ea=ea+dea
    goto 10
end subroutine kepnm








