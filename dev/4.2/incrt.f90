subroutine incrt(F,K,EL,ER,CSRN,XP,YP,YP2,np,nrn)

! Approximate Rayleigh scattering algorithm for the
! generation of a new propagation direction and Mueller
! matrix.


use common
implicit none


integer  :: j1,j2,np,nrn,mrn,irn
real(kind=dp) :: nt, c2psi, s2psi, rn, cthe, sthe, phi, nn, ran2, cphi, sphi
real(kind=dp), dimension(3) :: K, K1, EL, EL1, ER, ER1, T1, T2, N
real (kind=dp), dimension(4,4) :: F, F1, P
real (kind=dp), dimension(0:1000) :: CSRN
real (kind=dp), dimension(361) :: XP
real (kind=dp), dimension(361,4,4) :: YP, YP2
common irn



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
