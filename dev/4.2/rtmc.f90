subroutine rtmc(F,KE,KF,HL,HR,X,XP,YP,YP2,CSRN,qabs,len,abscf,omg,np,nrn)
 
! Monte Carlo radiative transfer.
!
! Author: Karri Muinonen, modified to f90 by Julia Martikainen
! Version: 2018 September 13

    use common
    implicit none
  
      
    integer  :: j1, j2, irn, np, nrn    
    real (kind=dp), dimension(4,4) :: F
    real (kind=dp), dimension(3) :: KE, KF, EL, ER, X
    real (kind=dp) :: qabs, len, abscf, omg, kekf, lenabs
    real (kind=dp), dimension(0:1000) :: CSRN
    real (kind=dp), dimension(361) :: XP
    real (kind=dp), dimension(361,4,4) :: YP, YP2
    complex (kind=dp), dimension(3) :: HL, HR


    do j1 = 1, 3
        X(j1) = X(j1)+len*KE(j1)
    end do


    call prosca(kekf,KE,KF)
    lenabs=len*abs(kekf)

    call absorb(F,qabs,lenabs,abscf)
    call absorbrt(F,qabs,omg)

    call ehk(EL,ER,HL,HR,KE)
    call incrt(F,KE,EL,ER,CSRN,XP,YP,YP2,np,nrn)

    do j1 = 1, 3
        KF(j1)=KE(j1)
        HL(j1)=cmplx(EL(j1),0.0_dp,kind=dp)
        HR(j1)=cmplx(ER(j1),0.0_dp,kind=dp)
    end do


end subroutine rtmc
