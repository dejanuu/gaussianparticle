subroutine scatter(S,FOUT,KOUT,ELOUT,ER,qsca,bin,bin0,nbin)
 
! Stores the scattered Mueller matrix into the scattering phase matrix.

  implicit none
       
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))   
!    integer (kind=sp) :: koko      
    integer  :: kbin, nbin
    real (kind=dp) :: qsca, bin, bin0, c2psi, s2psi, psi, pi, ran2
    real (kind=dp) , dimension(3) :: KOUT, ELOUT, ER
    real (kind=dp), dimension(4,4) :: FOUT
    real (kind=dp), dimension(361,4,4) :: S
    

    pi = 3.1415926535898_dp

  call random_number(ran2)

! Scattering efficiency:

    qsca = qsca + FOUT(1,1)

! First and last bins, i.e., forward and backward directions,
! accounting for random orientation:

    if (KOUT(3).gt.bin0) then ! forward scattering

        c2psi = ER(2)**2-ELOUT(2)**2
        s2psi = 2.0_dp*ELOUT(2)*ER(2)
        call frotr(FOUT,c2psi,-s2psi)

        call random_number(ran2)
        psi = 2.0_dp*pi*ran2
        c2psi = cos(psi)**2-sin(psi)**2
        s2psi = 2.0_dp*sin(psi)*cos(psi)
        call frotl(FOUT,c2psi, s2psi)
        call frotr(FOUT,c2psi,-s2psi)

        !kbin = 0
        kbin = 1
        
        call raysca(S,FOUT,kbin)

    else if (KOUT(3).lt.-bin0) then ! backscattering

        c2psi = ER(2)**2-ELOUT(2)**2
        s2psi = -2.0_dp*ELOUT(2)*ER(2)
        call frotr(FOUT,c2psi,-s2psi)
        call random_number(ran2)        
        psi = 2.0_dp*pi*ran2
        c2psi = cos(psi)**2-sin(psi)**2
        s2psi = 2.0_dp*sin(psi)*cos(psi)
        call frotl(FOUT,c2psi,-s2psi)
        call frotr(FOUT,c2psi,-s2psi)

        !kbin = nbin
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
        !kbin = int(0.5_dp+acos(KOUT(3))/bin)
        kbin = int(0.5_dp+acos(KOUT(3))/bin)+1        
        call raysca(S,FOUT,kbin)
    endif

end subroutine scatter       

