subroutine raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,NISBOX, &
           F,KE,KF,HL,HR,X,MA,pout,nis,jbox)

! Returns the Mueller matrix, propagation direction, parallel and 
! perpendicular coordinate axes, ray location, and the chord number.

  implicit none
  
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
       
  integer , dimension(50) :: PBOX, NISBOX 
  integer :: j1, j2, jbox, pout, nis
  real (kind=dp), dimension(3,50) :: KEBOX, KFBOX, XBOX
  real (kind=dp), dimension(2,50) :: MABOX
  real (kind=dp), dimension(4,4,50) :: FBOX
  real (kind=dp), dimension(4,4) :: F
  real (kind=dp), dimension(3) :: KE, KF, X 
  real (kind=dp), dimension(2) :: MA
  complex (kind=dp), dimension(3) :: HL, HR
  complex (kind=dp), dimension(3,50) :: HLBOX, HRBOX       

    KE = KEBOX(1:3,jbox)
    KF = KFBOX(1:3,jbox) ! onko KFBOX vai KEBOX (Karrin koodissa); testasin, ei nayta vaikuttavan tuloksiin...
    HL = HLBOX(1:3,jbox)
    HR = HRBOX(1:3,jbox)
    X = XBOX(1:3,jbox)    
    MA = MABOX(1:2,jbox)
    
    do j1 = 1, 4
      do j2 = 1, 4
          F(j1,j2) = FBOX(j1,j2,jbox)
        end do
  end do
       
    pout = PBOX(jbox)
    nis = NISBOX(jbox)
       
end subroutine raygetd1c  

