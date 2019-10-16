subroutine rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX, &
                  XBOX,MABOX,PBOX,NISBOX,F,KE,KF,HL,HR,X,MA, &
                  pout,nis,jbox)

! Stores the Mueller matrix, propagation direction, parallel and 
! perpendicular coordinate axes, ray location, chord number,
! and the triangle number.
 
  implicit none
       
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))      
    integer  :: j1,j2,jbox,pout,nis
    integer , dimension(50) :: PBOX, NISBOX   
  real (kind=dp), dimension(2) :: MA
  real (kind=dp), dimension(3) :: KE, KF, X
  real (kind=dp), dimension(4,4) :: F
  real (kind=dp), dimension(2,50) :: MABOX
  real (kind=dp), dimension(3,50) :: KEBOX, KFBOX, XBOX            
  real (kind=dp), dimension(4,4,50) :: FBOX
  complex (kind=dp), dimension(3) :: HL, HR
  complex (kind=dp), dimension(3,50) :: HLBOX, HRBOX  

    do j1 = 1, 3
      KEBOX(j1,jbox) = KE(j1)
        KFBOX(j1,jbox) = KF(j1)
        HLBOX(j1,jbox) = HL(j1)
        HRBOX(j1,jbox) = HR(j1)
        XBOX(j1,jbox) = X(j1)
  end do
    
    MABOX(1,jbox) = MA(1)
    MABOX(2,jbox) = MA(2)
    
    do j1 = 1, 4
      do j2 = 1, 4
          FBOX(j1,j2,jbox)=F(j1,j2)
    end do
  end do
       
    PBOX(jbox) = pout
    NISBOX(jbox) = nis
    
end subroutine rayputd1c
