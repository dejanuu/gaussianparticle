subroutine RAYINIC(F,KE,KF,HL,HR)

! Initial unit Mueller matrix, propagation direction, and parallel 
! and perpendicular coordinate axes.

  implicit none

    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
  
    real (kind=dp) :: KE(3), KF(3)
  complex (kind=dp) :: HL(3), HR(3)
  real (kind=dp), dimension(4,4) :: F

    KE = (/0.0_dp, 0.0_dp, 1.0_dp/)
    KF = (/0.0_dp, 0.0_dp, 1.0_dp/)
    HL = dcmplx((/1.0_dp, 0.0_dp, 0.0_dp/), 0.0_dp)
    HR = dcmplx((/0.0_dp, 1.0_dp, 0.0_dp/), 0.0_dp)    
  F  = reshape((/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp, &
                 0.0_dp, 0.0_dp, 1.0_dp, 0.0_dp, &
                 0.0_dp, 0.0_dp, 0.0_dp, 1.0_dp/), (/4, 4/))
    
end subroutine RAYINIC
