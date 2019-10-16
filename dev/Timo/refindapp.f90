subroutine REFINDAPP(MA,m,kekf)

! Apparent refractive index.
! MA = apparent ref. index
! m = complex ref. index

  implicit none
  
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
  
    real (kind=dp) :: MA(2),kekf,n,k,d,ktol
    complex (kind=dp) :: m
    
    ktol = 1.0_dp*10**(-8)

    n = dreal(m)
    k = dimag(m)
    d = n**2-k**2


    MA(1) = dsqrt(0.5_dp*(d+dsqrt(d**2+4.0_dp*(n*k/kekf)**2)))
    if (abs(k).gt.ktol) then
      MA(2) = sqrt(0.5_dp*(-d+sqrt(d**2+4.0_dp*(n*k/kekf)**2)))
    else
      MA(2) = 0.0_dp
    end if
    
end subroutine REFINDAPP
