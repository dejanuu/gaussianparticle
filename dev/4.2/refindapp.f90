subroutine REFINDAPP(MA,m,kekf)

! Apparent refractive index.
! MA = apparent ref. index
! m = complex ref. index

    use common
    implicit none
  
    real (kind=dp) :: MA(2),kekf,n,k,d,ktol
    complex (kind=dp) :: m
    
    ktol = 1e-10_dp

    n = real(m)
    k = aimag(m)
    d = n**2-k**2


    MA(1) = sqrt(0.5_dp*(d+sqrt(d**2+4.0_dp*(n*k/kekf)**2)))
    if (abs(k).gt.ktol) then
      MA(2) = sqrt(0.5_dp*(-d+sqrt(d**2+4.0_dp*(n*k/kekf)**2)))
    else
      MA(2) = 0.0_dp
    end if
    
end subroutine REFINDAPP
