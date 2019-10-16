subroutine absorb(F,qabs,len,abscf)
      
! Attenuates the Mueller matrix, and computes the absorption 
! cross section. 

  implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))    
    real (kind=dp) :: qabs,abscf,attcf,len
    real (kind=dp), dimension(4,4) :: F

    if (abs(abscf*len) .lt. 20.0_dp) then 
      attcf = exp(abscf*len)
    else
      attcf=0.0_dp
    end if

    qabs = qabs+(1.0_dp-attcf)*F(1,1)

    F = attcf*F
    
end subroutine absorb
