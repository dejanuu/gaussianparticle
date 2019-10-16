subroutine proscac(XY,X,Y)

! Scalar product for two complex 3-vectors.

  implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))    
    complex (kind=dp) :: XY
    complex (kind=dp), dimension(3) :: X, Y

    XY = X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
           
end subroutine proscac
