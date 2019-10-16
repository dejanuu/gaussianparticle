subroutine proveccn(XY,X,Y)

! Vector product for two complex 3-vectors.

    implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))    
    complex (kind=dp) :: q1
    complex (kind=dp), dimension(3) :: XY, X, Y

    call provecc(XY,X,Y)
    q1 = sqrt(XY(1)**2+XY(2)**2+XY(3)**2)
    XY = XY/q1

end subroutine proveccn
