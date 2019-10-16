subroutine provecc(XY,X,Y)

! Vector product for two complex 3-vectors.

  implicit none
       
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))           
  complex (kind=dp), dimension(3) :: XY, X, Y       

    XY(1) = X(2)*Y(3)-X(3)*Y(2)
    XY(2) = X(3)*Y(1)-X(1)*Y(3)
    XY(3) = X(1)*Y(2)-X(2)*Y(1)    

end subroutine provecc
