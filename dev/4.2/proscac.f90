subroutine proscac(XY,X,Y)

! Scalar product for two complex 3-vectors.

  use common
  implicit none
    
  complex (kind=dp) :: XY
  complex (kind=dp), dimension(3) :: X, Y

  XY = X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
           
end subroutine proscac
