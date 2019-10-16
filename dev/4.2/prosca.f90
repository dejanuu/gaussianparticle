subroutine prosca(XY,X,Y)

! Scalar product for two 3-vectors.

  use common
  implicit none

  real (kind=dp):: XY
  real (kind=dp), dimension(3) :: X, Y

  XY=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)    
    
end subroutine prosca
