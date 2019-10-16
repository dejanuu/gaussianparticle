




subroutine reflect(F,r11,r12,r33,r34)

! Mueller matrix multiplication by Fresnel reflection matrix (F <- RF).

  use common
  implicit none
  
  integer  :: j1
  real (kind=dp) :: r11,r12,r33,r34,q
  real (kind=dp), dimension(4,4) :: F

    do j1 = 1, 4
        q = r11*F(1,j1)+r12*F(2,j1)
        F(2,j1) = r12*F(1,j1)+r11*F(2,j1)
        F(1,j1) = q
        q = r33*F(3,j1)+r34*F(4,j1)
        F(4,j1) = -r34*F(3,j1)+r33*F(4,j1)
        F(3,j1) = q
  end do

end subroutine reflect
