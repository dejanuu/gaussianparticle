subroutine frotl(F,c2psi,s2psi)

! Rotation of the Mueller matrix from the left (F <- KF).

  use common
  implicit none  

  integer :: j1
  real (kind=dp):: c2psi, s2psi, q
  real (kind=dp), dimension(4,4) :: F

    do j1 = 1, 4
      q = c2psi*F(2,j1)+s2psi*F(3,j1)
        F(3,j1) = -s2psi*F(2,j1)+c2psi*F(3,j1)
        F(2,j1) = q
  end do

end subroutine frotl
