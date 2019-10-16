subroutine frotr(F,c2psi,s2psi)

! Rotation of the Mueller matrix from the right (F <- FK).

  use common
  implicit none  
       
    integer :: j1
    real (kind=dp) :: q, s2psi, c2psi
    real (kind=dp), dimension(4,4) :: F

    do j1 = 1, 4
      q = c2psi*F(j1,2)-s2psi*F(j1,3)
        F(j1,3) = s2psi*F(j1,2)+c2psi*F(j1,3)
        F(j1,2) = q
  end do
  
end subroutine frotr
