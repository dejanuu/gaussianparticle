subroutine univec(U,d,X,Y)

! Computes the unit direction vector from two given positions.

  use common
  implicit none
    
  integer  :: j1
  real (kind=dp), dimension(3) :: U, X, Y
  real (kind=dp) :: d
       
  d = 0.0_dp
  do j1 = 1, 3 
    U(j1) = X(j1)-Y(j1)
    d = d+U(j1)**2
  end do
  if (d .eq. 0.0_dp) then
    write(*,*) 'Trouble in UNIVEC: zero vector.'
  end if
  d = sqrt(d)
  U(1:3) = U(1:3)/d
       
end subroutine univec
     