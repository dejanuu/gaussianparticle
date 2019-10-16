subroutine mmmm(M3,M2,M1)
 
! Mueller matrix operation on another Mueller matrix.

  use common
  implicit none
 
      
    integer  :: j1, j2, j3    
    real (kind=dp), dimension(4,4) :: M1, M2, M3

    do j3 = 1, 4
      do j2 = 1, 4
          M3(j2,j3) = 0.0_dp
           do j1 = 1, 4
              M3(j2,j3) = M3(j2,j3)+M2(j2,j1)*M1(j1,j3)
            end do
    end do
  end do                

end subroutine mmmm
