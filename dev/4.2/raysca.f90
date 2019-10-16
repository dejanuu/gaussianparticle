subroutine raysca(S,F,k)

! Stores the ray into the scattering phase matrix.

    use common
    implicit none

    integer  :: j1, j2, k
    real (kind=dp), dimension(4,4) :: F
    real (kind=dp), dimension(361,4,4) :: S   

    do j1 = 1, 4
    do j2 = 1,4
          S(k,j1,j2) = S(k,j1,j2)+F(j1,j2)
    end do
  end do

end subroutine raysca
