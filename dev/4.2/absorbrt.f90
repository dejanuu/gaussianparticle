subroutine absorbrt(F,qabs,ot)

! Attenuates the Mueller matrix, and computes the absorption
! cross section.


use common
implicit none


integer :: j1,j2
real (kind=dp) :: qabs, ot
real (kind=dp), dimension(4,4) :: F


qabs=qabs+(1.0_dp-ot)*F(1,1)
do j1 = 1, 4
    do j2 = 1, 4
        F(j1,j2)=ot*F(j1,j2)
    end do
end do
end subroutine absorbrt


