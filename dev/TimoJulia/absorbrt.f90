subroutine absorbrt(F,qabs,ot)

! Attenuates the Mueller matrix, and computes the absorption
! cross section.


implicit none


integer, parameter ::                              &
sp = kind(1.0),                                     &
dp = selected_real_kind(2*precision(1.0_sp)),       &
qp = selected_real_kind(2*precision(1.0_dp))
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


