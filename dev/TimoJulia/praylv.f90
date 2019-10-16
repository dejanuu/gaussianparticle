subroutine praylv(P,cthe)

! Rayleigh scattering phase matrix.
!
! Copyright by Karri Muinonen. All rights are reserved.
! Version 3.1, 2003 September 12.
implicit none


integer, parameter ::                              &
sp = kind(1.0),                                     &
dp = selected_real_kind(2*precision(1.0_sp)),       &
qp = selected_real_kind(2*precision(1.0_dp))
integer  :: j1,j2
real (kind=dp) :: cthe
real (kind=dp), dimension(4,4) :: P



! Reset:

do j1 = 1, 4
    do j2 = 1, 4
        P(j1,j2)=0.0_dp
    end do
end do

! Phase matrix:

P(1,1)=0.75_dp*(1.0_dp+cthe**2)
P(2,2)=P(1,1)
P(1,2)=-0.75_dp*(1.0_dp-cthe**2)
P(2,1)=P(1,2)
P(3,3)=1.5_dp*cthe
P(4,4)=P(3,3)

end subroutine praylv


