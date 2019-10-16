


subroutine mrotation(M,c2psi,s2psi)

! Mueller matrix for rotation. Version 2003-11-14.
!
! Copyright (C) 2003 Karri Muinonen

implicit none


integer, parameter ::                              &
sp = kind(1.0),                                     &
dp = selected_real_kind(2*precision(1.0_sp)),       &
qp = selected_real_kind(2*precision(1.0_dp))
integer  :: j1,j2
real (kind=dp) :: c2psi,s2psi
real (kind=dp), dimension(4,4) :: M



do j1 = 1, 4
    do j2 = 1, 4
        M(j1,j2)=0.0_dp
    end do
end do

M(1,1)=1.0_dp
M(2,2)=c2psi
M(2,3)=s2psi
M(3,3)=c2psi
M(3,2)=-s2psi
M(4,4)=1.0_dp

end subroutine mrotation




