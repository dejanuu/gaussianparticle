
subroutine mipro(I2,M,I1)

! Product of a Mueller matrix and a matrix of Stokes parameters.
! Version 2008-04-09.
!
! Copyright (C) 2008 Karri Muinonen


implicit none


integer, parameter ::                              &
sp = kind(1.0),                                     &
dp = selected_real_kind(2*precision(1.0_sp)),       &
qp = selected_real_kind(2*precision(1.0_dp))
real (kind=dp), dimension(4,4) :: M
real(kind=dp), dimension(4,4) :: I1,I2

I2=0.0_dp

I2 = matmul(M,I1)


!do j1 = 1, 4
!    I2(j1)=0.0_dp
!    do j2 = 1, 4
!        I2(j1)=I2(j1)+M(j1,j2)*I1(j2)
!    end do
!end do
end subroutine mipro




