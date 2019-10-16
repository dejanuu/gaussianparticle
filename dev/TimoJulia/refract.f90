subroutine refract(F,G,t11,t12,t33,t34)

! Mueller matrix multiplication by Fresnel refraction matrix (G <- TF).

  implicit none
  
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
  integer  :: j1
  real (kind=dp) :: t11,t12,t33,t34
  real (kind=dp), dimension(4,4) :: F, G

    do j1 = 1, 4
      G(1,j1) = t11*F(1,j1)+t12*F(2,j1)
        G(2,j1) = t12*F(1,j1)+t11*F(2,j1)
        G(3,j1) = t33*F(3,j1)+t34*F(4,j1)
        G(4,j1) = -t34*F(3,j1)+t33*F(4,j1)
    end do
    
end subroutine refract

