subroutine raypart0c(X,mu,nu,cphi,sphi)

! Rotation from the ray to the particle coordinate system (complex).

  implicit none
  
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))  
  real (kind=dp) :: mu, nu, cphi, sphi
    complex (kind=dp) :: q1, q2, q3
    complex (kind=dp), dimension(3) :: X

    q1 = X(1)
    q2 = X(2)
    q3 = X(3)
    X(1) = q1*mu*cphi-q2*sphi+q3*nu*cphi
    X(2) = q1*mu*sphi+q2*cphi+q3*nu*sphi
    X(3) = -q1*nu+q3*mu
       
end subroutine raypart0c
