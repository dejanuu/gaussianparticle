subroutine partray0(X,mu,nu,cphi,sphi)

! Rotation from the particle to the ray coordinate system.

  use common
  implicit none    
       
  real (kind=dp):: mu, nu, cphi, sphi, q1, q2, q3 
  real (kind=dp), dimension(3) :: X

  q1 = X(1)
  q2 = X(2)
  q3 = X(3)
  X(1) = q1*mu*cphi+q2*mu*sphi-q3*nu
  X(2) = -q1*sphi+q2*cphi
  X(3) = q1*nu*cphi+q2*nu*sphi+q3*mu
       
end subroutine partray0
