subroutine isplane(S,X,K,U3,D,X1,X2,X3)

! Computes the intersection point on the plane determined
! by the three given positions.

  use common
  implicit none    

  real (kind=dp), dimension(3) :: S, X, K, X1, X2, X3, U1, U2, U3, DX
  real (kind=dp), dimension(2) :: D
  real (kind=dp) :: uu, u, u1dx, u2dx, u3dx, u1k, u2k, u3k, v, w, tol

    tol = 1.0e-14_dp

    call univec(U1,D(1),X2,X1)
    call univec(U2,D(2),X3,X1)
       
    DX = X1 - X

    call prosca(u3dx,U3,DX)
    
    if (abs(u3dx) .lt. tol) then
        S(3) = 0.0_dp
    else
        call prosca(u3k,U3,K)
        if (abs(u3k) .lt. tol) then
           write(*,*) U3, K, u3k, u3dx
          write(*,*) 'Trouble in ISPLANE: direction parallel to plane.'
        end if
        S(3)=u3dx/u3k
    end if

    call prosca(uu,U1,U2)
    u = 1.0_dp - uu**2
    if (abs(u).lt.tol) then
        write(*,*) U1,U2,u
        write(*,*) 'Trouble in ISPLANE: plane ill-defined.'
    end if 
    call prosca(u1dx,U1,DX)
    call prosca(u2dx,U2,DX)
    call prosca(u1k,U1,K)
    call prosca(u2k,U2,K)
    v = -u1dx+S(3)*u1k
    w = -u2dx+S(3)*u2k
    S(1) = (v-uu*w)/u
    S(2) = (w-uu*v)/u

end subroutine isplane
