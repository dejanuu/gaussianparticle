subroutine istri(X,K,N,XN,NT,IT,nk,s3,ntri,nis,nie)

! Determines the ray intersection point on the discretized object
! from the given position and unit direction vectors. 

! X = ray coordinates, K = keout, 
! XN = node coords, NT = outer triangle normals, IT = indexes 
! nk = N dot K, nis = leikkaako (1) vai ei (0) 

  use common
  implicit none
  
  integer :: ntri, nis, j1, nie
  real (kind=dp), dimension(130000,3) :: XN ! ncoord
  integer , dimension(260000,3) :: IT ! nnod
  real (kind=dp), dimension(260000,3) :: NT ! trinorm
  real (kind=dp), dimension(3) :: X, K, N, X1, X2, X3, H, U3, S
  real (kind=dp), dimension(2) :: D  
  real (kind=dp) :: nk, s3, p1, p2, p3, q1, q2, q3, kx

! Initialize:

  nis = 0
  s3 = 1.0_dp*(10_dp**10)
  
! Check all triangles:

loop: do j1 = 1, ntri     
  
! Triangle orientation must allow interaction:
! N on ulkonormaali ja K sateen etenemissuunta

    N = nie*NT(j1,1:3)
    
    call prosca(nk,N,K) ! nk on N ja K vektorien pistetulo
  
    if (nk .ge. 0.0_dp) cycle loop ! jos nk ge 0.0 niin ollaan väärällä puolella. 
    
! Triangle located behind the ray path:

    X1 = XN(IT(j1,1), 1:3) - X(1:3)
    X2 = XN(IT(j1,2), 1:3) - X(1:3)
    X3 = XN(IT(j1,3), 1:3) - X(1:3)
  
    call prosca(p1, K, X1)
    call prosca(p2, K, X2)
    call prosca(p3, K, X3)
              
    if (p1 .lt. 0.0_dp .and. p2 .lt. 0.0_dp .and. p3 .lt. 0.0_dp) cycle loop                    
              
! Triangle located fully under/above the ray path: 

         call prosca(q1,X,X1)
         call prosca(q2,X,X2)
         call prosca(q3,X,X3)
         call prosca(kx,K,X)
         p1 = q1-kx*p1
         p2 = q2-kx*p2
         p3 = q3-kx*p3

         if (p1.lt.0.0_dp .and. p2.lt.0.0_dp .and. p3.lt.0.0_dp) cycle loop
         if (p1.gt.0.0_dp .and. p2.gt.0.0_dp .and. p3.gt.0.0_dp) cycle loop
    
! Triangle located fully aside the ray path: 

        call provec(H,K,X) 
        call prosca(p1,H,X1)
        call prosca(p2,H,X2)
        call prosca(p3,H,X3)

      if (p1.lt.0.0_dp .and. p2.lt.0.0_dp .and. p3.lt.0.0_dp) cycle loop
        if (p1.gt.0.0_dp .and. p2.gt.0.0_dp .and. p3.gt.0.0_dp) cycle loop        
        
! Triangle intersection:

        U3(1:3) = NT(j1,1:3)
        X1(1:3) = XN(IT(j1,1),1:3)
        X2(1:3) = XN(IT(j1,2),1:3)
        X3(1:3) = XN(IT(j1,3),1:3)

      call isplane(S,X,K,U3,D,X1,X2,X3)



! Intersection update:

         if (S(1) .ge. 0.0_dp .and. S(2) .ge. 0.0_dp .and. &
        D(1)*S(2)+D(2)*S(1).le.D(1)*D(2) .and. S(3).ge.0.0_dp &
          .and. S(3).le.s3) then
          nis = j1
         !write (*,*) 'nis = ', nis
         s3 = S(3)
        !write (*,*) 's3 = ', s3
      end if 
  end do  loop

! Compute intersection point and normal:

    if (nis.gt.0_dp) then
      X(1:3) = X(1:3)+s3*K(1:3)
        N(1:3) = nie*NT(nis,1:3)
        call prosca(nk,N,K)
    end if

end subroutine istri
