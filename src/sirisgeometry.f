MODULE SIRISGEOMETRY

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  use sirismath
  
  public
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Volume of triangle mesh object
function compute_volume(XT,IT,ntri) result(vol)

  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: ntri
  real(kind=dp) :: vol
  real(kind=dp) :: O(3), B(3), C(3), D(3)
  integer :: j1

  O=(/0.0_dp,0.0_dp,0.0_dp/)
  vol=0.0_dp
  do j1=1,ntri
    B=XT(IT(j1,1),:)
    C=XT(IT(j1,2),:)
    D=XT(IT(j1,3),:)
    vol=vol+simplex3volume(B,D,C)
  end do

end function compute_volume


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the intersection point on the plane determined
! by the three given positions.
subroutine isplane(S,X,K,U3,D,X1,X2,X3)

  real(kind=dp), dimension(3), intent(out) :: S
  real(kind=dp), dimension(3), intent(in) :: X, K, U3
  real(kind=dp), dimension(2), intent(out) :: D
  real(kind=dp), dimension(3), intent(in) :: X1, X2, X3
  real(kind=dp) :: uu, u, u1dx, u2dx, u3dx, u1k, u2k, u3k, v, w
  real(kind=dp), parameter :: tol=1.0e-14
  real(kind=dp), dimension(3) :: U1, U2, DX

  call univec(U1,D(1),X2,X1)
  call univec(U2,D(2),X3,X1)
     
  DX = X1 - X

  call prosca(u3dx,U3,DX)

  if(abs(u3dx) < tol) then
    S(3) = 0.0_dp
  else
    call prosca(u3k,U3,K)
    if (abs(u3k) < tol) then
      write(output_unit,*) U3, K, u3k, u3dx
      write(output_unit,*) 'Trouble in ISPLANE: direction parallel to plane.'
    end if
    S(3)=u3dx/u3k
  end if

  call prosca(uu,U1,U2)
  u = 1.0_dp - uu**2
  if(abs(u) < tol) then
    write(output_unit,*) U1,U2,u
    write(output_unit,*) 'Trouble in ISPLANE: plane ill-defined.'
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Determines the ray intersection point on the discretized object
! from the given position and unit direction vectors. 
! X = ray coordinates, K = keout, 
! XN = node coords, NT = outer triangle normals, IT = indexes 
! nk = N dot K, nis = leikkaako (1) vai ei (0) 
subroutine istri(X,K,N,XN,NT,IT,nk,s3,ntri,nis,nie)

  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), dimension(3), intent(in) :: K
  real(kind=dp), dimension(3), intent(out) :: N
  real (kind=dp), dimension(:,:), intent(in) :: XN ! ncoord
  real (kind=dp), dimension(:,:), intent(in) :: NT ! trinorm
  integer, dimension(:,:), intent(in) :: IT ! nnod
  real(kind=dp), intent(out) :: nk, s3
  integer, intent(in) :: ntri
  integer, intent(out) :: nis
  real(kind=dp), intent(in) :: nie
  integer :: j1
  real (kind=dp), dimension(2) :: D
  real(kind=dp), dimension(3) :: X1, X2, X3, H, U3, S
  real (kind=dp) :: p1, p2, p3, q1, q2, q3, kx

  ! Initialize:
  nis = 0
  s3 = 1.0e10_dp

  ! Check all triangles:
  loop: do j1 = 1, ntri     

    ! Triangle orientation must allow interaction:
    ! N on ulkonormaali ja K sateen etenemissuunta
    N = nie*NT(j1,1:3)
    call prosca(nk,N,K) ! nk on N ja K vektorien pistetulo
    if (nk >= 0.0_dp) cycle loop ! jos nk ge 0.0 niin ollaan vŠŠrŠllŠ puolella. 

    ! Triangle located behind the ray path:
    X1 = XN(IT(j1,1), 1:3) - X(1:3)
    X2 = XN(IT(j1,2), 1:3) - X(1:3)
    X3 = XN(IT(j1,3), 1:3) - X(1:3)

    call prosca(p1, K, X1)
    call prosca(p2, K, X2)
    call prosca(p3, K, X3)
              
    if (p1 < 0.0_dp .and. p2 < 0.0_dp .and. p3 < 0.0_dp) cycle loop                    
              
    ! Triangle located fully under/above the ray path: 
    call prosca(q1,X,X1)
    call prosca(q2,X,X2)
    call prosca(q3,X,X3)
    call prosca(kx,K,X)
    p1 = q1-kx*p1
    p2 = q2-kx*p2
    p3 = q3-kx*p3

    if (p1 < 0.0_dp .and. p2 < 0.0_dp .and. p3 < 0.0_dp) cycle loop
    if (p1 > 0.0_dp .and. p2 > 0.0_dp .and. p3 > 0.0_dp) cycle loop
    
    ! Triangle located fully aside the ray path: 
    call provec(H,K,X) 
    call prosca(p1,H,X1)
    call prosca(p2,H,X2)
    call prosca(p3,H,X3)

    if (p1 < 0.0_dp .and. p2 < 0.0_dp .and. p3 < 0.0_dp) cycle loop
    if (p1 > 0.0_dp .and. p2 > 0.0_dp .and. p3 > 0.0_dp) cycle loop        

    ! Triangle intersection:
    U3(1:3) = NT(j1,1:3)
    X1(1:3) = XN(IT(j1,1),1:3)
    X2(1:3) = XN(IT(j1,2),1:3)
    X3(1:3) = XN(IT(j1,3),1:3)

    call isplane(S,X,K,U3,D,X1,X2,X3)

    ! Intersection update:
    if(S(1) >= 0.0_dp .and. S(2) >= 0.0_dp .and. &
      D(1)*S(2)+D(2)*S(1) <= D(1)*D(2) .and. S(3) >= 0.0_dp &
      .and. S(3) <= s3) then
      nis = j1
      s3 = S(3)
    end if 
  end do  loop

  ! Compute intersection point and normal:
  if(nis > 0) then
    X(1:3) = X(1:3)+s3*K(1:3)
    N(1:3) = nie*NT(nis,1:3)
    call prosca(nk,N,K)
  end if

end subroutine istri


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Random particle orientation using Euler angles.
! Version: 2015-02-05
subroutine pranor(XN,NT,XN0,NT0,nnod,ntri)

  real(kind=dp), dimension(:,:), intent(out) :: XN
  real(kind=dp), dimension(:,:), intent(out) :: NT
  real(kind=dp), dimension(:,:), intent(in) :: XN0
  real(kind=dp), dimension(:,:), intent(in) :: NT0
  integer, intent(in) :: nnod, ntri
  integer :: j1, j2
  real(kind=dp) :: gamma, alpha, ran2
  real(kind=dp), dimension(3) :: X, CEU, SEU

  ! Euler angles in the laboratory reference frame (K) for
  ! expressing a vector given in K in the particle reference
  ! frame (K'):
  call random_number(ran2)
  gamma=2.0_dp*pi*ran2
  alpha=2.0_dp*pi*ran2
  CEU(1)=cos(gamma)
  CEU(2)=1.0_dp-2.0_dp*ran2
  CEU(3)=cos(alpha)
  SEU(1)=sin(gamma)
  SEU(2)=sqrt(1.0_dp-CEU(2)**2)
  SEU(3)=sin(alpha)

  ! Nodes:
  do j1 = 1, nnod
    do j2 = 1, 3
      X(j2)=XN0(j1,j2)
    end do
    call vproteut(X,CEU,SEU)
    do j2 = 1, 3
      XN(j1,j2)=X(j2)
    end do
  end do

  ! Normals:
  do j1 = 1, ntri
    do j2 = 1, 3
      X(j2)=NT0(j1,j2)
    end do
    call vproteut(X,CEU,SEU)
    do j2 = 1, 3
      NT(j1,j2)=X(j2)
    end do
  end do

end subroutine pranor


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function simplex3volume(u,v,b) result(A)

  real(kind=dp), intent(in) :: u(3), v(3), b(3)
  real(kind=dp) :: A

  A=(u(2)*v(3)-u(3)*v(2))*b(1)
  A=A+(u(3)*v(1)-u(1)*v(3))*b(2)
  A=A+(u(1)*v(2)-u(2)*v(1))*b(3)
  A=1.0_dp/6.0_dp*abs(A)

end function simplex3volume


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! TRIDS discretizes the spherical surface into altogether ntri=8*ntr**2
! triangles. It stores the nnod=4*ntr**2+2 nodes and right-handed node
! addresses for each triangle. ntr is the number of triangle rows in an
! octant.
subroutine trids(MU,PHI,IT,nnod,ntri,ntr)

  real(kind=dp), dimension(:), pointer, intent(out) :: MU, PHI
  integer, dimension(:,:), pointer, intent(out) :: IT
  integer, intent(out) :: nnod, ntri
  integer, intent(in) :: ntr
  integer :: j0, j1, j2, j3
  integer, dimension(0:360,0:720) :: NJJ
  real(kind=dp) :: the, fi, ct, st, cf, sf
  real(kind=dp), dimension(:,:), allocatable :: U
  
  ! Allocate tables
  allocate(MU(4*ntr**2+2), PHI(4*ntr**2+2), IT(8*ntr**2,3), U(4*ntr**2+2,3))

  ! NODES:

  ! Upper hemisphere including equator:
  nnod=1
  U(nnod,1)=0.0_dp
  U(nnod,2)=0.0_dp
  U(nnod,3)=1.0d0
  MU(nnod)=1.0_dp
  PHI(nnod)=0.0_dp
  NJJ(0,0)=nnod
  do j1=1,ntr
    the=j1*pi/(2*ntr)
    ct=cos(the)
    st=sin(the)
    do j2=0,4*j1-1
      fi=j2*pi/(2*j1)
      cf=cos(fi)
      sf=sin(fi)
      nnod=nnod+1
      U(nnod,1)=st*cf
      U(nnod,2)=st*sf
      U(nnod,3)=ct
      MU(nnod)=ct
      PHI(nnod)=fi
      NJJ(j1,j2)=nnod
      if (j2 == 0) NJJ(j1,4*j1)=nnod
    end do
  end do

  ! Lower hemisphere excluding equator:
  do j1=ntr-1,1,-1
    the=(2*ntr-j1)*pi/(2*ntr)
    ct=cos(the)
    st=sin(the)
    do j2=0,4*j1-1
      fi=j2*pi/(2*j1)
      cf=cos(fi)
      sf=sin(fi)
      nnod=nnod+1
      U(nnod,1)=st*cf
      U(nnod,2)=st*sf
      U(nnod,3)=ct
      MU(nnod)=ct
      PHI(nnod)=fi
      NJJ(2*ntr-j1,j2)=nnod
      if (j2 == 0) NJJ(2*ntr-j1,4*j1)=nnod
    end do
  end do

  nnod=nnod+1
  U(nnod,1)=0.0_dp
  U(nnod,2)=0.0_dp
  U(nnod,3)=-1.0_dp
  MU(nnod)=-1.0_dp
  PHI(nnod)=0.0_dp
  NJJ(2*ntr,0)=nnod

  if (nnod /= 4*ntr**2+2) stop 'Trouble in TRIDS: number of nodes inconsistent.'

  ! TRIANGLES:

  ! Upper hemisphere:
  ntri=0
  do j1=1,ntr
    do j3=1,4
      j0=(j3-1)*j1
      ntri=ntri+1
      IT(ntri,1)=NJJ(j1-1,j0-(j3-1))
      IT(ntri,2)=NJJ(j1,  j0       )
      IT(ntri,3)=NJJ(j1,  j0+1     )
      do j2=j0+1,j0+j1-1
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1,  j2         )
        IT(ntri,2)=NJJ(j1-1,j2  -(j3-1))
        IT(ntri,3)=NJJ(j1-1,j2-1-(j3-1))
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1-1,j2-(j3-1)  )
        IT(ntri,2)=NJJ(j1,  j2         )
        IT(ntri,3)=NJJ(j1,  j2+1       )
      end do
    end do
  end do

  ! Lower hemisphere:
  do j1=ntr+1,2*ntr
    do j3=1,4
      j0=(j3-1)*(2*ntr-j1)
      ntri=ntri+1
      IT(ntri,1)=NJJ(j1,  j0         )
      IT(ntri,2)=NJJ(j1-1,j0+1+(j3-1))
      IT(ntri,3)=NJJ(j1-1,j0  +(j3-1))
      do j2=j0+1,j0+(2*ntr-j1)
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1,  j2         )
        IT(ntri,2)=NJJ(j1-1,j2+(j3-1)  )
        IT(ntri,3)=NJJ(j1,  j2-1       )
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1,  j2         )
        IT(ntri,2)=NJJ(j1-1,j2+1+(j3-1))
        IT(ntri,3)=NJJ(j1-1,j2  +(j3-1))
      end do
    end do
  end do

  if (ntri /= 8*ntr**2) stop 'Trouble in TRIDS: number of triangles inconsistent.'

end subroutine trids



END MODULE SIRISGEOMETRY
