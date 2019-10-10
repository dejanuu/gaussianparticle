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
