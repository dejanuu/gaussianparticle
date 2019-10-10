MODULE SIRISGAUSSIANSPHERE

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  use sirismath
  use sirisgeometry
  
  public
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
! Returns the Legendre coefficients for the correlation
! function with power-law Legendre coefficients.
subroutine cs1cf(CSCF,nu,lmin,lmax)

  real(kind=dp), dimension(0:), intent(out) :: CSCF
  real(kind=dp), intent(in) :: nu
  integer, intent(in) :: lmin, lmax
  integer :: l
  real(kind=dp) :: norm

  do l=0,lmin-1
    CSCF(l)=0.0_dp
  end do

  norm=0.0_dp
  do l=lmin,lmax
    CSCF(l)=1.0_dp/l**nu
    norm=norm+CSCF(l)
  end do

  do l=lmin,lmax
    CSCF(l)=CSCF(l)/norm
  end do

end subroutine cs1cf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the correlation length, and second and fourth derivatives
! for an autocorrelation function expressed in Legendre series with
! unnormalized coefficients. The coefficients are properly normalized
! in output.
subroutine csini(CL,crlen,cs2d,cs4d,lmax)

  real(kind=dp), dimension(0:), intent(inout) :: CL
  real(kind=dp), intent(out) :: crlen, cs2d, cs4d
  integer, intent(in) :: lmax
  integer :: l
  real (kind=dp) :: norm

  ! Normalization:
  norm=0.0_dp
  do l=0,lmax
    norm=norm+CL(l)
  end do
  do l=0,lmax
    CL(l)=CL(l)/norm
  end do

  ! Derivatives and correlation length:
  cs2d=0.0_dp
  cs4d=0.0_dp
  do l=1,lmax
    cs2d=cs2d-CL(l)*l*(l+1.0_dp)/2.0_dp
    cs4d=cs4d+CL(l)*l*(l+1.0_dp)*(3.0_dp*l**2+3.0_dp*l-2.0_dp)/8.0_dp
  end do
  if (cs2d == 0.0_dp) then
    crlen=2.0_dp
  else
    crlen=1.0_dp/sqrt(-cs2d)
  endif

end subroutine csini


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_idl(fbn,XT,IT,nnod,ntri)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri
  integer :: j1, j2, fu
  character(len=file_name_length) :: fn

  ! File name
  write(fn, '(A,A)') trim(fbn), ".idf"
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='write', status='replace')

  ! Write number of vertices and triangles
  write(fu,*) nnod, ntri
  ! Vertices
  do j1=1,nnod
    write(fu,*) (XT(j1,j2),j2=1,3)
  end do
  ! Triangle indices
  do j1=1,ntri
    write (fu,*) 3
    write (fu,*) (IT(j1,j2),j2=1,3)
  end do

  close(fu)

end subroutine save_idl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_matlab(fbn,XT,nnod)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, intent(in) :: nnod
  integer :: j1, fux, fuy, fuz
  character(len=file_name_length) :: fnx, fny, fnz

  ! File names
  write(fnx, '(A,A)') trim(fbn), "x.out"
  write(fny, '(A,A)') trim(fbn), "y.out"
  write(fnz, '(A,A)') trim(fbn), "z.out"
  
  ! File units and unit opening
  open(newunit=fux, file=trim(fnx), action='write', status='replace')
  open(newunit=fuy, file=trim(fny), action='write', status='replace')
  open(newunit=fuz, file=trim(fnz), action='write', status='replace')

  ! Write surface vertices
  do j1=1,nnod
    write (fux,*) XT(j1,1)
    write (fuy,*) XT(j1,2)
    write (fuz,*) XT(j1,3)
  end do

  close(fux)
  close(fuy)
  close(fuz)
   
end subroutine save_matlab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_vtk(fbn,XT,IT,nnod,ntri)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri
  integer :: j1, j2, fu
  character(len=file_name_length) :: fn

  ! File name
  write(fn, '(A,A)') trim(fbn), ".vtk"
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='write', status='replace')

  write(fu,'(A)') '# vtk DataFile Version 2.0'
  write(fu,'(A)') 'gsphere output            '
  write(fu,'(A)') 'ASCII                     '
  write(fu,'(A)') 'DATASET POLYDATA          '
  write(fu,'(A,I0,A)') 'POINTS ',nnod,' float'

  ! Write surface vertices
  do j1=1,nnod
    write (fu,*) (XT(j1,j2),j2=1,3)
  end do
  
  ! Triangle indices
  write (fu,'(A,I7,I7)') 'POLYGONS ',ntri,4*ntri
  do j1=1,ntri
    write (fu,*) 3,(IT(j1,j2)-1,j2=1,3)
  end do
  close(fu)

end subroutine save_vtk


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes associated Legendre functions from degree l=m
! up to l=lmax.
!
! Version 3.1, 2003 September 12.
subroutine plmg(PLM,x,lmax,m)

  real(kind=dp), dimension(0:,0:), intent(out) :: PLM
  real(kind=dp), intent(in) :: x
  integer, intent(in) :: lmax, m
  integer :: l, temp
  complex(kind=dp), dimension(0:lmax,0:lmax,-2:2) :: PLMM
  complex(kind=dp) :: i

  i=cmplx(0.0_dp,1.0_dp)

  ! Check degree, orders, and argument:
  if (lmax < 0) stop 'Trouble in PLMG: degree negative.'
  if (m > lmax .or. m < 0) stop 'Trouble in PLMG: order out of range.'
  if (abs(x) > 1.0_dp) stop 'Trouble in PLMG: argument out of range.'

  ! Compute associated Legendre functions with the help of
  ! the generalized spherical functions:
  temp=0
  call plmmg(PLMM,x,lmax,m,temp)

  do l=m,lmax
    PLM(l,m)=real(i**m*sqrt(FACTRL(l+m)/FACTRL(l-m))*PLMM(l,m,0), dp)
  end do

end subroutine plmg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes generalized spherical functions from degree max(abs(m1),abs(m2))
! up to lmax.
!
! Version 3.1, 2003 September 12.
subroutine plmmg(PLMM,x,lmax,m1,m2)

  complex(kind=dp), dimension(0:,0:,-2:), intent(out) :: PLMM
  real(kind=dp), intent(in) :: x
  integer, intent(in) :: lmax, m1, m2
  integer :: l, m0, m12, p12
  complex(kind=dp) :: i

  i=cmplx(0.0_dp,1.0_dp)

  ! Check degree, orders, and argument:
  if (lmax < 0) stop 'Trouble in PLMMG: degree negative.'
  if (abs(m1) > lmax .or. abs(m2) > min(2,lmax) .or. m1 < 0) stop 'Trouble in PLMMG: order out of range.'
  if (abs(x) > 1.0_dp) stop 'Trouble in PLMMG: argument out of range.'

  ! Compute generalized spherical functions:
  m0=max(abs(m1),abs(m2))
  m12=abs(m1-m2)
  p12=abs(m1+m2)

  if (m0 > 0) then
    if (m12 /= 0 .and. p12 /= 0) then
      PLMM(m0,m1,m2)=(-i)**m12/2.0_dp**m0*sqrt(factrl(2*m0)/(factrl(m12)*factrl(p12))*(1.0_dp-x)**m12*(1.0_dp+x)**p12)
    elseif (m12 == 0) then
      PLMM(m0,m1,m2)=1.0_dp/2.0_dp**m0*sqrt(FACTRL(2*m0)/FACTRL(p12)*(1.0_dp+x)**p12)
    else
      PLMM(m0,m1,m2)=(-i)**m12/2.0_dp**m0*sqrt(FACTRL(2*m0)/FACTRL(m12)*(1.0_dp-x)**m12)
    endif

    if (m0 == lmax) return

    PLMM(m0+1,m1,m2)=(2*m0+1)*(m0*(m0+1)*x-m1*m2)*PLMM(m0,m1,m2)/&
      (m0*sqrt(real((m0+1)**2-m2**2,dp))*sqrt(real((m0+1)**2-m1**2,dp)))

    if (m0+1 == lmax) return

    do l=m0+1,lmax-1
      PLMM(l+1,m1,m2)=((2*l+1)*(l*(l+1)*x-m1*m2)*PLMM(l,m1,m2)&
        -(l+1)*sqrt(real((l**2-m1**2)*(l**2-m2**2),dp))*PLMM(l-1,m1,m2))/&
        (l*sqrt(real(((l+1)**2-m1**2)*((l+1)**2-m2**2),dp)))
    end do

  else

    PLMM(0,0,0)=1.0_dp
    if (lmax == 0) return
    PLMM(1,0,0)=x
    if (lmax == 1) return

    do l=m0+1,lmax-1
      PLMM(l+1,0,0)=((2*l+1)*x*PLMM(l,0,0)-l*PLMM(l-1,0,0))/(l+1)
    end do

  endif

end subroutine plmmg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Radial distance in a given direction for a sample G-sphere.
real(kind=dp) function rgs(ACF,BCF,mu,phi,beta,lmin,lmax)

  real(kind=dp), dimension(0:,0:), intent(in) :: ACF, BCF
  real(kind=dp), intent(in) :: mu, phi, beta
  integer, intent(in) :: lmin, lmax

  rgs=exp(sgs(ACF,BCF,mu,phi,lmin,lmax)-0.5_dp*beta**2)

end function rgs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Discrete triangle representation for a sample G-sphere.
subroutine rgstd(X,N,MU,PHI,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)

  real(kind=dp), dimension(:,:), intent(out) :: X
  real(kind=dp), dimension(:,:), intent(out) :: N
  real(kind=dp), dimension(:), intent(in) :: MU, PHI
  real(kind=dp), dimension(0:,0:), intent(in) :: ACF, BCF
  real(kind=dp), intent(out) :: rmax
  real(kind=dp), intent(in) :: beta
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri, lmax
  integer :: j1, j2
  real(kind=dp) :: r, nu
  real(kind=dp), dimension(3) :: X1, X2, X3

  ! Node coordinates:
  rmax=0.0_dp
  do j1=1,nnod
    r=rsph1(ACF,BCF,MU(j1),PHI(j1),beta,lmax)
    nu=sqrt(1.0_dp-MU(j1)**2)
    X(j1,1)=r*nu*cos(PHI(j1))
    X(j1,2)=r*nu*sin(PHI(j1))
    X(j1,3)=r*MU(j1)
    if (r >= rmax) rmax=r
  end do

  ! Outer unit triangle normals:
  do j1=1,ntri
    do j2=1,3
      r=X(IT(j1,1),j2)
      X1(j2)=X(IT(j1,2),j2)-r
      X2(j2)=X(IT(j1,3),j2)-r
    end do
    call provec(X3,X1,X2)
    r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
    do j2=1,3
      N(j1,j2)=X3(j2)/r
    end do
  end do

end subroutine rgstd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the radius of the lognormal particle (spherical coordinate
! system).
! Version 3.1, 2003 September 12.
real(kind=dp) function rsph1(ALM,BLM,mu,phi,beta,lmax)

  real(kind=dp), dimension(0:,0:), intent(in) :: ALM,BLM
  real(kind=dp), intent(in) :: mu, phi, beta
  integer, intent(in) :: lmax

  rsph1=exp(ssph1(ALM,BLM,mu,phi,lmax)-0.5_dp*beta**2)

end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Logarithmic radial distance in a given direction for a sample G-sphere.
reaL(kind=dp) function sgs(ACF,BCF,mu,phi,lmin,lmax)

  real(kind=dp), dimension(0:,0:), intent(in) :: ACF, BCF
  real(kind=dp), intent(in) :: mu,phi
  integer, intent(in) :: lmin, lmax
  integer :: l, m
  real(kind=dp), dimension(lmax) :: CPHI, SPHI
  real(kind=dp), dimension(0:lmax,0:lmax) :: LEGP

  if (lmax == 0) then
    sgs=ACF(0,0)
    return
  endif

  ! Precomputation of sines, cosines, and associated Legendre functions:
  call LEGA(LEGP,mu,lmax,0)
  do m=1,lmax
    call LEGA(LEGP,mu,lmax,m)
    CPHI(m)=cos(m*phi)
    SPHI(m)=sin(m*phi)
  end do
  LEGP(0,0)=1.0_dp

  ! Sum up:
  sgs=0.0_dp
  do l=lmin,lmax
    sgs=sgs+LEGP(l,0)*ACF(l,0)
  end do
  do m=1,lmax
    do l=max(m,lmin),lmax
      sgs=sgs+LEGP(l,m)*(ACF(l,m)*CPHI(m)+BCF(l,m)*SPHI(m))
    end do
  end do

end function sgs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generates the sample spherical harmonics coefficients for the
! logarithmic radial distance of the G-sphere.
! Version 2002-12-16
subroutine sgscf(ACF,BCF,SCFSTD,lmax)

  real(kind=dp), dimension(0:,0:), intent(out) :: ACF, BCF
  real(kind=dp), dimension(0:,0:), intent(in) :: SCFSTD
  integer, intent(in) :: lmax
  integer :: l, m
  real(kind=dp) :: rn

  do l=0,lmax
    call rang2(rn)
    ACF(l,0)=rn*SCFSTD(l,0)
    BCF(l,0)=0.0_dp
  end do

  do m=1,lmax
    do l=m,lmax
      call rang2(rn)
      ACF(l,m)=rn*SCFSTD(l,m)
      call rang2(rn)
      BCF(l,m)=rn*SCFSTD(l,m)
    end do
  end do

end subroutine sgscf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Generates the standard deviations for the spherical harmonics
! coefficients of the logarithmic radial distance.
! Version 2002-12-16
subroutine sgscfstd(SCFSTD,CSCF,beta,lmax)

  real(kind=dp), dimension(0:,0:), intent(out) :: SCFSTD
  real(kind=dp), dimension(0:), intent(in) :: CSCF
  real(kind=dp), intent(in) :: beta
  integer, intent(in) :: lmax
  integer :: l, m

  SCFSTD(0,0)=beta*sqrt(CSCF(0))

  do l=1,lmax
    SCFSTD(l,0)=beta*sqrt(CSCF(l))
    do m=1,l
      SCFSTD(l,m)=SCFSTD(l,0)*sqrt(2.0_dp*factrl(l-m)/factrl(l+m))
    end do
  end do

end subroutine sgscfstd


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the logradius in a given direction (spherical coordinate
! system).
!
! Version 3.1, 2003 September 12.
real(kind=dp) function ssph1(ALM,BLM,mu,phi,lmax)

  real(kind=dp),dimension(0:,0:), intent(in) :: ALM, BLM
  real(kind=dp), intent(in) :: mu, phi
  integer, intent(in) :: lmax
  real(kind=dp), dimension(lmax) :: SPHI, CPHI
  integer :: l, m, temp
  real(kind=dp), dimension(0:lmax,0:lmax) :: PLM


  ! Instant return if lmax=0:
  if (lmax == 0) then
    ssph1=ALM(0,0)
    return
  endif

  ! Precomputation of sines, cosines, and associated Legendre functions:
  temp=0
  call plmg(PLM,mu,lmax,temp)
  do m=1,lmax
    call plmg(PLM,mu,lmax,m)
    CPHI(m)=cos(m*phi)
    SPHI(m)=sin(m*phi)
  end do
  PLM(0,0)=1.0_dp

  ! Sum for the logradius:
  ssph1=0.0_dp
  do l=0,lmax
    SSPH1=SSPH1+PLM(l,0)*ALM(l,0)
  end do
  do m=1,lmax
    do l=m,lmax
      ssph1=ssph1+PLM(l,m)*(ALM(l,m)*CPHI(m)+BLM(l,m)*SPHI(m))
    end do
  end do
  
end function ssph1


END MODULE SIRISGAUSSIANSPHERE
