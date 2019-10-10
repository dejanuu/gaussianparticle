MODULE SIRISMATH

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  
  public
  
  interface prosca
    module procedure proscar, proscac
  end interface
  
  interface provec
    module procedure provecr, provecc
  end interface
  
  interface provecn
    module procedure provecrn, proveccn
  end interface
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real(kind=dp) function factrl(n)

  integer, intent(in) :: n
  integer :: j
  integer, save :: ntop=0
  real(kind=dp), dimension(33), save :: a = 1.0_dp
  
  if (n < 0) then
    stop 'Trouble in FACTRL: negative factorial.'
  else if (n <= ntop) then
    factrl=a(n+1)
  else if (n <= 32) then
    do j=ntop+1,n
      a(j+1)=j*a(j)
    end do
    ntop=n
    factrl=a(n+1)
  else
    factrl=exp(gammln(n+1.0_dp))
  endif

end function factrl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
real(kind=dp) function gammln(x)

  real(kind=dp), intent(in) :: x
  integer :: j
  real(kind=dp) :: ser,tmp,y
  real(kind=dp), parameter :: stp = 2.5066282746310005_dp
  real(kind=dp), dimension(6), parameter :: cof = (/ 76.18009172947146_dp, &
    -86.50532032941677_dp, &
    24.01409824083091_dp, &
    -1.231739572450155_dp, &
    0.001208650973866179_dp, &
    -0.000005395239384953_dp /)

  ser = 1.000000000190015_dp
  y=x
  tmp=x+5.5_dp
  tmp=(x+0.5_dp)*log(tmp)-tmp
  do j=1,6
      y=y+1.d0
      ser=ser+cof(j)/y
  end do
  gammln=tmp+log(stp*ser/x)

end function gammln


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Init random number generator, either from
! a given seed or from the system clock
subroutine init_random(seed)

  integer, optional, intent(in) :: seed
  integer :: seeddim, clock, i
  integer, dimension(:), allocatable :: ranseed
  logical :: clockseed
  
  if(present(seed)) then
    if(seed <= 0) then
      clockseed = .true.
    else
      clockseed = .false.
    endif
  else
    clockseed = .true.
  endif

  call random_seed(size=seeddim)
  allocate(ranseed(seeddim))
  
  if(clockseed) then
    call system_clock(count=clock)
    ranseed = clock + 37 * (/ (i - 1, i = 1, seeddim) /)
  else
    ranseed = seed + 37 * (/ (i - 1, i = 1, seeddim) /)
  endif
  
  call random_seed(put=ranseed)
  deallocate(ranseed)

end subroutine init_random


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The subroutine has been changed to agree with the sign convention
! in G. Arfken, Mathematical Methods for Physicists. It computes a sequence
! of associated Legendre functions with given m.
subroutine LEGA(LEGPLM,x,lmax,m)

  real (kind=dp), dimension(0:256,0:256), intent(out) :: LEGPLM
  real (kind=dp), intent(in) :: x
  integer, intent(in) :: lmax, m
  integer :: i, l
  real(kind=dp) :: somx2, fact

  if(m < 0 .or. m > lmax .or. abs(x) > 1.0_dp) return

  LEGPLM(m,m)=1.0_dp
  if(m > 0) then
    somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
    fact=1.0_dp
    do i=1,m
      LEGPLM(m,m)=LEGPLM(m,m)*fact*somx2
      fact=fact+2.0_dp
    end do
  endif

  if(lmax == m) then
    LEGPLM(lmax,m)=LEGPLM(m,m)
  else
    LEGPLM(m+1,m)=x*(2*m+1)*LEGPLM(m,m)
    if(lmax == m+1) then
      LEGPLM(lmax,m)=LEGPLM(m+1,m)
    else
      do l=m+2,lmax
        LEGPLM(l,m)=(x*(2*l-1)*LEGPLM(l-1,m)-(l+m-1)*LEGPLM(l-2,m))/(l-m)
      end do
    endif
  endif
end subroutine LEGA


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scalar product for two 3-vectors.
! For interface prosca.
subroutine proscar(XY,X,Y)

  real(kind=dp), intent(out):: XY
  real (kind=dp), dimension(3), intent(in) :: X, Y

  XY=X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
    
end subroutine proscar


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Scalar product for two complex 3-vectors.
! For interface prosca.
subroutine proscac(XY,X,Y)

  complex(kind=dp), intent(out) :: XY
  complex(kind=dp), dimension(3), intent(in) :: X, Y

  XY = X(1)*Y(1)+X(2)*Y(2)+X(3)*Y(3)
           
end subroutine proscac


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector product for two 3-vectors.
! For interface provec
subroutine provecr(XY,X,Y)

  real(kind=dp), dimension(3), intent(out) :: XY
  real(kind=dp), dimension(3), intent(in) :: X, Y

  XY(1)=X(2)*Y(3)-X(3)*Y(2)
  XY(2)=X(3)*Y(1)-X(1)*Y(3)
  XY(3)=X(1)*Y(2)-X(2)*Y(1)  
  
end subroutine provecr


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Normalized vector product for two 3-vectors.
! For interface provecn
subroutine provecrn(XY,X,Y)

  real(kind=dp), dimension(3), intent(out) :: XY
  real(kind=dp), dimension(3), intent(in) :: X, Y
  real(kind=dp) :: norm
  
  XY(1) = X(2)*Y(3)-X(3)*Y(2)
  XY(2) = X(3)*Y(1)-X(1)*Y(3)
  XY(3) = X(1)*Y(2)-X(2)*Y(1)    

  norm = sqrt(XY(1)**2+XY(2)**2+XY(3)**2)
  XY = XY/norm

end subroutine provecrn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector product for two complex 3-vectors.
! For interface provec
subroutine provecc(XY,X,Y)

  complex(kind=dp), dimension(3), intent(out) :: XY
  complex(kind=dp), dimension(3), intent(in) :: X, Y

  XY(1) = X(2)*Y(3)-X(3)*Y(2)
  XY(2) = X(3)*Y(1)-X(1)*Y(3)
  XY(3) = X(1)*Y(2)-X(2)*Y(1)    

end subroutine provecc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector product for two complex 3-vectors.
! For interface provecn
subroutine proveccn(XY,X,Y)

  complex(kind=dp), dimension(3), intent(out) :: XY
  complex(kind=dp), dimension(3), intent(in) :: X, Y
  complex(kind=dp) :: norm

  XY(1) = X(2)*Y(3)-X(3)*Y(2)
  XY(2) = X(3)*Y(1)-X(1)*Y(3)
  XY(3) = X(1)*Y(2)-X(2)*Y(1)  
  norm = sqrt(XY(1)**2+XY(2)**2+XY(3)**2)
  XY = XY/norm

end subroutine proveccn


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Returns a normally distributed random deviate with zero mean and
! unit variance.
subroutine rang2(r1)

  real(kind=dp), intent(out) :: r1
  logical, save :: flag = .false.
  real(kind=dp) :: ran2, q1, q2
  real(kind=dp), save :: r2

  if (flag) then
    r1=r2
    flag=.false.
    return
  endif

  do
    call random_number(ran2)
    r1=2.0_dp*ran2-1.0_dp
    call random_number(ran2)
    r2=2.0_dp*ran2-1.0_dp
    q1=r1**2+r2**2
    if (0.0_dp < q1 .and. q1 < 1.0_dp) exit
  end do

  q2=sqrt(-2.0_dp*log(q1)/q1)
  r1=r1*q2
  r2=r2*q2
  flag=.true.

end subroutine rang2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Computes the unit direction vector from two given positions.
subroutine univec(U,d,X,Y)

  real (kind=dp), dimension(3), intent(inout) :: U
  real (kind=dp), dimension(3), intent(in) :: X, Y
  real (kind=dp), intent(out) :: d
  integer  :: j1
       
  d = 0.0_dp
  do j1 = 1, 3 
    U(j1) = X(j1)-Y(j1)
    d = d+U(j1)**2
  end do
  if (d == 0.0_dp) then
    write(*,*) 'Trouble in UNIVEC: zero vector.'
  end if
  d = sqrt(d)
  U(1:3) = U(1:3)/d
       
end subroutine univec


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Transpose of vector rotation using Euler angles. Version 2003-11-07.
subroutine vproteut(X,CA,SA)

  real(kind=dp), dimension(3), intent(inout) ::  X
  real(kind=dp), dimension(3), intent(in) :: CA, SA
  
  call vrotz(X,CA(3),-SA(3))
  call vroty(X,CA(2),-SA(2))
  call vrotz(X,CA(1),-SA(1))

end subroutine vproteut


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector rotation about the y-axis. Version 2002-12-16.
subroutine vroty(X,c,s)

  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), intent(in) :: c,s
  real(kind=dp) :: q

  q = c*X(3)+s*X(1)
  X(1)=-s*X(3)+c*X(1)
  X(3)=q

end subroutine vroty


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vector rotation about the z-axis. Version 2002-12-16.
subroutine vrotz(X,c,s)

  real(kind=dp), dimension(3), intent(inout) :: X
  real(kind=dp), intent(in) :: c,s
  real(kind=dp) :: q

  q = c*X(1)+s*X(2)
  X(2)=-s*X(1)+c*X(2)
  X(1)=q

end subroutine vrotz


END MODULE SIRISMATH
