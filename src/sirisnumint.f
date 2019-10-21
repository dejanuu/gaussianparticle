MODULE SIRISNUMINT

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  
  public
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gauleg(x1,x2,x,w,n)

  real(kind=dp), intent(in) :: x1, x2
  integer, intent(in) :: n
  real(kind=dp), dimension(n), intent(out) :: x, w
  integer :: i, j, m
  real(kind=dp) :: p1, p2, p3, pp, xl, xm, z, z1
  real(kind=dp), parameter :: EPS=3.0e-14

  m=(n+1)/2
  xm=0.5_dp*(x2+x1)
  xl=0.5_dp*(x2-x1)
  do i = 1, m
    z=cos(PI*(i-0.25_dp)/(n+0.5_dp))
1   continue
    p1=1.0_dp
    p2=0.0_dp
    do j = 1, n
      p3=p2
      p2=p1
      p1=((2.0_dp*j-1.0_dp)*z*p2-(j-1.0_dp)*p3)/j
    end do
    pp=n*(z*p1-p2)/(z*z-1.0_dp)
    z1=z
    z=z1-p1/pp
    if(abs(z-z1) > EPS) goto 1
    x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
    w(i)=2.0_dp*xl/((1.0_dp-z*z)*pp*pp)
    w(n+1-i)=w(i)
  end do

end subroutine gauleg


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine spline(x,y,n,yp1,ypn,y2)

  integer, intent(in) :: n
  real(kind=dp), dimension(n), intent(in) :: x, y
  real(kind=dp), intent(in) :: yp1, ypn
  real(kind=dp), dimension(n), intent(inout) :: y2
  integer :: i, k
  real(kind=dp) :: p, qn, sig, un
  integer, parameter :: NMAX=500
  real(kind=dp), dimension(NMAX) :: u

  if(yp1 > 0.99d30) then
    y2(1)=0.0_dp
    u(1)=0.0_dp
  else
    y2(1)=-0.5_dp
    u(1)=(3.0_dp/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do i = 2, n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.0_dp
    y2(i)=(sig-1.0_dp)/p
    u(i)=(6.0_dp*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  end do
  if(ypn > 0.99d30) then
    qn=0.0_dp
    un=0.0_dp
  else
    qn=0.50_dp
    un=(3.0_dp/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0_dp)
  do k = n-1, 1, -1
    y2(k)=y2(k)*y2(k+1)+u(k)
  end do

end subroutine spline


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine splint(xa,ya,y2a,n,x,y)

  integer, intent(in) :: n
  real(kind=dp), dimension(n), intent(in) :: xa, y2a, ya
  real(kind=dp), intent(in) :: x
  real(kind=dp), intent(out) :: y
  real(kind=dp) :: a, b, h
  integer :: k, khi, klo

  klo=1
  khi=n
  do while(khi-klo > 1)
    k=(khi+klo)/2
    if(xa(k) > x) then
      khi=k
    else
      klo=k
    endif
  end do
  h=xa(khi)-xa(klo)
  if(h == 0.0_dp) stop 'Trouble in SPLINT: bad xa input.'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.0_dp

end subroutine splint

END MODULE SIRISNUMINT
