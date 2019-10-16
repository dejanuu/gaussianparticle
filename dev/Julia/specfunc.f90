
subroutine LEGA(LEGPLM,x,lmax,m)

! The subroutine has been changed to agree with the sign convention
! in G. Arfken, Mathematical Methods for Physicists. It computes a sequence
! of associated Legendre functions with given m.

implicit none
integer, parameter ::                              &
sp = kind(1.0),                                  &
dp = selected_real_kind(2*precision(1.0_sp)),    &
qp = selected_real_kind(2*precision(1.0_dp))


integer  :: i,l,m,lmax
real (kind=dp) :: x,somx2,fact
real (kind=dp), dimension(0:256,0:256) :: LEGPLM

if(m.lt.0 .or. m.gt.lmax .or. abs(x).gt.1.0_dp) return

LEGPLM(m,m)=1.0_dp
if(m.gt.0) then
    somx2=sqrt((1.0_dp-x)*(1.0_dp+x))
    fact=1.0_dp
    do i=1,m
        LEGPLM(m,m)=LEGPLM(m,m)*fact*somx2
        fact=fact+2.0_dp
    end do
endif

if(lmax.eq.m) then
    LEGPLM(lmax,m)=LEGPLM(m,m)
else
    LEGPLM(m+1,m)=x*(2*m+1)*LEGPLM(m,m)
    if(lmax.eq.m+1) then
        LEGPLM(lmax,m)=LEGPLM(m+1,m)
    else
        do l=m+2,lmax
            LEGPLM(l,m)=(x*(2*l-1)*LEGPLM(l-1,m)-(l+m-1)*LEGPLM(l-2,m))/(l-m)
        end do
    endif
endif
END



DOUBLE PRECISION FUNCTION factrl(n)

implicit none
integer, parameter ::                              &
sp = kind(1.0),                                  &
dp = selected_real_kind(2*precision(1.0_sp)),    &
qp = selected_real_kind(2*precision(1.0_dp))

integer  :: n
integer  :: j,ntop
real (kind=dp) :: gammln
real (kind=dp), dimension(33) :: a
SAVE ntop,a
DATA ntop,a(1)/0,1.0_dp/
if (n.lt.0) then
    stop 'Trouble in FACTRL: negative factorial.'
else if (n.le.ntop) then
    factrl=a(n+1)
else if (n.le.32) then
    do j=ntop+1,n
        a(j+1)=j*a(j)
    end do
    ntop=n
    factrl=a(n+1)
else
    factrl=exp(gammln(n+1.0_dp))
endif
return
END



DOUBLE PRECISION FUNCTION gammln(xx)

implicit none
integer, parameter ::                              &
sp = kind(1.0),                                  &
dp = selected_real_kind(2*precision(1.0_sp)),    &
qp = selected_real_kind(2*precision(1.0_dp))

real (kind=dp) :: xx
integer  :: j
real (kind=dp) :: ser,stp,tmp,x,y
real (kind=dp), dimension(6) :: cof


cof(1)=76.18009172947146_dp
cof(2)=-86.50532032941677_dp
cof(3)=24.01409824083091_dp
cof(4)=-1.231739572450155_dp
cof(5)=0.001208650973866179_dp
cof(6)=-0.000005395239384953_dp
stp=2.5066282746310005_dp

x=xx
y=x
tmp=x+5.5_dp
tmp=(x+0.5_dp)*log(tmp)-tmp
ser=1.000000000190015_dp
do j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
end do
gammln=tmp+log(stp*ser/x)
return
END











