module parameters
use constants
! Technical parameters for gsphere
implicit none

integer,parameter :: maxpar=2000, &
                     nmax=6000, & !for some reason this value cannot be changed. investigate!
                     nedmax=500
                     
real(kind=dp),parameter ::  tiny=1e-8, epsilon0=0.00001_dp
                           
end module parameters
