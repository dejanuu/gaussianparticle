MODULE SIRISMATERIAL

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
! Apparent refractive index.
! MA = apparent ref. index
! m = complex ref. index
subroutine refindapp(MA,m,kekf)

  real(kind=dp), dimension(2), intent(out) :: MA
  complex(kind=dp), intent(in) :: m
  real(kind=dp), intent(in) :: kekf
  real(kind=dp) :: n, k, d
  real(kind=dp), parameter :: ktol=1e-8

  n = real(m,dp)
  k = aimag(m)
  d = n**2-k**2

  MA(1) = sqrt(0.5_dp*(d+sqrt(d**2+4.0_dp*(n*k/kekf)**2)))
  if(abs(k) > ktol) then
    MA(2) = sqrt(0.5_dp*(-d+sqrt(d**2+4.0_dp*(n*k/kekf)**2)))
  else
    MA(2) = 0.0_dp
  end if
    
end subroutine refindapp


END MODULE SIRISMATERIAL
