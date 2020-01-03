MODULE SIRISCONSTANTS

! Notes.
!
! v2019-10-08
!
! Karri Muinonen, Timo Väisänen, Hannakaisa Lindqvist, Julia Martikainen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use, intrinsic :: iso_fortran_env
  
  public
  
  integer, parameter :: dp = REAL64
  integer, parameter :: qp = REAL128
  
  real(kind=dp), parameter :: pi = 3.1415926535897932_dp
  
  integer, parameter :: file_name_length = 256, line_str_length = 1024

END MODULE SIRISCONSTANTS
