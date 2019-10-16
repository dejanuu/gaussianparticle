#ifndef DEBUG
#define    DEBUGGING .false.
#else
#define    DEBUGGING .true.
#endif


module constants
    !!Defines precisions and other constants.
    !!Copyright (C) 2016 Timo Väisänen and University of Helsinki
    !!All rights reserved.
    !!The new BSD License is applied to this software, see LICENSE.txt
    use, intrinsic :: iso_fortran_env
    implicit none
    integer, parameter :: dp = REAL64         
    integer, parameter :: qp = REAL128                      !!global precision(quad)
    real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)    !!pi
    logical, parameter :: debugging = DEBUGGING             !!debugging
    character(*), parameter :: str_line ="//////////////////////////////////////////////////////////////////"
    !!String line in output
    integer, parameter :: fnlen = 128

!    real(dp), parameter :: mu = 4.0*pi*10.0**(-7.0)         !!some constant

end module

