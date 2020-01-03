MODULE SIRISUTILS

! Utility functions for SIRIS4 package
!
! v2019-12-31
!
! Karri Muinonen, Timo Väisänen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  
  public
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Cleaning a string from comments and empty spaces in beginning or end
subroutine strip_string(strin,strout,strlen)

  character(*), intent(in) :: strin
  character(*), intent(out) :: strout
  integer, intent(out) :: strlen
  character, dimension(*), parameter :: comchars = (/ '!', '#' /)
  integer, parameter :: cclen = size(comchars)
  integer :: osl, i, ci
  
  osl = len(strin)
  write(strout,'(A)') strin
  
  do i=1,cclen
    ci = index(strout, comchars(i))
    if(ci == 0) cycle
    write(strout,'(A)') strout(1:ci-1)
  end do
  
  write(strout,'(A)') trim(adjustl(strout))
  strlen = len_trim(strout)
    
end subroutine strip_string



END MODULE SIRISUTILS
