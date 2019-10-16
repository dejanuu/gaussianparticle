subroutine ehk(EL,ER,HL,HR,KE)
      
! Real-valued basis vectors from complex-valued basis vectors.
!
! Author: Karri Muinonen
! Version: 2018 September 13

  use common
  implicit none

  integer :: j1
  real (kind=dp), dimension(3) :: EL,ER,KE
  complex (kind=dp), dimension(3) :: HL, HR

  do j1 = 1, 3
      EL(j1)=real(HL(j1))
  end do
  call provecn(ER,KE,EL)
  call provecn(EL,ER,KE)

    
end subroutine ehk
