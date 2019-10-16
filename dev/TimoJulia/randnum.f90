! Init random number generator, either from
! a given seed or from the system clock
! 
! Antti Penttilä, 2017

    
SUBROUTINE init_random(seed)
  INTEGER, OPTIONAL, INTENT(INOUT) :: seed
  INTEGER :: seeddim, clock
  INTEGER, DIMENSION(:), ALLOCATABLE :: ranseed
  LOGICAL :: clockseed
  
  IF(PRESENT(seed)) THEN
    IF(seed <= 0) THEN
      clockseed = .TRUE.
    ELSE
      clockseed = .FALSE.
    ENDIF
  ELSE
    clockseed = .TRUE.
  ENDIF
  
  CALL RANDOM_SEED(SIZE=seeddim)
  ALLOCATE(ranseed(seeddim))
  
  IF(clockseed) THEN
    CALL SYSTEM_CLOCK(COUNT=clock)
    ranseed = clock + 37 * (/ (i - 1, i = 1, seeddim) /)
  ELSE
    ranseed = seed + 37 * (/ (i - 1, i = 1, seeddim) /)
  ENDIF
  
  CALL RANDOM_SEED(PUT=ranseed)
  DEALLOCATE(ranseed)

END SUBROUTINE init_random