MODULE rng

! Random number generator initialization
! Antti Penttilä, 2018

CONTAINS

SUBROUTINE init_random_seed(inseed)

  INTEGER, OPTIONAL, INTENT(IN) :: inseed
  INTEGER :: i, n, clock
  INTEGER :: seed1
  INTEGER, DIMENSION(:), ALLOCATABLE :: seed
  LOGICAL :: clockseed = .TRUE.

  CALL RANDOM_SEED(size = n)
  ALLOCATE(seed(n))
  
  IF(PRESENT(inseed)) THEN
    IF(inseed > 0) THEN
      seed = inseed
      clockseed = .FALSE.
    END IF
  END IF
  
  IF(clockseed) THEN
    CALL SYSTEM_CLOCK(COUNT=clock)
    seed1 = abs( mod((clock*181)*((getpid()-83)*359), 104729)) 
    seed = seed1 + 37 * (/ (i - 1, i = 1, n) /)
  END IF

  CALL RANDOM_SEED(PUT = seed)
  DEALLOCATE(seed)
  
END SUBROUTINE init_random_seed

END MODULE rng