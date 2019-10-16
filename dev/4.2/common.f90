MODULE common

  PUBLIC
  
  INTEGER, PARAMETER ::               &
    sp = SELECTED_REAL_KIND(6, 37),   &
    dp = SELECTED_REAL_KIND(15, 307), &
    qp = SELECTED_REAL_KIND(33, 4931)

  REAL(KIND=dp), PARAMETER :: pi=3.1415926535897932_dp

END MODULE common
