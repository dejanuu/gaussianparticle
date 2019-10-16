! Gsphere main program for running from command line with an input file

PROGRAM main
  USE constants
  USE ginterface
  IMPLICIT NONE

  REAL(KIND=dp) :: rmax, volume
  CHARACTER(len=fnlen) :: infile
  
  ! Input file specified in command line argument:

  CALL GETARG(1, infile)
  
  CALL init_ginterface(infile)

  CALL generate_gsphere(rmax)
  CALL compute_volume(volume)
  WRITE(*,*) "max radius",rmax,"volume",volume
  CALL save_gsphere()

END PROGRAM main
