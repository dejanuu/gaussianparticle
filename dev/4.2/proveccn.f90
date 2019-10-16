subroutine proveccn(XY,X,Y)

! Vector product for two complex 3-vectors.

    use common
    implicit none
     
    complex (kind=dp) :: q1
    complex (kind=dp), dimension(3) :: XY, X, Y

    call provecc(XY,X,Y)
    q1 = sqrt(XY(1)**2+XY(2)**2+XY(3)**2)
    XY = XY/q1

end subroutine proveccn
