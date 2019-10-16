subroutine frotlc(F2,F1,H1,H2,G1,G2)

! Rotation of the Mueller matrix from the left (F2 =K*F1).

    implicit none
    
    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
  real (kind=dp), dimension(4,4) :: F1, F2, K
  real (kind=dp), dimension(2,4) :: A      
  complex (kind=dp) :: img
  complex (kind=dp), dimension(3) :: H1, H2, G1, G2
  complex (kind=dp), dimension(4) :: B
  complex (kind=dp), dimension(2,2) :: Q

    img = dcmplx(0.0_dp,1.0_dp)

    call proscac(Q(1,1),H1,G1)
    call proscac(Q(1,2),H1,G2)
    call proscac(Q(2,1),H2,G1)
    call proscac(Q(2,2),H2,G2)

    A(1,1) = 0.5_dp*(abs(Q(1,2))**2+abs(Q(1,1))**2)
    A(1,2) = 0.5_dp*(abs(Q(1,2))**2-abs(Q(1,1))**2)
    A(1,3) = dreal(Q(1,2)*dconjg(Q(1,1)))
    A(1,4) = dimag(Q(1,2)*dconjg(Q(1,1)))

    A(2,1) = 0.5_dp*(abs(Q(2,2))**2+abs(Q(2,1))**2)
    A(2,2) = 0.5_dp*(abs(Q(2,2))**2-abs(Q(2,1))**2)
    A(2,3) = dreal(Q(2,2)*dconjg(Q(2,1)))
    A(2,4) = dimag(Q(2,2)*dconjg(Q(2,1)))

    B(1) = 0.5_dp*(Q(2,2)*dconjg(Q(1,2))+Q(2,1)*dconjg(Q(1,1)))
    B(2) = 0.5_dp*(Q(2,2)*dconjg(Q(1,2))-Q(2,1)*dconjg(Q(1,1)))
    B(3) = 0.5_dp*(Q(2,2)*dconjg(Q(1,1))+Q(2,1)*dconjg(Q(1,2)))
    B(4) = -0.5_dp*img*(Q(2,2)*dconjg(Q(1,1))-Q(2,1)*dconjg(Q(1,2)))

    K(1,1) = A(1,1)+A(2,1)
    K(1,2) = A(1,2)+A(2,2)
    K(1,3) = A(1,3)+A(2,3)
    K(1,4) = A(1,4)+A(2,4)

    K(2,1) = -A(1,1)+A(2,1)
    K(2,2) = -A(1,2)+A(2,2)
    K(2,3) = -A(1,3)+A(2,3)
    K(2,4) = -A(1,4)+A(2,4)

    K(3,1) = dreal(2.0_dp*B(1))
    K(3,2) = dreal(2.0_dp*B(2))
    K(3,3) = dreal(2.0_dp*B(3))
    K(3,4) = dreal(2.0_dp*B(4))
       
    K(4,1) = -dimag(2.0_dp*B(1))
    K(4,2) = -dimag(2.0_dp*B(2))
    K(4,3) = -dimag(2.0_dp*B(3))
    K(4,4) = -dimag(2.0_dp*B(4))

    call MMMM(F2,K,F1)

end subroutine frotlc

