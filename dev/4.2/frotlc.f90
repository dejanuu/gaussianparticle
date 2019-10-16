subroutine frotlc(F2,F1,H1,H2,G1,G2)

! Rotation of the Mueller matrix from the left (F2 =K*F1).

    use common
    implicit none
    
    integer :: j1,j2
    real (kind=dp) :: norm
    real (kind=dp), dimension(4,4) :: F1, F2, K
    complex (kind=dp) :: D
    complex (kind=dp), dimension(3) :: H1, H2, G1, G2
    complex (kind=dp), dimension(2,2) :: Q
    complex (kind=dp), dimension(2) :: U, V


    call PROSCAC(Q(1,1),H1,G1)
    call PROSCAC(Q(1,2),H1,G2)
    call PROSCAC(Q(2,1),H2,G1)
    call PROSCAC(Q(2,2),H2,G2)

    D=Q(1,1)*Q(2,2)-Q(1,2)*Q(2,1)
    U(1)=Q(2,2)/D
    U(2)=-Q(1,2)/D
    V(1)=-Q(2,1)/D
    V(2)=Q(1,1)/D

    K(1,1)=.5_dp*(abs(U(1))**2+abs(V(1))**2+abs(U(2))**2+abs(V(2))**2)
    K(1,2)=.5_dp*(abs(U(1))**2-abs(V(1))**2+abs(U(2))**2-abs(V(2))**2)
    K(1,3)=real(U(1)*conjg(V(1)))+real(U(2)*conjg(V(2)))
    K(1,4)=aimag(U(1)*conjg(V(1)))+aimag(U(2)*conjg(V(2)))

    K(2,1)=.5_dp*(abs(U(1))**2+abs(V(1))**2-abs(U(2))**2-abs(V(2))**2)
    K(2,2)=.5_dp*(abs(U(1))**2-abs(V(1))**2-abs(U(2))**2+abs(V(2))**2)
    K(2,3)=real(U(1)*conjg(V(1)))-real(U(2)*conjg(V(2)))
    K(2,4)=aimag(U(1)*conjg(V(1)))-aimag(U(2)*conjg(V(2)))

    K(3,1)=real(U(1)*conjg(U(2)))+real(V(1)*conjg(V(2)))
    K(3,2)=real(U(1)*conjg(U(2)))-real(V(1)*conjg(V(2)))
    K(3,3)=real(U(1)*conjg(V(2)))+real(U(2)*conjg(V(1)))
    K(3,4)=aimag(U(1)*conjg(V(2)))+aimag(U(2)*conjg(V(1)))

    K(4,1)=-aimag(U(1)*conjg(U(2)))-aimag(V(1)*conjg(V(2)))
    K(4,2)=-aimag(U(1)*conjg(U(2)))+aimag(V(1)*conjg(V(2)))
    K(4,3)=-aimag(U(1)*conjg(V(2)))+aimag(U(2)*conjg(V(1)))
    K(4,4)=real(U(1)*conjg(V(2)))-real(U(2)*conjg(V(1)))

    call mmmm(F2,K,F1)

    norm=F1(1,1)/F2(1,1)
    do j1 = 1, 4
        do j2 = 1, 4
            F2(j2,j1) = norm*F2(j2,j1)
        end do
    end do
       
end subroutine frotlc

