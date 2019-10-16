module objsupport
use iso_fortran_env, only : REAL64
use obj_io


implicit none
integer, parameter :: rk = REAL64

private rk
private simplex3volume,compute_volume
contains


function simplex3volume(u,v,b) result(A)
    real(kind=rk),intent(in) :: u(3),v(3),b(3)
    real(kind=rk) :: A
    A=(u(2)*v(3)-u(3)*v(2))*b(1)
    A=A+(u(3)*v(1)-u(1)*v(3))*b(2)
    A=A+(u(1)*v(2)-u(2)*v(1))*b(3)
    A=1.0_rk/6.0_rk*abs(A)
end function



function compute_volume(XT,IT,ntri) result(volume)
    real(kind=rk), intent(in) :: XT(:,:)
    integer, intent(in) :: IT(:,:),ntri
    real(kind=rk) :: volume
    real(kind=rk) :: O(3),B(3),C(3),D(3)
    integer :: j1,j2
    O=(/0.0_rk,0.0_rk,0.0_rk/)
    volume=0.0_rk
    do j1=1,ntri
        B=XT(IT(j1,1),:)
        C=XT(IT(j1,2),:)
        D=XT(IT(j1,3),:)
        volume=volume+simplex3volume(B,D,C)
    end do
end function



subroutine readMesh2(XN,NT,IT,ntri,nnod,rmax,fname,newV)
    integer, intent(out) :: ntri,nnod
    integer, intent(out) :: IT(:,:)
    real(kind=rk), intent(in) :: newV
    real(kind=rk), intent(out) :: XN(:,:),NT(:,:)
    real(kind=rk), intent(out) :: rmax
    real(kind=rk), allocatable :: XN0(:,:),NT0(:,:)
    integer, allocatable :: IT0(:,:)
    real(kind=rk) :: volume1,volume2,rScale
    character(32), intent(in) :: fname
    character(64) :: cmd
    integer :: vertices,normals,faces
    integer,allocatable :: face_order(:),vertex_normal(:,:)
    integer :: tmp1,tmp2,tmp3,i
    
    real(kind=rk) :: r,X1(3),X2(3),X3(3)
    integer :: j1,j2


    write(cmd,*) "python test.py ",trim(adjustl(fname))," > meshstat"
    write(6,*) cmd
    call execute_command_line(cmd) 
    open(unit=1,file="meshstat")
    read(1,*) tmp1,tmp2,tmp3
    close(1)

    
    XN = 0.0_rk
    NT = 0.0_rk
    IT = 0
    allocate(face_order(tmp3))
    allocate(XN0(3,tmp1),NT0(3,tmp2),IT0(3,tmp3))
    allocate(vertex_normal(3,tmp3))

    call obj_read(trim(adjustl(fname)),tmp1, tmp3, tmp2, &
    3, XN0, face_order, IT0, NT0, vertex_normal )
    deallocate(face_order)
    deallocate(vertex_normal)
    
    

    
    XN(1:tmp1,:)=transpose(XN0)
    IT(1:tmp3,:)=transpose(IT0)
 


   




    ntri=tmp3
    nnod=tmp1
   
    do j1=1,ntri
        do j2=1,3
            r=XN(IT(j1,1),j2)
            X1(j2)=XN(IT(j1,2),j2)-r
            X2(j2)=XN(IT(j1,3),j2)-r
        end do
        call provec(X3,X1,X2)
        r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
        do j2=1,3
            NT(j1,j2)=X3(j2)/r
        end do
    end do 

    volume1 = compute_volume(XN,IT,ntri)
    rScale = (newV/volume1)**(1.0_rk/3.0_rk)
    XN(:,:) = XN(:,:)*rScale
    volume2 = compute_volume(XN,IT,ntri)


    rmax = -1;
    do i=1,tmp1
        rmax = max(rmax,dot_product(XN(i,:),XN(i,:)))
    enddo
    rmax = sqrt(rmax)
    !write(6,*) rmax




    write(6,*) volume1,volume2,rmax


end subroutine





subroutine readMesh(XN,NT,IT,ntri,nnod,rmax)
    integer, intent(out) :: ntri,nnod
    integer, intent(out) :: IT(:,:)
    real(kind=rk), intent(out) :: XN(:,:),NT(:,:)
    real(kind=rk), intent(out) :: rmax
    real(kind=rk), allocatable :: XN0(:,:),NT0(:,:)
    integer, allocatable :: IT0(:,:)

    character(32) :: fname
    integer :: vertices,normals,faces
    integer,allocatable :: face_order(:),vertex_normal(:,:)
    integer :: tmp1,tmp2,tmp3,i
    
    real(kind=rk) :: r,X1(3),X2(3),X3(3)
    integer :: j1,j2


    write(6,*) "fname"
    read(5,*) fname    
    write(6,*) "vertices"
    read(5,*) tmp1
    write(6,*) "normals"
    read(5,*) tmp2
    write(6,*) "faces"
    read(5,*) tmp3 
    XN = 0.0_rk
    NT = 0.0_rk
    IT = 0
    allocate(face_order(tmp3))
    allocate(XN0(3,tmp1),NT0(3,tmp2),IT0(3,tmp3))
    allocate(vertex_normal(3,tmp3))

    call obj_read(trim(adjustl(fname)),tmp1, tmp3, tmp2, &
    3, XN0, face_order, IT0, NT0, vertex_normal )
    deallocate(face_order)
    deallocate(vertex_normal)
    
    

    
    XN(1:tmp1,:)=transpose(XN0)
    IT(1:tmp3,:)=transpose(IT0)
 


   



    rmax = -1;
    do i=1,tmp1
        rmax = max(rmax,dot_product(XN(i,:),XN(i,:)))
    enddo
    rmax = sqrt(rmax)
    write(6,*) rmax

    ntri=tmp3
    nnod=tmp1
   
    do j1=1,ntri
        do j2=1,3
            r=XN(IT(j1,1),j2)
            X1(j2)=XN(IT(j1,2),j2)-r
            X2(j2)=XN(IT(j1,3),j2)-r
        end do
        call provec(X3,X1,X2)
        r=sqrt(X3(1)**2+X3(2)**2+X3(3)**2)
        do j2=1,3
            NT(j1,j2)=X3(j2)/r
        end do
    end do



 
end subroutine



end module
