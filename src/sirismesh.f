MODULE SIRISMESH

! Big parts of this module are functions copied from John Burkardt's obj_io.f90 library.
!
! v2019-12-31
!
! Karri Muinonen, Timo Väisänen, Antti Penttilä
! Department of Physics, University of Helsinki, Finland

  use sirisconstants
  use sirisutils
  use sirismath
  use sirisgeometry
  
  public
  
contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read mesh geometry from external OFF file
subroutine read_off(fn,XT,IT,NT,nnod,ntri,outstat)

  character(*), intent(in) :: fn
  real(kind=dp), dimension(:,:), pointer, intent(out) :: XT
  integer, dimension(:,:), pointer, intent(out) :: IT
  real(kind=dp), dimension(:,:), pointer, intent(out) :: NT
  integer, intent(out) :: nnod, ntri
  integer, intent(out) :: outstat
  integer :: j1, j2, fu, slen, cstat
  real(kind=dp) :: r
  real(kind=dp), dimension(3) :: X1, X2
  character(len=line_str_length) :: line, sl

  outstat = -1
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='read', status='old')
  
  ! Scan for keyword 'OFF'
  do
    read(fu,'(A)',IOSTAT=cstat) line
    if(cstat /= 0) then
      write(error_unit,*) "File ended while scanning keyword 'OFF'"
      stop
    end if
    call strip_string(line,sl,slen)
    j2 = index(sl,"OFF")
    if(j2==1) then
      exit
    end if
  end do
  
  ! Scan for no. of vertices, and faces
  do
    read(fu,'(A)',IOSTAT=cstat) line
    if(cstat /= 0) then
      write(error_unit,*) "File ended while scanning for no. of vertices and faces"
      stop
    end if
    call strip_string(line,sl,slen)
    if(slen>0) then
      read(sl,*,IOSTAT=cstat) nnod, ntri
      if(cstat /= 0 .or. nnod <= 0 .or. ntri <= 0) then
        write(error_unit,*) "Wrong format when scanning for no. of vertices and faces"
        stop
      end if
      exit
    end if
  end do
  
  ! Allocate arrays
  allocate(XT(nnod,3), IT(ntri,3), NT(ntri,3), STAT=cstat)
  if(cstat /= 0) then
    write(error_unit,*) "Error in memory allocation of vertices and faces"
    stop
  end if
  
  ! Read vertices
  do j1=1,nnod
    read(fu,*,IOSTAT=cstat) XT(j1,:)
    if(cstat /= 0) then
      write(error_unit,*) "Error reading vertice no. ", j1
      stop
    end if
  end do
    
  ! Read faces. Note that we only read triangles, and discard other vertices
  do j1=1,ntri
    read(fu,*,IOSTAT=cstat) j2, IT(j1,:)
    if(cstat /= 0) then
      write(error_unit,*) "Error reading face no. ", j1
      stop
    end if
    IT(j1,1) = IT(j1,1)+1
    IT(j1,2) = IT(j1,2)+1
    IT(j1,3) = IT(j1,3)+1
  end do
  
  ! Outer unit triangle normals:
  do j1=1,ntri
    do j2=1,3
      r=XT(IT(j1,1),j2)
      X1(j2)=XT(IT(j1,2),j2)-r
      X2(j2)=XT(IT(j1,3),j2)-r
    end do
    call provecn(NT(j1,:),X1,X2)
  end do
  
  close(fu)
  
  outstat = 0

end subroutine read_off


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_idl(fbn,XT,IT,nnod,ntri)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri
  integer :: j1, fu
  character(len=file_name_length) :: fn

  ! File name
  write(fn, '(A,A)') trim(fbn), ".idf"
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='write', status='replace')

  ! Write number of vertices and triangles
  write(fu,*) nnod, ntri
  ! Vertices
  do j1=1,nnod
    write(fu,*) XT(j1,:)
  end do
  ! Triangle indices
  do j1=1,ntri
    write (fu,*) 3
    write (fu,*) IT(j1,:)
  end do

  close(fu)

end subroutine save_idl


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_matlab(fbn,XT,nnod)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, intent(in) :: nnod
  integer :: j1, fux, fuy, fuz
  character(len=file_name_length) :: fnx, fny, fnz

  ! File names
  write(fnx, '(A,A)') trim(fbn), "x.out"
  write(fny, '(A,A)') trim(fbn), "y.out"
  write(fnz, '(A,A)') trim(fbn), "z.out"
  
  ! File units and unit opening
  open(newunit=fux, file=trim(fnx), action='write', status='replace')
  open(newunit=fuy, file=trim(fny), action='write', status='replace')
  open(newunit=fuz, file=trim(fnz), action='write', status='replace')

  ! Write surface vertices
  do j1=1,nnod
    write (fux,*) XT(j1,1)
    write (fuy,*) XT(j1,2)
    write (fuz,*) XT(j1,3)
  end do

  close(fux)
  close(fuy)
  close(fuz)
   
end subroutine save_matlab


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_off(fbn,XT,IT,nnod,ntri)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri
  integer :: j1, fu
  character(len=file_name_length) :: fn

  ! File name
  write(fn, '(A,A)') trim(fbn), ".off"
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='write', status='replace')

  write(fu,'(A)') "OFF"
  write(fu,'(I0,1X,I0,1X,I0)') nnod, ntri, 0

  ! Write surface vertices
  do j1=1,nnod
    write (fu,*) XT(j1,:)
  end do
  ! Triangle indices
  do j1=1,ntri
    write (fu,*) 3, IT(j1,:)-1
  end do
  close(fu)

end subroutine save_off


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
subroutine save_vtk(fbn,XT,IT,nnod,ntri)

  character(*), intent(in) :: fbn
  real(kind=dp), dimension(:,:), intent(in) :: XT
  integer, dimension(:,:), intent(in) :: IT
  integer, intent(in) :: nnod, ntri
  integer :: j1, fu
  character(len=file_name_length) :: fn

  ! File name
  write(fn, '(A,A)') trim(fbn), ".vtk"
  
  ! File unit and unit opening
  open(newunit=fu, file=trim(fn), action='write', status='replace')

  write(fu,'(A)') '# vtk DataFile Version 2.0'
  write(fu,'(A)') 'gsphere output            '
  write(fu,'(A)') 'ASCII                     '
  write(fu,'(A)') 'DATASET POLYDATA          '
  write(fu,'(A,I0,A)') 'POINTS ',nnod,' float'

  ! Write surface vertices
  do j1=1,nnod
    write (fu,*) XT(j1,:)
  end do
  
  ! Triangle indices
  write (fu,'(A,I7,I7)') 'POLYGONS ',ntri,4*ntri
  do j1=1,ntri
    write (fu,*) 3,IT(j1,:)-1
  end do
  close(fu)

end subroutine save_vtk


END MODULE SIRISMESH
