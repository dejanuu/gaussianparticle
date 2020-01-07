MODULE SIRISMESH

! For dealing with mesh-based geometries in SIRIS4-package.
!
! OBJ-import related functions are copied from John Burkardt's obj_io.f90 library.
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
  
  private :: obj_read
  
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  From John Burkardt's obj_io.f90 library
!
!! OBJ_READ reads graphics information from a Wavefront OBJ file.
!
!  Discussion:
!
!    It is intended that the information read from the file can
!    either start a whole new graphics object, or simply be added
!    to a current graphics object via the '<<' command.
!
!    This is controlled by whether the input values have been zeroed
!    out or not.  This routine simply tacks on the information it
!    finds to the current graphics object.
!
!  Example:
!
!    #  magnolia.obj
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 1.0 0.0
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    13 April 2001
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the name of the input file.
!
!    Input, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Input, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Input, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Input, integer ( kind = 4 ) ORDER_MAX, the maximum number of vertices
!    per face.
!
!    Output, real ( kind = 8 ) NODE_XYZ(3,NODE_NUM), the coordinates of points.
!
!    Output, integer ( kind = 4 ) FACE_ORDER(FACE_NUM), the number of vertices
!    per face.
!
!    Output, integer ( kind = 4 ) FACE_NODE(ORDER_MAX,FACE_NUM), the nodes
!    making faces.
!
!    Output, real ( kind = 8 ) NORMAL_VECTOR(3,NORMAL_NUM), normal vectors.
!
!    Output, integer ( kind = 4 ) VERTEX_NORMAL(ORDER_MAX,FACE_NUM), the indices
!    of normal vectors per vertex.
subroutine obj_read(input_filename, node_num, face_num, normal_num, &
  order_max, node_xyz, face_order, face_node, normal_vector, vertex_normal)

  character(len=*), intent(in) :: input_filename 
  integer, intent(in) :: node_num, face_num, normal_num, order_max
  real(kind=dp), dimension(3,node_num), intent(out) ::  node_xyz
  integer, dimension(face_num), intent(out) :: face_order
  integer, dimension(order_max,face_num), intent(out) :: face_node
  real(kind=dp), dimension(3,normal_num), intent(out) :: normal_vector
  integer, dimension(order_max,face_num), intent(out) :: vertex_normal
  
  logical :: done
  integer :: face, i, i1, i2, ierror, input_file_unit, &
    input_file_status, itemp, lchar, node, normal, text_num, &
    vertex, word_index
  real(kind=dp) :: temp
  character(len=255) :: line, word, word_one

  ierror = 0
  face = 0
  node = 0
  normal = 0
  text_num = 0

  face_node(1:order_max,1:face_num) = 0
  face_order(1:face_num) = 0
  node_xyz(1:3,1:node_num) = 0.0_dp
  normal_vector(1:3,1:normal_num) = 0.0_dp
  vertex_normal(1:order_max,1:face_num) = 0

  open(newunit=input_file_unit, file=input_filename, status='old', iostat=input_file_status)

  if(input_file_status /= 0) then
    write(error_unit, '(A)') ' '
    write(error_unit, '(A)') 'OBJ_READ - Fatal error!'
    write(error_unit, '(A,A,A)') "  Could not open the file '", trim(input_filename), "'."
    stop
  end if

  word = ' '
  !
  !  Read a line of text from the file.
  !
  do
    read(input_file_unit, '(A)', iostat=input_file_status) line

    if(input_file_status /= 0) then
      exit
    end if

    text_num = text_num + 1
    !
    !  Replace any control characters (in particular, TAB's) by blanks.
    !
    call s_control_blank(line)

    done = .true.
    word_index = 0
    !
    !  Read a word from the line.
    !
    call word_next_read(line, word, done)
    !
    !  If no more words in this line, read a new line.
    !
    if(done) then
      cycle
    end if
    !
    !  If this word begins with '#' or '$', then it's a comment.  Read a new line.
    !
    if(word(1:1) == '#' .or. word(1:1) == '$') then
      cycle
    end if

    word_index = word_index + 1

    if(word_index == 1) then
      word_one = word
    end if
    !  BEVEL
    !  Bevel interpolation.
    if(s_eqi(word_one, 'BEVEL')) then
    !  BMAT
    !  Basis matrix.
    else if(s_eqi( word_one, 'BMAT')) then
    !  C_INTERP
    !  Color interpolation.
    else if(s_eqi(word_one, 'C_INTERP')) then
    !  CON
    !  Connectivity between free form surfaces.
    else if(s_eqi(word_one, 'CON')) then
    !  CSTYPE
    !  Curve or surface type.
    else if(s_eqi(word_one, 'CSTYPE')) then
    !  CTECH
    !  Curve approximation technique.
    else if(s_eqi(word_one, 'CTECH')) then
    !  CURV
    !  Curve.
    else if(s_eqi(word_one, 'CURV')) then
    !  CURV2
    !  2D curve.
    else if(s_eqi(word_one, 'CURV2')) then
    !  D_INTERP
    !  Dissolve interpolation.
    else if(s_eqi(word_one, 'D_INTERP')) then
    !  DEG
    !  Degree.
    else if(s_eqi(word_one, 'DEG')) then
    !  END
    !  End statement.
    else if(s_eqi(word_one, 'END')) then
    !  F V1 V2 V3 ...
    !    or
    !  F V1/VT1/VN1 V2/VT2/VN2 ...
    !    or
    !  F V1//VN1 V2//VN2 ...
    !
    !  Face.
    !  A face is defined by the vertices.
    !  Optionally, slashes may be used to include the texture vertex
    !  and vertex normal indices.
    else if(s_eqi(word_one, 'F')) then

      face = face + 1

      vertex = 0

      do
        call word_next_read(line, word, done)
        if(done) exit

        vertex = vertex + 1

        !  Locate the slash characters in the word, if any.
        i1 = index(word, '/')
        if(0 < i1) then
          i2 = index(word(i1+1:), '/') + i1
        else
          i2 = 0
        end if

        !  Read the vertex index.
        call s_to_i4(word, itemp, ierror, lchar)

        if(ierror /= 0) then
          itemp = -1
          ierror = 0
          write(error_unit, '(A)') ' '
          write(error_unit, '(A)') 'OBJ_READ - Error!'
          write(error_unit, '(A)') '  Bad FACE field.'
          write(error_unit, '(A)') trim(word)
        end if

        face_node(vertex,face) = itemp
        face_order(face) = face_order(face) + 1

        !  If there are two slashes, then read the data following the second one.
        if(0 < i2) then
          call s_to_i4(word(i2+1:), itemp, ierror, lchar)
          vertex_normal(vertex,face) = itemp
        end if
        
      end do

    !  G
    !  Group name.
    else if(s_eqi(word_one, 'G')) then
    !  HOLE
    !  Inner trimming loop.
    else if(s_eqi(word_one, 'HOLE')) then
    !  L
    !  A line, described by a sequence of vertex indices.
    !  Are the vertex indices 0 based or 1 based?
    else if(s_eqi(word_one, 'L')) then
    !  LOD
    !  Level of detail.
    else if(s_eqi(word_one, 'LOD')) then
    !  MG
    !  Merging group.
    else if(s_eqi(word_one, 'MG')) then
    !  MTLLIB
    !  Material library.
    else if(s_eqi(word_one, 'MTLLIB')) then
    !  O
    !  Object name.
    else if(s_eqi(word_one, 'O')) then
    !  P
    !  Point.
    else if(s_eqi(word_one, 'P')) then
    !  PARM
    !  Parameter values.
    else if(s_eqi(word_one, 'PARM')) then
    !  S
    !  Smoothing group.
    else if(s_eqi(word_one, 'S')) then
    !  SCRV
    !  Special curve.
    else if(s_eqi(word_one, 'SCRV')) then
    !  SHADOW_OBJ
    !  Shadow casting.
    else if(s_eqi(word_one, 'SHADOW_OBJ')) then
    !  SP
    !  Special point.
    else if(s_eqi(word_one, 'SP')) then
    !  STECH
    !  Surface approximation technique.
    else if(s_eqi(word_one, 'STECH')) then
    !  STEP
    !  Stepsize.
    else if(s_eqi(word_one, 'STEP')) then
    !  SURF
    !  Surface.
    else if(s_eqi(word_one, 'SURF')) then
    !  TRACE_OBJ
    !  Ray tracing.
    else if(s_eqi(word_one, 'TRACE_OBJ')) then
    !  TRIM
    !  Outer trimming loop.
    else if(s_eqi(word_one, 'TRIM')) then
    !  USEMTL
    !  Material name.
    else if(s_eqi(word_one, 'USEMTL')) then
    !  V X Y Z
    !  Geometric vertex.
    else if(s_eqi(word_one, 'V')) then
      node = node + 1

      do i = 1, 3
        call word_next_read(line, word, done)
        call s_to_r8(word, temp, ierror, lchar)
        node_xyz(i,node) = temp
      end do
    !  VN
    !  Vertex normals.
    else if(s_eqi(word_one, 'VN')) then
      normal = normal + 1

      do i = 1, 3
        call word_next_read(line, word, done)
        call s_to_r8(word, temp, ierror, lchar)
        normal_vector(i,normal) = temp
      end do
    !  VT
    !  Vertex texture.
    else if(s_eqi(word_one, 'VT')) then
    !  VP
    !  Parameter space vertices.
    else if(s_eqi(word_one, 'VP')) then
    !  Unrecognized keyword.
    else
    end if

  end do

  close(input_file_unit)

end subroutine obj_read


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  From John Burkardt's obj_io.f90 library
!
!! OBJ_SIZE determines sizes of graphics objects in an Alias OBJ file.
!
!  Discussion:
!
!    The only items of interest to this routine are vertices,
!    faces, and normal vectors.
!
!  Example:
!
!    #  magnolia.obj
!
!    v -3.269770 -39.572201 0.876128
!    v -3.263720 -39.507999 2.160890
!    ...
!    v 0.000000 -9.988540 0.000000
!
!    vn 1.0 0.0 0.0
!    ...
!    vn 0.0 1.0 0.0
!
!    f 8 9 11 10
!    f 12 13 15 14
!    ...
!    f 788 806 774
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    15 December 2005
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, character ( len = * ) INPUT_FILENAME, the input file name.
!
!    Output, integer ( kind = 4 ) NODE_NUM, the number of points.
!
!    Output, integer ( kind = 4 ) FACE_NUM, the number of faces.
!
!    Output, integer ( kind = 4 ) NORMAL_NUM, the number of normal vectors.
!
!    Output, integer ( kind = 4 ) ORDER_MAX, the maximum face order.
subroutine obj_size(input_filename, node_num, face_num, normal_num, order_max)

  character(len=*), intent(in) :: input_filename
  integer, intent(out) :: node_num, face_num, normal_num, order_max

  logical :: done
  integer  :: i1, i2, ierror, input_file_status, input_file_unit, &
    itemp, lchar, text_num, vertex, word_index
  character(len=255) :: line, word, word_one

  ierror=0
  face_num=0
  node_num=0
  normal_num=0
  order_max=0
  text_num=0

  open(newunit=input_file_unit, file=input_filename, status='old', &
    iostat=input_file_status)

  if(input_file_status /= 0) then
    write(error_unit, '(A)') ' '
    write(error_unit, '(A)') 'OBJ_SIZE - Fatal error!'
    write(error_unit, '(A,A,A)') "  Could not open the file '",  &
      trim(input_filename), "'."
    stop
  end if

  word=' '

  !  Read a line of text from the file.
  do
    read(input_file_unit, '(A)', iostat=input_file_status) line
    if(input_file_status /= 0) then
      exit
    end if

    text_num=text_num + 1

    !  Replace any control characters (in particular, TAB's) by blanks.
    call s_control_blank(line)
    done=.true.
    word_index=0

    !  Read a word from the line.
    call word_next_read(line, word, done)

    !  If no more words in this line, read a new line.
    if(done) then
      cycle
    end if

    !  If this word begins with '#' or '$', then it's a comment.  Read a new line.
    if(word(1:1) == '#' .or. word(1:1) == '$') then
      cycle
    end if

    word_index=word_index + 1

    if(word_index == 1) then
      word_one=word
    end if

    !  F V1 V2 V3 ...
    !    or
    !  F V1/VT1/VN1 V2/VT2/VN2 ...
    !    or
    !  F V1//VN1 V2//VN2 ...
    !
    !  Face.
    !  A face is defined by the vertices.
    !  Optionally, slashes may be used to include the texture vertex
    !  and vertex normal indices.
    if(s_eqi(word_one, 'F')) then

      face_num=face_num + 1
      vertex=0
      do
        call word_next_read(line, word, done)
        if(done) then
          exit
        end if
        vertex=vertex + 1
        order_max=max(order_max, vertex)
        
        !  Locate the slash characters in the word, if any.
        i1=index(word, '/')
        if(0 < i1) then
          i2=index(word(i1+1:), '/') + i1
        else
          i2=0
        end if

        !  Read the vertex index.
        call s_to_i4(word, itemp, ierror, lchar)

        if(ierror /= 0) then
          itemp=-1
          ierror=0
          write(*, '(a)') ' '
          write(*, '(a)') 'OBJ_SIZE - Error!'
          write(*, '(a)') '  Bad FACE field.'
          write(*, '(a)') trim(word)
        end if

        !  If there are two slashes, then read the data following the second one.
        if(0 < i2) then
          call s_to_i4(word(i2+1:), itemp, ierror, lchar)
        end if

      end do

    !  V X Y Z W
    !  Geometric vertex.
    else if(s_eqi(word_one, 'V')) then
      node_num=node_num + 1
      cycle

    !  VN
    !  Vertex normals.
    else if(s_eqi(word_one, 'VN')) then
      normal_num=normal_num + 1
      cycle

    end if
  end do

  close(input_file_unit)

end subroutine obj_size


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Read mesh geometry from external OBJ file
subroutine read_obj(fn,XT,IT,NT,nnod,ntri,outstat)

  character(*), intent(in) :: fn
  real(kind=dp), dimension(:,:), pointer, intent(out) :: XT
  integer, dimension(:,:), pointer, intent(out) :: IT
  real(kind=dp), dimension(:,:), pointer, intent(out) :: NT
  integer, intent(out) :: nnod, ntri
  integer, intent(out) :: outstat
  integer :: nnorm, omax, cstat, j1, j2
  integer, dimension(:), allocatable :: FACE_ORDER
  integer, dimension(:,:), allocatable :: FACE_NODE, VERTEX_NORMAL
  real(kind=dp) :: r
  real(kind=dp), dimension(3) :: X1, X2
  real(kind=dp), dimension(:,:), allocatable :: NODE_XYZ, NORMAL_VECTOR

  outstat = -1

  ! Read mesh object size
  call obj_size(fn, nnod, ntri, nnorm, omax)
  if(omax > 3) then
    write(error_unit, '(A,A,I0)') "OBJ-file does not consist of triangles only, ", &
      "max no. of nodes is ", omax
    stop
  end if
  
  ! Allocate arrays
  allocate(XT(nnod,3), IT(ntri,3), NT(ntri,3), &
    NODE_XYZ(3,nnod), FACE_ORDER(ntri), FACE_NODE(3,ntri), &
    NORMAL_VECTOR(3,nnod), VERTEX_NORMAL(3,ntri), STAT=cstat)
  if(cstat /= 0) then
    write(error_unit,*) "Error in memory allocation of vertices and faces"
    stop
  end if
  
  call obj_read(fn, nnod, ntri, nnod, omax, NODE_XYZ, FACE_ORDER, &
    FACE_NODE, NORMAL_VECTOR, VERTEX_NORMAL)
    
  ! Data to SIRIS4-structures
  XT = transpose(NODE_XYZ)
  IT = transpose(FACE_NODE)
  
  ! Discard OBJ-normals and compute owns, similar to OFF-file read
  do j1=1,ntri
    do j2=1,3
      r=XT(IT(j1,1),j2)
      X1(j2)=XT(IT(j1,2),j2)-r
      X2(j2)=XT(IT(j1,3),j2)-r
    end do
    call provecn(NT(j1,:),X1,X2)
  end do
  
  deallocate(NODE_XYZ, FACE_ORDER, FACE_NODE, NORMAL_VECTOR, VERTEX_NORMAL)
  
  outstat = 0
  
end subroutine read_obj


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
    
  ! Read faces. Note that we only read triangles, and discard other nodes
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
