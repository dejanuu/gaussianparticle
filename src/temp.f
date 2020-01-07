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
  integer  :: i, i1, i2, ierror, input_file_status, input_file_unit, &
    itemp, lchar, text_num, vertex, word_index
  real(kind=dp) temp
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
!
!  VN
!  Vertex normals.
!
    else if(s_eqi(word_one, 'VN')) then

      normal_num=normal_num + 1
      cycle

    end if

  end do

  close(unit=input_file_unit)

  return
end