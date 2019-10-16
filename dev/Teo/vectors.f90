module vectors

! Defines a simple vector data type and some basic operations

use parameters

implicit none

type vector
    real(kind=dp) :: x,y,z
end type vector


interface operator(+)
    procedure vadd
end interface operator(+)

interface operator(-)
    procedure vdiff
end interface operator(-)

interface operator(*)
    procedure vmult
end interface operator(*)

interface operator(/)
    procedure vdiv
end interface operator(/)

interface operator(.cross.)
    procedure vprod
end interface operator(.cross.)

interface operator(.dot.)
    procedure sprod
end interface operator(.dot.)


contains

! The addition of two vectors

type(vector) function vadd(v1,v2)

implicit none
type(vector),intent(in) :: v1,v2

vadd%x = v1%x+v2%x
vadd%y = v1%y+v2%y
vadd%z = v1%z+v2%z

end function


! The difference between two vectors.

type(vector) function vdiff(v1,v2)

implicit none
type(vector),intent(in) :: v1,v2

vdiff%x = v1%x-v2%x
vdiff%y = v1%y-v2%y
vdiff%z = v1%z-v2%z

end function

! The division of a vector by a real number.

type(vector) function vdiv(v1,r1)

implicit none
type(vector),intent(in) :: v1
real(kind=dp),intent(in) :: r1

vdiv%x = v1%x/r1
vdiv%y = v1%y/r1
vdiv%z = v1%z/r1

end function

! The multiplication of a vector by a real number.

type(vector) function vmult(v1,r1)

implicit none
type(vector),intent(in) :: v1
real(kind=dp),intent(in) :: r1

vmult%x = v1%x*r1
vmult%y = v1%y*r1
vmult%z = v1%z*r1

end function

! The angle between two vectors in [0,pi]

real(kind=dp) function angle1(v1,v2)

implicit none
type(vector),intent(in) :: v1,v2
real(kind=dp) :: sprn

sprn = v1 .dot. v2
angle1 = acos(sprn/(length(v1)*length(v2)))

end function


! The angle between two vectors in [0,pi]

real(kind=dp) function angle2(v1,v2)

implicit none
type(vector),intent(in) :: v1,v2
real(kind=dp) :: vprn

vprn = length(v1 .cross. v2)
angle2 = asin(vprn/(length(v1)*length(v2)))

end function


! Returns the scalar product of a vector

real(kind=dp) function sprod(v1,v2)

implicit none
type(vector),intent(in) :: v1,v2

sprod = v1%x*v2%x + v1%y*v2%y + v1%z*v2%z

end function


! Returns the vector product of two vectors

type(vector) function vprod(v1,v2)

type(vector),intent(in) :: v1,v2

vprod%x = v1%y*v2%z - v1%z*v2%y
vprod%y = v1%z*v2%x - v1%x*v2%z
vprod%z = v1%x*v2%y - v1%y*v2%x

end function


! Returns the length of a vector

real(kind=dp) function length(v1)

implicit none
type(vector) :: v1

length = sqrt(v1%x**2+v1%y**2+v1%z**2)

end function length



end module
