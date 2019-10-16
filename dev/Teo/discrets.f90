module discrets
! Discretization:
!
! SPHDS: spherical coordinates 
! TRIDSA: triangles, allocatable
! TRIDS:  triangles, static
!
! Required modules:
!
use parameters          ! Parameters
        
        
contains

subroutine SPHDS(MU,PHI,nthe,nphi)

! SPHDS discretizes the spherical surface into a polar angle -azimuth angle
! grid. Version 2002-12-16.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
integer,intent(in) :: nthe,nphi
real(kind=dp),intent(inout) :: MU(0:180),PHI(0:360)
integer :: j1,j2
real(kind=dp) :: dthe,dphi,pi

pi=4.0d0*atan(1.0d0)
dthe=pi/nthe
dphi=2.0d0*pi/nphi

do j1=0,nthe
 MU(j1)=cos(j1*dthe)
 do j2=0,nphi
  PHI(j2)=(j2+0.5d0)*dphi
 end do
end do
end subroutine SPHDS



subroutine TRIDSA(MU,PHI,IT,ntr)

! TRIDSA discretizes the spherical surface into altogether ntri=8*ntr**2 
! triangles. It stores the nnod=4*ntr**2+2 nodes and right-handed node 
! addresses for each triangle. ntr is the number of triangle rows in an
! octant. Version 2017-05-23. Uses dynamically allocated arrays.
!
! Copyright (C) 2002 Karri Muinonen

implicit none
integer,intent(in) :: ntr
integer :: jtri
integer,intent(inout),allocatable :: IT(:,:)
real(kind=dp),intent(inout),allocatable :: MU(:),PHI(:)
integer :: NJJ(0:360,0:720),j0,j1,j2,j3,jnod          ! NJJ maybe into dyn-array
real(kind=dp) :: the,fi,ct,st,cf,sf,pi
!real(kind=dp),allocatable :: U(:,:) !U uses max (130000,3)

pi=4.0d0*atan(1.0d0)

!allocate(U(4*ntr**2+2,3))

! NODES:

! Upper hemisphere including equator:
if(.not. (allocated(IT) .or. allocated(MU) .or. allocated(PHI))) &
   stop 'trouble in TRIDS: arrays not allocated'

jnod=1
MU(jnod)=1.0d0
PHI(jnod)=0.0d0
NJJ(0,0)=jnod

do j1=1,ntr
 the=j1*pi/(2*ntr)
 ct=cos(the)
 st=sin(the)

 do j2=0,4*j1-1
  fi=j2*pi/(2*j1)
  cf=cos(fi)
  sf=sin(fi)

  jnod=jnod+1

  MU(jnod)=ct
  PHI(jnod)=fi
  NJJ(j1,j2)=jnod
  if (j2.eq.0) NJJ(j1,4*j1)=jnod
 end do
end do

! Lower hemisphere excluding equator:

do j1=ntr-1,1,-1
 the=(2*ntr-j1)*pi/(2*ntr)
 ct=cos(the)
 st=sin(the)

 do j2=0,4*j1-1
  fi=j2*pi/(2*j1)
  cf=cos(fi)
  sf=sin(fi)

  jnod=jnod+1
  MU(jnod)=ct
  PHI(jnod)=fi
  NJJ(2*ntr-j1,j2)=jnod
  if (j2.eq.0) NJJ(2*ntr-j1,4*j1)=jnod
 end do
end do

jnod=jnod+1

MU(jnod)=-1.0d0
PHI(jnod)=0.0d0
NJJ(2*ntr,0)=jnod
if (jnod.ne.4*ntr**2+2) &
stop 'Trouble in TRIDS: number of nodes inconsistent.'


! TRIANGLES:

! Upper hemisphere: 

jtri=0
do j1=1,ntr
 do j3=1,4
  j0=(j3-1)*j1
  
  jtri=jtri+1
  IT(jtri,1)=NJJ(j1-1,j0-(j3-1))
  IT(jtri,2)=NJJ(j1,  j0       )
  IT(jtri,3)=NJJ(j1,  j0+1     )
                     
  do j2=j0+1,j0+j1-1
   jtri=jtri+1
   IT(jtri,1)=NJJ(j1,  j2         )
   IT(jtri,2)=NJJ(j1-1,j2  -(j3-1))
   IT(jtri,3)=NJJ(j1-1,j2-1-(j3-1))

   jtri=jtri+1
   IT(jtri,1)=NJJ(j1-1,j2-(j3-1)  )
   IT(jtri,2)=NJJ(j1,  j2         )
   IT(jtri,3)=NJJ(j1,  j2+1       )
  end do
 end do
end do

! Lower hemisphere: 

do j1=ntr+1,2*ntr
 do j3=1,4
  j0=(j3-1)*(2*ntr-j1)
  
  jtri=jtri+1
  IT(jtri,1)=NJJ(j1,  j0         )
  IT(jtri,2)=NJJ(j1-1,j0+1+(j3-1))
  IT(jtri,3)=NJJ(j1-1,j0  +(j3-1))
                     
  do j2=j0+1,j0+(2*ntr-j1)
   jtri=jtri+1
   IT(jtri,1)=NJJ(j1,  j2         )
   IT(jtri,2)=NJJ(j1-1,j2+(j3-1)  )
   IT(jtri,3)=NJJ(j1,  j2-1       )
   
   jtri=jtri+1
   IT(jtri,1)=NJJ(j1,  j2         )
   IT(jtri,2)=NJJ(j1-1,j2+1+(j3-1))
   IT(jtri,3)=NJJ(j1-1,j2  +(j3-1))
  end do
 end do
end do

if (jtri.ne.8*ntr**2) &
stop 'Trouble in TRIDSA: number of triangles inconsistent.'
end subroutine TRIDSA



subroutine TRIDS(MU,PHI,IT,nnod,ntri,ntr)

! TRIDS discretizes the spherical surface into altogether ntri=8*ntr**2 
! triangles. It stores the nnod=4*ntr**2+2 nodes and right-handed node 
! addresses for each triangle. ntr is the number of triangle rows in an
! octant. 

! Author: Johanna Torppa (based on code by Karri Muinonen)
! Version: 2015-07-13

implicit none
integer::IT(260000,3),NJJ(0:360,0:720),nnod,ntri,ntr,j0,j1,j2,j3 !,j4
real(kind=dp)::MU(130000),PHI(130000),the,fi,ct,st,cf,sf
!real(kind=dp)::U(130000,3)

! NODES:

! Upper hemisphere including equator:

nnod=1
!U(nnod,1)=0.0d0
!U(nnod,2)=0.0d0
!U(nnod,3)=1.0d0
MU(nnod)=1.0d0
PHI(nnod)=0.0d0
NJJ(0,0)=nnod

do j1=1,ntr
    the=j1*pi/(2*ntr)
    ct=cos(the)
    st=sin(the)

    do j2=0,4*j1-1
        fi=j2*pi/(2*j1)
        cf=cos(fi)
        sf=sin(fi)

        nnod=nnod+1
!        U(nnod,1)=st*cf
!        U(nnod,2)=st*sf
!        U(nnod,3)=ct
        MU(nnod)=ct
        PHI(nnod)=fi
        NJJ(j1,j2)=nnod
        if (j2.eq.0) NJJ(j1,4*j1)=nnod
    end do
end do

! Lower hemisphere excluding equator:

do j1=ntr-1,1,-1
    the=(2*ntr-j1)*pi/(2*ntr)
    ct=cos(the)
    st=sin(the)

    do j2=0,4*j1-1
        fi=j2*pi/(2*j1)
        cf=cos(fi)
        sf=sin(fi)

        nnod=nnod+1
!        U(nnod,1)=st*cf
!        U(nnod,2)=st*sf
!        U(nnod,3)=ct
        MU(nnod)=ct
             PHI(nnod)=fi
        NJJ(2*ntr-j1,j2)=nnod
        if (j2.eq.0) NJJ(2*ntr-j1,4*j1)=nnod
    end do
end do

nnod=nnod+1
!U(nnod,1)=0.0d0
!U(nnod,2)=0.0d0
!U(nnod,3)=-1.0d0
MU(nnod)=-1.0d0
PHI(nnod)=0.0d0
NJJ(2*ntr,0)=nnod
       
if (nnod.ne.4*ntr**2+2) stop 'Trouble in TRIDS: number of nodes inconsistent.'


! TRIANGLES:

! Upper hemisphere: 

ntri=0
do j1=1,ntr
    do j3=1,4
        j0=(j3-1)*j1
         
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1-1,j0-(j3-1))
        IT(ntri,2)=NJJ(j1,  j0       )
        IT(ntri,3)=NJJ(j1,  j0+1     )
                            
        do j2=j0+1,j0+j1-1
            ntri=ntri+1
            IT(ntri,1)=NJJ(j1,  j2         )
            IT(ntri,2)=NJJ(j1-1,j2  -(j3-1))
            IT(ntri,3)=NJJ(j1-1,j2-1-(j3-1))

            ntri=ntri+1
            IT(ntri,1)=NJJ(j1-1,j2-(j3-1)  )
            IT(ntri,2)=NJJ(j1,  j2         )
            IT(ntri,3)=NJJ(j1,  j2+1       )
        end do
    end do
end do

! Lower hemisphere: 

do j1=ntr+1,2*ntr
    do j3=1,4
        j0=(j3-1)*(2*ntr-j1)
         
        ntri=ntri+1
        IT(ntri,1)=NJJ(j1,  j0         )
        IT(ntri,2)=NJJ(j1-1,j0+1+(j3-1))
        IT(ntri,3)=NJJ(j1-1,j0  +(j3-1))
                            
        do j2=j0+1,j0+(2*ntr-j1)
            ntri=ntri+1
            IT(ntri,1)=NJJ(j1,  j2         )
            IT(ntri,2)=NJJ(j1-1,j2+(j3-1)  )
            IT(ntri,3)=NJJ(j1,  j2-1       )
          
            ntri=ntri+1
            IT(ntri,1)=NJJ(j1,  j2         )
            IT(ntri,2)=NJJ(j1-1,j2+1+(j3-1))
            IT(ntri,3)=NJJ(j1-1,j2  +(j3-1))
        end do
    end do
end do

if (ntri.ne.8*ntr**2) stop 'Trouble in TRIDS: number of triangles inconsistent.'

end subroutine TRIDS

end module





