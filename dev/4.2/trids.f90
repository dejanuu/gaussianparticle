subroutine trids(MU,PHI,IT,nnod,ntri,ntr)

! TRIDS discretizes the spherical surface into altogether ntri=8*ntr**2
! triangles. It stores the nnod=4*ntr**2+2 nodes and right-handed node
! addresses for each triangle. ntr is the number of triangle rows in an
! octant.


use common
implicit none

integer  :: nnod,ntri,ntr,j0,j1,j2,j3
integer , dimension(260000,3) :: IT
integer , dimension(0:360,0:720) :: NJJ
real (kind=dp) :: the,fi,ct,st,cf,sf
real (kind=dp), dimension(130000,3) :: U
real (kind=dp), dimension(130000) :: MU, PHI

! NODES:

! Upper hemisphere including equator:

nnod=1
U(nnod,1)=0.0_dp
U(nnod,2)=0.0_dp
U(nnod,3)=1.0_dp
MU(nnod)=1.0_dp
PHI(nnod)=0.0_dp
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
        U(nnod,1)=st*cf
        U(nnod,2)=st*sf
        U(nnod,3)=ct
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
        U(nnod,1)=st*cf
        U(nnod,2)=st*sf
        U(nnod,3)=ct
        MU(nnod)=ct
        PHI(nnod)=fi
        NJJ(2*ntr-j1,j2)=nnod
        if (j2.eq.0) NJJ(2*ntr-j1,4*j1)=nnod
    end do
end do

nnod=nnod+1
U(nnod,1)=0.0_dp
U(nnod,2)=0.0_dp
U(nnod,3)=-1.0_dp
MU(nnod)=-1.0_dp
PHI(nnod)=0.0_dp
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


end subroutine trids






