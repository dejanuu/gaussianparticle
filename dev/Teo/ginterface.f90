module ginterface
use corrfunc     ! Correlation functions.
use discrets     ! Discretization.
use gsaxg        ! Axisymmetric G-sphere generator.
use gsg          ! G-sphere generator.
use minko       ! Solves Minkowski problem, used for convex spheres.
use chull       ! Finds the convex hull of a sphere using the 3-D gift wrapping algorithm.
use parameters 
use rng
implicit none


real(kind=dp),allocatable :: CSCF(:),SCFSTD(:,:),ACF(:,:),BCF(:,:)
real(kind=dp) :: CEU(3),SEU(3),a,sig,beta,gami,elli,nuc,rd,x
integer :: nss,lmin,lmax
integer :: gflg,ntr,nnod,ntri
real(kind=dp),allocatable :: MUT(:),PHIT(:),XT(:,:),NT(:,:)
integer,allocatable :: IT(:,:)
integer :: indx
character(32) :: outfile
private indx

contains

subroutine init_ginterface(infile)
    integer :: j0,j1,j2,j3,dflg,cflg,fflg,convflg,seed

    integer :: nthe,nphi

    character :: ca*2
    character(*), intent(in) :: infile

    logical :: there
    integer :: irnd
    common irnd

    ! Initializations:

    rd=pi/180.0d0

    irnd=5
           
    ! Input file specified in command line argument:
    inquire(file=infile,exist=there)
    if (.not. there) stop &
     'Trouble in GSPHERE: no input file specified in argument!'

    ! Input parameters from option file:

    open(unit=1, file=trim(infile), status='old')

    a=1.0d0
    read (1,30) ca
    read (1,20) gflg     ! General (1) or axisymmetric spheres (2).
    read (1,20) dflg     ! Spherical-coord. (1) or triangle (2) or spherical
                         ! harmonic coefficient (3) representation.
    read (1,20) convflg  ! 1): Nonconvex sphere, 2) the convex hull of the sphere, 3)
                         ! convex sphere with Minkowski solver, 4) Minko with coefficients from
                         ! mcmc1.out.
    read (1,20) seed     ! Seed for the random number generator.
    read (1,20) cflg     ! Cor. function (C_1=power law, C_2=Gauss, C_3=file).
    read (1,10) sig      ! Relative standard deviation of radial distance.
    read (1,10) nuc      ! Power law index for C_1 correlation.
    read (1,10) gami     ! Input angle for C_2 correlation.
    read (1,20) lmin     ! Minimum degree in C_1, C_2, C_3.
    read (1,20) lmax     ! Maximum degree in C_1, C_2, C_3.
    read (1,20) nthe     ! Discretization: number of polar angles.
    read (1,20) nphi     ! Discretization: number of azimuths.
    read (1,20) ntr      ! Discretization: number of triangle rows per octant.
    read (1,20) nss      ! Sphere identification number.
    read (1,20) fflg     ! Matlab (1), vtk (2), or idf (3) output.
    read (1,40) outfile  ! Output file name.

    10     format (E12.6)
    20     format (I12)
    30     format (/A2/)
    40     format (A32)
    close(unit=1)


    ! Input check:
           
    if (gflg.ne.1 .and. gflg.ne.2) stop &
     'Trouble in GSPHERE: general or axisymmetri! spheres.'
    if (dflg.ne.1 .and. dflg.ne.2 .and. dflg.ne.3 .and. dflg.ne.4) stop &
     'Trouble in GSPHERE: spherical or triangle discretization.'
    if (cflg.ne.1 .and. cflg.ne.2 .and. cflg.ne.3 .and. cflg.ne.4) stop &
     'Trouble in GSPHERE: correlation function unknown.'

    if (convflg .lt. 1 .or. convflg .gt. 4) stop &
     'Trouble in GSPHERE: need to specify convexity (1,2, or 3).'

    if (sig.le.0.0d0) stop &
     'Trouble in GSPHERE: standard deviation .le. 0.'

    if (cflg.eq.2) then
     if (gami.le.0.0d0 .or. gami.gt.180.0d0) stop &
      'Trouble in GSPHERE: input angle .le. 0. .or.  .gt. 180'
     if (lmin.gt.0 .or. lmax.lt.int(300.0d0/gami)) then
      print*,'Warning in GSPHERE: correlation angle will differ '
      print*,'from input value. Set minimum degree to 0 and '
      print*,'maximum degree .gt. (300 deg)/(input value).'
     endif
    endif

    if (lmax.gt.256) stop &
     'Trouble in GSPHERE: maximum degree .gt. 256.'
    if (lmin.lt.0) stop &
     'Trouble in GSPHERE: minimum degree .lt. 0.'
    if (lmin.gt.lmax) stop &
     'Trouble in GSPHERE: minimum degree .lt. maximum degree.'
    if (cflg.eq.1 .and. lmin.lt.2) stop &
     'Trouble in GSPHERE: minimum degree .lt.2.'

    if (nthe.gt.180) stop &
     'Trouble in GSPHERE: number of polar angles .gt.180.'
    if (nphi.gt.360) stop &
     'Trouble in GSPHERE: number of azimuths .gt.360.'
    if (ntr.gt.180) stop &
     'Trouble in GSPHERE: number of triangle rows .gt.180.'

    if (nss.le.0) stop &
     'Trouble in GSPHERE: sphere identification number .lt. 0.'

    if (fflg.le.0 .or. fflg.ge.4) stop &
     'Trouble in GSPHERE: output format not specified properly'

    ! Miscellaneous:

    gami=gami*rd
    elli=2.0d0*sin(0.5d0*gami)
    
    ! Init random number generator
    CALL init_random_seed(seed)

    ! Initialization of the Gaussian random sphere:

    beta=sqrt(log(sig**2+1.0d0))
    allocate(CSCF(0:lmax))


    if (cflg.eq.1) then
     call CS1CF(CSCF,nuc,lmin,lmax)
    elseif (cflg.eq.2) then
     call CS2CF(CSCF,elli,lmin,lmax)
    elseif (cflg.eq.3) then
     call CS3CF(CSCF,lmin,lmax)
    elseif (cflg.eq.4) then
     CSCF = 1.
    endif


    do j1=lmin,lmax
     if (CSCF(j1).lt.0.0d0) stop &
      'Trouble in GSPHERE: negative Legendre coefficient.'
    end do

           
    allocate(SCFSTD(0:lmax,0:lmax))
    call SGSCFSTD(SCFSTD,CSCF,beta,lmin,lmax)
    nnod=1
    do j1=1,ntr
        do j2=0,4*j1-1
            nnod=nnod+1
        end do
    end do
    do j1=ntr-1,1,-1
        do j2=0,4*j1-1
            nnod=nnod+1
        end do
    end do
    nnod=nnod+1


    ntri=0
    do j1=1,ntr
     do j3=1,4
      ntri=ntri+1           
      do j2=j0+1,j0+j1-1
       ntri=ntri+1
       ntri=ntri+1
      end do
     end do
    end do

    ! Lower hemisphere: 

    do j1=ntr+1,2*ntr
     do j3=1,4
      ntri=ntri+1
      do j2=j0+1,j0+(2*ntr-j1)
       ntri=ntri+1
       ntri=ntri+1
      end do
     end do
    end do
    


    
    ! Triangle representation for general and axisymmetric shapes:
    allocate(IT(8*ntr**2,3),PHIT(4*ntr**2+2),MUT(4*ntr**2+2))
    allocate(ACF(0:lmax,0:lmax),BCF(0:lmax,0:lmax))
    allocate(XT(nnod,3),NT(ntri,3))
    indx=0

end subroutine


subroutine generate_gsphere(rmax)
    ! Generate a sample Gaussian sphere with identification number
    ! nss, then move to discretize and output:
    real(kind=dp), intent(out) :: rmax
    integer :: j0

    do j0=1,nss
     if (gflg==1) then
      call SGSCF(ACF,BCF,SCFSTD,lmin,lmax)
     else
      call SGSAXCF(ACF,CEU,SEU,SCFSTD,lmin,lmax)
     endif
    end do

    call TRIDSA(MUT,PHIT,IT,ntr)
    if (gflg==1) then
        call RGSTD(XT,NT,MUT,PHIT,ACF,BCF,rmax,beta, &
        IT,nnod,ntri,lmin,lmax)

    else
        call RGSAXTD(XT,NT,MUT,PHIT,ACF,CEU,SEU,rmax,beta, &
                   IT,nnod,ntri,lmin,lmax)
    endif

end subroutine


function simplex3volume(u,v,b) result(A)
    real(kind=dp),intent(in) :: u(3),v(3),b(3)
    real(kind=dp) :: A
    A=(u(2)*v(3)-u(3)*v(2))*b(1)
    A=A+(u(3)*v(1)-u(1)*v(3))*b(2)
    A=A+(u(1)*v(2)-u(2)*v(1))*b(3)
    A=1.0_dp/6.0_dp*abs(A)
end function



subroutine compute_volume(volume)
    real(kind=dp), intent(out) :: volume
    real(kind=dp) :: O(3),B(3),C(3),D(3)
    integer :: j1
    O=(/0.0_dp,0.0_dp,0.0_dp/)
    volume=0.0_dp
    do j1=1,ntri
        B=XT(IT(j1,1),:)
        C=XT(IT(j1,2),:)
        D=XT(IT(j1,3),:)
        volume=volume+simplex3volume(B,D,C)
    end do
end subroutine


subroutine save_gsphere()
    integer :: j2,j1
    character(32) :: outfile2
    indx=indx+1

    outfile2=outfile
    write(outfile2,'(A,I1)') trim(adjustl(outfile)),indx

    open(unit=1, file=trim(outfile2) // '.vtk')               ! VTK
    write (1,150) '# vtk DataFile Version 2.0'
    write (1,150) 'gsphere output            '
    write (1,150) 'ASCII                     '
    write (1,150) 'DATASET POLYDATA          '
    write (1,160) 'POINTS ',nnod,' float'
150     format(a26)
160     format(a7,I7,A7)
    do j1=1,nnod
        write (1,*) (XT(j1,j2),j2=1,3)
    end do
    write (1,180) 'POLYGONS ',ntri,4*ntri
180     format(a9,I7,I7)
    do j1=1,ntri
        write (1,*) 3,(IT(j1,j2)-1,j2=1,3)
    end do
    close(unit=1)

    print *,'Wrote g-sphere in ' // trim(outfile2) // '.vtk'

end subroutine

end module 
