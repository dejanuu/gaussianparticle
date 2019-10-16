program GSPHERE

! GSPHERE generates sample Gaussian spheres.
! Version 2017-05-19, based on the version 2003-12-10.
!
!
! Free software license 
!
! The code is available under the GNU General Public License.
!
! Contact addresses for authors:
!
! Karri Muinonen
! E-mail: Karri.Muinonen@helsinki.fi
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Publi! License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!
! This version of GSPHERE has been updated and refactored in 2017 to 
! modern Fortran standard from the original F77 version by Teo Korhonen.
!
!
! Required modules:
!
use corrfunc     ! Correlation functions.
use discrets     ! Discretization.
use gsaxg        ! Axisymmetric G-sphere generator.
use gsg          ! G-sphere generator.
use minko       ! Solves Minkowski problem, used for convex spheres.
use chull       ! Finds the convex hull of a sphere using the 3-D gift wrapping algorithm.
       
implicit none
       
       
integer :: nss,lmin,lmax, &
j0,j1,j2,gflg,dflg,cflg,fflg,convflg,seed(12),seedsize=1

integer :: nthe,nphi,ntr,nnod,ntri
integer,allocatable :: IT(:,:)
real(kind=dp) :: XS(0:180,0:360,3), &  
 MUS(0:180),PHIS(0:360)
real(kind=dp),allocatable :: MUT(:),PHIT(:),XT(:,:),NT(:,:)
real(kind=dp), dimension(100) :: P

real(kind=dp) :: CEU(3),SEU(3),a,sig,beta,gami,elli,nuc,rmax,rd
real(kind=dp),allocatable :: CSCF(:),SCFSTD(:,:),ACF(:,:),BCF(:,:) !,testminkonor(:,:),testminkoar(:)
                                                                   ! test Minkowski with a unit cube
integer :: j,k,l,lb,nacf=30

character :: ca*2
character(32) :: infile,outfile,matlabx,matlaby,matlabz, &
                 letterx='x',lettery='y',letterz='z'

logical :: there
integer :: irnd
common irnd

! Initializations:

rd=pi/180.0d0

irnd=5
       
! Input file specified in command line argument:

call getarg(1, infile)
WRITE(*,*) TRIM(infile)
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
read (1,20) seed(1)  ! Seed for the random number generator.
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


! Initialize random number generator

do j0=2,12
    seed(j0)=j0*seed(1)
end do

call random_seed(size=seedsize)
call random_seed(put=seed)


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

! Generate a sample Gaussian sphere with identification number
! nss, then move to discretize and output:

do j0=1,nss
 if (gflg.eq.1) then
  allocate(ACF(0:lmax,0:lmax),BCF(0:lmax,0:lmax))
  call SGSCF(ACF,BCF,SCFSTD,lmin,lmax)
 else
  allocate(ACF(0:lmax,0:lmax))
  call SGSAXCF(ACF,CEU,SEU,SCFSTD,lmin,lmax)
 endif
end do


! Output a convex sphere (EXPERIMENTAL)
if(convflg.eq.3) then
 ! Add file flag!
 call MINKO_GSPHERE(ACF,BCF,nnod,ntr,lmin,lmax,outfile,fflg,dflg,3)
 
 deallocate(ACF,BCF,SCFSTD,CSCF)

! Tests the Minkowski solver with given spherical harmonics

elseif(convflg.eq.4) then
 ! Add file flag!
 ACF = 0.
 BCF = 0.
 open(unit=1,file='mcmc1.out')
 
 nacf = 9+(lmax+1)*lmax/2
 read(1,*) (P(j),j=1,9+(lmax+1)**2)
 l=9
 lb=nacf+1
 ACF=0.0
 do j=0,lmax
	do k=0,j
		l=l+1
		ACF(j,k)=P(l)
		if (k>0) then
			lb=lb+1
			BCF(j,k)=P(lb)
		endif
	end do
 end do
 
 close(1)
 
 call MINKO_GSPHERE(ACF,BCF,nnod,ntr,lmin,lmax,outfile,fflg,dflg,4)
 
 deallocate(ACF,BCF,SCFSTD,CSCF)
 
! Spherical harmonic coefficient output


elseif (dflg.eq.3) then
    open(unit=1,file=trim(outfile))
    do j1=lmin,lmax
        write(1,112) (ACF(j1,j2),j2=0,j1)
    end do
    
    do j1=lmin,lmax
        write(1,112) (BCF(j1,j2),j2=0,j1)
    end do
    
112 format(256(E14.8))

    close(unit=1)

elseif (dflg.eq.1) then

! Spherical-coordinate representation for general and axisymmetric shapes:

  call SPHDS(MUS,PHIS,nthe,nphi)
  if (gflg.eq.1) then
   call RGSSD(XS,MUS,PHIS,ACF,BCF,rmax,beta, &
             nthe,nphi,lmin,lmax)
  else
   call RGSAXSD(XS,MUS,PHIS,ACF,CEU,SEU,rmax,beta, &
               nthe,nphi,lmin,lmax)
  endif

  open(unit=1, file=outfile // 'x.out')            ! Matlab
  open(unit=2, file=outfile // 'y.out')
  open(unit=3, file=outfile // 'z.out')
  do j1=0,nthe
   write (1,115) (XS(j1,j2,1),j2=0,nphi)
   write (2,115) (XS(j1,j2,2),j2=0,nphi)
   write (3,115) (XS(j1,j2,3),j2=0,nphi)
  end do
115     format(500(E14.8,1X))
 close(unit=3)
 close(unit=2)
 close(unit=1)

else

! Triangle representation for general and axisymmetric shapes:
        
  allocate(IT(8*ntr**2,3),PHIT(4*ntr**2+2),MUT(4*ntr**2+2))
  
  call TRIDSA(MUT,PHIT,IT,nnod,ntri,ntr)
  if (gflg.eq.1) then
   call RGSTD(XT,NT,MUT,PHIT,ACF,BCF,rmax,beta, &
    IT,nnod,ntri,lmin,lmax)
   deallocate(ACF,BCF)
  
  else
   call RGSAXTD(XT,NT,MUT,PHIT,ACF,CEU,SEU,rmax,beta, &
               IT,nnod,ntri,lmin,lmax)
   deallocate(ACF)
  endif
  
! Output the convex hull of the sphere
  
  if(convflg.eq.2) then
    call GSPGW(XT,nnod,outfile,fflg)
  end if


  if(fflg.eq.1) then
   matlabx= trim(outfile) // trim(letterx) // '.out'
   matlaby= trim(outfile) // trim(lettery) // '.out'
   matlabz= trim(outfile) // trim(letterz) // '.out' 
   open(unit=1, file=matlabx)           ! Matlab
   open(unit=2, file=matlaby)
   open(unit=3, file=matlabz)
   do j1=1,nnod
    write (1,*) XT(j1,1)
    write (2,*) XT(j1,2)
    write (3,*) XT(j1,3)
   end do
125     format(130000(E14.8,1X))
   close(unit=3)
   close(unit=2)
   close(unit=1)

   print *,'Wrote x-coordinates in ' // trim(outfile) // trim(letterx) // '.out'
   print *,'Wrote y-coordinates in ' // trim(outfile) // trim(letterx) // '.out'
   print *,'Wrote z-coordinates in ' // trim(outfile) // trim(letterx) // '.out'

  elseif(fflg.eq.3) then
   open(unit=1, file=trim(outfile) // '.idf')               ! IDL
   write (1,*) nnod,ntri
   do j1=1,nnod
    write (1,*) (XT(j1,j2),j2=1,3)
   end do
   do j1=1,ntri
    write (1,*) 3
    write (1,*) (IT(j1,j2),j2=1,3)
   end do
  
   print *,'Wrote g-sphere in ' // trim(outfile) // '.idf'

   close(unit=1)

  elseif(fflg.eq.2) then
   open(unit=1, file=trim(outfile) // '.vtk')               ! VTK
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
   
   print *,'Wrote g-sphere in ' // trim(outfile) // '.vtk'

  endif
  
 deallocate(IT,MUT,PHIT)

endif
end program GSPHERE
