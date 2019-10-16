module Minko

use parameters
use specfunc
use chn2
use lm
use ellips
use discrets
use voper
use chull

implicit none


contains


subroutine MINKO_GSPHERE(ACFIN,BCFIN,nnod,ntr,lmin,lmax,outfile,fflg,dflg,convflg)

implicit none

integer,intent(in) :: ntr,nnod,lmin,lmax,fflg,dflg,convflg
integer::IT(260000,3)
real(kind=dp),dimension(130000)::MUN,PHIN
real(kind=dp),intent(inout),allocatable :: ACFIN(:,:),BCFIN(:,:)
real(kind=dp), dimension(0:256,0:256) :: ACF,BCF


integer::i,j
real(kind=dp),dimension(nmax)::initarea
real(kind=dp)::nor(nmax,3)!,P(npar) !nor nmax 3, P npar
real(kind=dp) :: XN(130000,3),NT(260000,3)
real(kind=dp),dimension(260000)::AR,AT,MUT,PHIT
real(kind=dp)::AX(6),rmax,atot,coef,nut,DAR(7200,256)

 character(30)::outfile

integer::nacf,npar,ntri!,numedges(10),edge(10,10)!,ntr !,NCF(maxpar,2)
integer::j1,j2
real(kind=dp) :: norm1,norm2
integer :: numfaces

ntri=8*ntr**2
numfaces = ntri

! Convert dynamic coefficient arrays into static arrays

ACF = 0.
BCF = 0.

norm1=1.

do j1=0,lmax
  if(convflg .eq. 3) norm1=1.0d0/(2.0d0*pi)
  ACF(j1,0)=ACFIN(j1,0) !*norm1
end do

do j1=1,lmax
    if(convflg.eq.3) norm1=2.0d0/(2.0d0*pi)
    do j2=1,j1
        norm2=norm1
        ACF(j1,j2)=ACFIN(j1,j2) !*norm2
        BCF(j1,j2)=BCFIN(j1,j2) !*norm2
    end do
end do

! Revert original scaling; Minko normalizes the polynomials
! by itself, no need for scaling like this

if (convflg .eq. 3) then

do j=0,lmax

 if (j.eq.0) then
    do i=0,lmax
        ACF(i,0)=ACF(i,0)/sqrt((2*i+1)/2.0d0)
   end do
 else
    do i=j,lmax
       ACF(i,j)=ACF(i,j)/ &
       sqrt(2.0d0*FACTI(i-j)/FACTI(i+j))
       BCF(i,j)=BCF(i,j)/ &
       sqrt(2.0d0*FACTI(i-j)/FACTI(i+j))
   end do
endif

end do

end if

! Create a discretization for a unit sphere

call TRIDS(MUN,PHIN,IT,nnod,ntri,ntr)

nacf=size(ACFIN)
npar=4+(lmax+1)**2-lmin**2+5

! MUN,PHIN=node radius spherical coordinates
do i=1,6
        AX(i)=1.0d0 ! Dimensions for a sphere (sphere good)
end do

! Generate normals for the sphere

call RELLTD(XN,NT,MUN,PHIN,AX,rmax,IT,nnod,ntri)

! NT unit normals of each triangle, XN coordinates of triangle vertices,
! MUN, PHIN spherical coordinate angles,

! Generate areas scaled for a unit sphere

call ARELLTD(AT,MUN,PHIN,AX,rmax,IT,nnod,ntri)

do i=1,ntri
    MUT(i)=NT(i,3)
    nut=sqrt(1.0d0-MUT(i)**2)
    PHIT(i)=acos(NT(i,1)/nut)
    if (NT(i,2).lt.0.0d0) PHIT(i)=2.0d0*pi-PHIT(i)
    do j=1,3
        nor(i,j)=NT(i,j)
    end do
end do
! MUT,PHIT = facet normal spherical coordinates
           
! Use spherical harmonics to scale the log-areas

call ARDAR(AR,DAR,ACF,BCF,AT,MUT,PHIT,ntri,nacf,lmin,lmax)


! Scale the areas for a reasonably sized object   
    
atot =0.0d0
do j=1,ntri
    atot=atot+AR(j)
end do
coef=1000./atot
    
do j=1,ntri
    initarea(j)=coef*AR(j)
end do


 open(unit=1,file='acf.out')
 write(1,*) ACF(0:6,0:6)
 write(1,*) BCF(0:6,0:6)
! write(1,*) initarea(1:20)
 close(1)  

! Find the polyhedron corresponding to the 'EGI' given by the areas
! and the normals
call Minkowski(initarea,nor,ntri,numfaces,outfile,fflg,dflg)

end subroutine MINKO_GSPHERE




!----------------------------------------------------

! A program for determining a convex polyhedron from its Extended
! Gaussian Image by Minkowski. 
! Writes every 20th result to file minkotmp.

! Original author: Johanna Torppa, updated in 2017 for a more portable spec
! by Teo Korhonen.

!---------------------------------------------------

subroutine Minkowski(initarea,nor,imcmc,numfaces,outfile,fflg,dflg)

!use constants_m

!implicit none
integer,intent(in) :: fflg, dflg
integer::i,j,k,ncoef,posroute,numedges(nmax),r(nmax,nedmax),numvertices,maxind,first,iter,nwrite,&
	 imcmc,edge(nmax,nedmax),vert,n,l,temp,numvert(nmax)

integer :: ntr,triangles(nmax,3)
real(kind=dp),dimension(nmax)::initarea,d,pxk,area,xk,prearea,ek,newx,w
real(kind=dp), dimension(0:nmax) :: x,y,z
real(kind=dp),dimension(3)::anorm,bnorm,p1,p2,p0,cross,normal,rv,sv,addvect,ra,rb,cmass,avect,bvect,&
			      av,bv,cv,norm
real(kind=dp)::dhelp,initm,m,prevolume,eta,initareasum,initarealen,&
		 lower,rx(3,nmax),sqlength,maxpxk,h(1200),max,volume,scoef,retemp,premaxpxk,&
		 arealen,alphak,tk,coeff,f,wmin,nor(nmax,3)
integer :: numfaces
 character(30)::outfile,matlabx,matlaby,matlabz
integer :: polygonsandvertices

real(kind=dp) :: xmax,ymax,zmax

nwrite=0
iter=0

!write(ifile,'(I5)') imcmc
!ifile=adjustl(ifile)


! An auxiliary value initm (related to the maximum value of dot(a(i),b(j)),
! discarding values above 0.9999, corresponding to the normals of the i'th
! and j'th facet), is computed.

initm=0.
!write(2,*) numfaces
do i=1,numfaces
	do k=1,3
        	anorm(k)=nor(i,k)
        end do
        do j=1,numfaces
        	do k=1,3
        		bnorm(k)=nor(j,k)
         	end do
          	dhelp=1-(dot(anorm,bnorm))**2
          	if (dhelp.gt.0.0001) then
           		m=1./sqrt(dhelp)
           		if (m.gt.initm) initm=m
          	end if
        end do
end do
ncoef=0
posroute=1
prevolume=0.
eta=1.3
! The vector of the distances of the faces from the origin is initialized.
initareasum=0
do i=1,numfaces
	initareasum=initareasum+initarea(i)
end do

initarealen=sqrt(gdot(initarea,initarea,numfaces))
do i=1,numfaces
 	d(i)=initarealen**2/initareasum
end do
! The origin is x(0),y(0),z(0)
x(0)=0.
y(0)=0.
z(0)=0.

! The polytope in R^3 is transformed into a polytope in dual space, i.e., 
! the faces are transformed into points x(i),y(i),z(i).
10 do i=1,numfaces
	x(i)=nor(i,1)/d(i)
 	y(i)=nor(i,2)/d(i)
 	z(i)=nor(i,3)/d(i)
end do

!if(iter .eq. 0) write(2,*) numfaces, initarea

call convhull(numedges,edge,nmax,nedmax,x,y,z,numfaces)

! The convex hull of the dual polytope is transformed into a polytope in 
! R^3, i.e., the faces surrounding a vertex in the dual space become 
! vertices surrounding a face in R^3.

k=1
! k counts the total number of vertices, vert the vertices for face i
do i=1,numfaces
	vert=0
        do n=1,numedges(i)
          	if(numedges(i).gt.nedmax) then
			print*,'size of edge exceeded',nedmax
             		stop
          	endif
          	
          	call mvector(i,edge(i,n),p1,x,y,z)
          	call mvector(i,edge(i,modlo(n+1,numedges(i))),p2,x,y,z)
          	call mvector(0,i,p0,x,y,z)
          	call crossproduct(p1,p2,cross)
          	lower=dot(p0,cross)
! Coordinates of a vertex
		do l=1,3
           		if(k.gt.nmax) then
              			print*,'size of rx exceeded',nmax
              			stop
           		endif
           		rx(l,k)=cross(l)/lower
          	end do
! Check if this vertex has been computed earlier
          	do j=1,k-1
          		sqlength=0.
           		do l=1,3
            			sqlength=sqlength+(rx(l,j)-rx(l,k))**2
           		end do
! The test allows for some rounding errors etc.
           		if (sqlength.lt.tiny) then
            			temp=j
            			if (n.gt.1) then
! Ignore this vertex if met earlier at this face
             				if(vert.gt.nedmax) then
                				print*,'size of r exceeded',nedmax
                				stop
             				endif
             				if ((temp.eq.r(i,vert)).or.(temp.eq.r(i,1))) goto 70
            			end if
! If met earlier but not at this face, increase vert but not k
            			goto 60
           		end if
! End for j
          	end do
          	temp=k
          	k=k+1
!          	 if (k.gt.(nmax-10)) write(6,*) k
60        	vert=vert+1
! r holds the vertices of facet i
          	if(vert.gt.nedmax) then
            		print*,'size of r exceeded',nedmax
            		stop
          	endif
          	r(i,vert)=temp
! End for n
70      end do
	if(i.gt.nmax) then
            	print*,'size of numvert exceeded'
            	stop
        endif
        numvert(i)=vert
! End for i
end do

numvertices=k-1
! Now the points of the dual space can be erased, and x,y,z will represent
! coordinates in R^3.
do i=1,numvertices
	x(i)=rx(1,i)
        y(i)=rx(2,i)
        z(i)=rx(3,i)
end do
maxpxk=0.
do i=1,numfaces
        h(i)=0.
! Compute 'vertices' and auxiliary values h for faces that have no vertices
        if (numvert(i).eq.0) then
          	do j=1,3
           		normal(j)=nor(i,j)
          	end do
          	max=0.
          	do l=1,numvertices
           		call mvector(0,l,rv,x,y,z)
           		if (dot(rv,normal).gt.max) then
            			max=dot(rv,normal)
            			maxind=l
           		end if
          	end do
          	numvert(i)=1
          	r(i,1)=maxind
          	do l=1,numvertices
           		call mvector(0,l,rv,x,y,z)
           		if ((dot(rv,normal).eq.max).and.(l.ne.maxind)) then
            			numvert(i)=2
            			r(i,2)=l
           		end if
          	end do
          	h(i)=d(i)-max
         end if
! The circuits of the edges of each face are computed and the largest 
! circuit is seeked.
         pxk(i)=0
         do j=1,numvert(i)
          	call mvector(r(i,j),r(i,modlo((j+1),numvert(i))),sv,x,y,z)
          	pxk(i)=pxk(i)+sqrt(dot(sv,sv))
         end do
         if (maxpxk.lt.pxk(i)) maxpxk=pxk(i)
end do

! Some auxiliary values are computed
volume=0.
do i=1,numfaces
        do j=1,3
          	addvect(j)=0.
        end do
        do n=1,numvert(i)
          	call mvector(0,r(i,n),ra,x,y,z)
          	call mvector(0,r(i,modlo(n+1,numvert(i))),rb,x,y,z)
          	call crossproduct(ra,rb,cross)
          	do j=1,3
           		addvect(j)=addvect(j)+cross(j)
          	end do
        end do
        area(i)=sqrt(dot(addvect,addvect))/2.
        volume=volume+d(i)*area(i)/3.
end do
iter=iter+1
if (volume.ge.prevolume) then
        do j=1,3
          	cmass(j)=0.
        end do
        do i=1,numfaces
          	do n=1,numvert(i)-2
           		call mvector(r(i,1),r(i,n+1),avect,x,y,z)
           		call mvector(r(i,1),r(i,n+2),bvect,x,y,z)
           		call crossproduct(avect,bvect,cross)
           		call mvector(0,r(i,1),av,x,y,z)
           		call mvector(0,r(i,n+1),bv,x,y,z)
           		call mvector(0,r(i,n+2),cv,x,y,z)
           		do j=1,3
            			cmass(j)=cmass(j)+d(i)*sqrt(dot(cross,cross))*(av(j)+bv(j)+cv(j))
           		end do
          	end do
        end do
        do j=1,3
          	cmass(j)=cmass(j)/(24*volume)
        end do
        do i=1,numfaces
          	do j=1,3
           		norm(j)=nor(i,j)
          	end do
          	xk(i)=d(i)-dot(norm,cmass)-h(i)
        end do
        scoef=initarealen**2/gdot(xk,initarea,numfaces)
        do i=1,numfaces
          	xk(i)=scoef*xk(i)
          	area(i)=scoef**2*area(i)
        end do
        maxpxk=scoef*maxpxk
        volume=scoef**3*volume
        if (prevolume.ne.0) then
          	retemp=eta**ncoef*(arealen*sin(alphak))**2/(6*initm*premaxpxk)
        else
          	retemp=0.
        end if
        if (((volume-prevolume).ge.retemp).and.(posroute.eq.1)) then
          	ncoef=ncoef+1
        else
          	ncoef=ncoef-1
        end if
! else, i.e., if volume.lt.prevolume
else
        ncoef=ncoef-1
! No changes are made to xk:s.
        volume=prevolume
        maxpxk=premaxpxk
        do i=1,numfaces
          	area(i)=prearea(i)
        end do
end if

!print *,numvertices

if((volume-prevolume.lt.1.0d-12.and.volume-prevolume.gt.1.0d-15)) goto 100
prevolume=volume
premaxpxk=maxpxk
do i=1,numfaces
        prearea(i)=area(i)
end do

! The iteration step and the angle alphak are computed.
arealen=sqrt(gdot(area,area,numfaces))
alphak=acos(gdot(area,initarea,numfaces)/(arealen*initarealen))

! Test whether to stop iterating
if (alphak.lt.epsilon0) goto 100
! Writes the temporary results into the temporaryfile 'minkotmp'

if(nwrite*20+1.eq.iter) then
         nwrite=nwrite+1         
         open(3,file='minkotmp')
         write(3,*) numvertices, numfaces
         scoef=scoef*sqrt(initarealen/arealen)
         do i=1,numvertices
           	x(i)=x(i)-cmass(1)
           	y(i)=y(i)-cmass(2)
           	z(i)=z(i)-cmass(3)
           	write(3,110) scoef*x(i),scoef*y(i),scoef*z(i)
         end do
         110 format(3(E23.16,1X))
         do i=1,numfaces
            	write(3,*) numvert(i)
            	write(3,*) (r(i,j),j=1,numvert(i))
         end do
         close(3)
endif
tk=eta**ncoef*arealen*sin(alphak)/(3*initm*maxpxk)
if(tk.lt.1.d-8) goto 100
coeff=gdot(area,initarea,numfaces)/(initarealen**2)
do i=1,numfaces
        f=area(i)-coeff*initarea(i)
        ek(i)=f/(arealen*sin(alphak))
end do
do i=1,numfaces
        newx(i)=xk(i)+tk*ek(i)
        if (newx(i).le.0) goto 90
end do
do i=1,numfaces
        d(i)=newx(i)
end do
posroute=1
goto 10
90 posroute=-1
first=0
do i=1,numfaces
        w(i)=xk(i)/(coeff*initarea(i)-area(i))
        if ((w(i).gt.0).and.(first.eq.0)) then
          	wmin=w(i)
          	first=1
        end if
        if ((w(i).gt.0).and.(w(i).lt.wmin)) wmin=w(i)
end do
tk=0.9*arealen*wmin*sin(alphak)
do i=1,numfaces
        d(i)=xk(i)+tk*ek(i)
end do
goto 10

! Printing the output:


100 scoef=0 ! scoef*sqrt(initarealen/arealen)

! The final vertices.

do i=1,numfaces
  if(d(i) .gt. scoef) scoef=d(i)
end do

do i=1,numvertices
        x(i)=(x(i)-cmass(1))/scoef
        y(i)=(y(i)-cmass(2))/scoef
        z(i)=(z(i)-cmass(3))/scoef
end do

! Output as triangles instead of arbitrary polygons
ntr = 0
if(dflg.eq.2) then
    call triangulate(r,triangles,numvert,numfaces,ntr)
end if


! NOTE: Matlab format only prints the vertices.

if(fflg.eq.1) then
   matlabx= trim(outfile) // 'x.out'
   matlaby= trim(outfile) // 'y.out'
   matlabz= trim(outfile) // 'z.out' 
   open(unit=1, file=matlabx)           ! Matlab
   open(unit=2, file=matlaby)
   open(unit=3, file=matlabz)

! do i=1,numfaces
!	write(1,125) x(j),j=1,numvert(i)
!	write(2,125) y(j),j=1,numvert(i)
!	write(3,125) z(j),j=1,numvert(i)
! end do

  do i=1,numvertices
   write (1,*) x(i)
   write (2,*) y(i)
   write (3,*) z(i)
  end do
   close(unit=3)
   close(unit=2)
   close(unit=1)

   print *,'Wrote x-coordinates in ' // trim(outfile) // 'x.out'
   print *,'Wrote y-coordinates in ' // trim(outfile) // 'y.out'
   print *,'Wrote z-coordinates in ' // trim(outfile) // 'z.out'

elseif(fflg.eq.3) then
   open(unit=1, file=trim(outfile) // '.idf')               ! IDL
   write (1,*) numvertices,numfaces
   do i=1,numvertices
    write (1,*) x(i),y(i),z(i)
   end do
   
   if(dflg.eq.2) then
    do i=1,ntr
       write(1,*) 3
       write(1,*) (triangles(i,j),j=1,3)
    end do
   else
    do i=1,numfaces
     write (1,*) numvert(i)
     write (1,*) (r(i,j),j=1,numvert(i))
    end do
   end if
   print *,'Wrote g-sphere in ' // trim(outfile) // '.idf'

   close(unit=1)

elseif(fflg.eq.2) then

   open(unit=2,file=trim(outfile)//'.vtk')

   polygonsandvertices = 0
   do i=1,numfaces
    polygonsandvertices = polygonsandvertices + numvert(i) + 1
   end do

   write (2,150) '# vtk DataFile Version 2.0'   
   write (2,150) 'gsphere output            '
   write (2,150) 'ASCII                     '
   write (2,150) 'DATASET POLYDATA          '
   write (2,160) 'POINTS ',numvertices,' float'
150     format(a26)
160     format(a7,I7,A7)
   do i=1,numvertices
    write (2,*) x(i),y(i),z(i)
   end do
   

   if(dflg.eq.2) then
    write (2,180) 'POLYGONS ',ntr,4*ntr
    do i=1,ntr
       write(2,*) 3,(triangles(i,j)-1,j=1,3)
    end do
    
   else
    write (2,180) 'POLYGONS ',numfaces,polygonsandvertices
180     format(a9,I7,I7)
    do i=1,numfaces
         write(2,*) numvert(i),(r(i,j)-1,j=1,numvert(i))
    end do
   end if
   close(2)

   print *,'Wrote g-sphere in ' // trim(outfile) // '.vtk'

end if

! Calculates some stats.
open(unit=3,file=trim(outfile) // 'stat.out')

xmax=0.
ymax=0.
zmax=0.

do i=1,numvertices
  if(x(i) .gt. xmax) xmax = x(i)
  if(y(i) .gt. ymax) ymax = y(i)
  if(z(i) .gt. zmax) zmax = z(i)
end do

write(3,*) xmax,ymax,zmax

close(3)

end subroutine


! Uses fan triangulation to convert an arbitrary polytope into a triangle
! polygonal model. Teo Korhonen, 2017.

subroutine triangulate(facets,triangles,numvert,numfaces,ntr)

implicit none
integer ::  facets(nmax,nedmax),triangles(nmax,3),numvert(nmax),numfaces
integer :: i,j,ntr

! ntr counts the number of triangles.

ntr = 0

open(unit=4,file='trilog.out')

do i=1,numfaces
    
    ! If facet has 1 or 2 points, no need to triangulate.
    if (numvert(i) .lt. 3) then
        write(4,*) 'no additions'
        cycle
        
    ! If facet is already a triangle, add it as such.
    elseif (numvert(i) .eq. 3) then
        ntr = ntr+1
        triangles(ntr,1:3) = facets(i,1:3)
        write(4,*) ntr
        write(4,*) 'added', triangles(ntr,1:3)
        cycle 
    end if
    
    ! If not, turn it into triangles using fan triangulation.    
    do j=3,numvert(i)
        ntr = ntr+1        
        triangles(ntr,1) = facets(i,1)
        triangles(ntr,2) = facets(i,j-1)
        triangles(ntr,3) = facets(i,j)   
        
        write(4,*) ntr
        write(4,*) 'added',triangles(ntr,1:3)     
    end do    
end do

close(4)

end subroutine


! The dot product in gradient space. Johanna Torppa.

function gdot(avect,bvect,numfaces)
!use constants_m,only:rk,nmax,numfaces
implicit none
integer::i,numfaces
real(kind=dp)::avect(nmax),bvect(nmax),gdot

gdot=0.
do i=1,numfaces
        gdot=gdot+avect(i)*bvect(i)
end do
end

end module minko

!----------------------------------------------

!include 'Lib/chn2.f90'
!include 'Lib/vectprocs.f90'
!include 'Lib/Ellips.f90'
!include 'Lib/Voper.f90'
!include 'Lib/LM.f90'
!include 'Lib/Lditd.f90'
!include 'Lib/Scatlaw.f90'  ! Scattering law.
!include 'Lib/Discrets.f90' ! Particle discretization.          
!include 'Lib/Randev.f90'   ! Random deviates.                  
!!include 'Lib/Pario.f90'    ! Parameter input/output.
!include 'Lib/Viewgeo.f90'  ! Viewing geometry.
!include 'Lib/HG1G2.f90'    ! H, G1, G2 phase function.
!include 'Lib/Insectcx.f90' ! Particle-line intersections.     
!include 'Lib/Specfunc.f90' ! Special functions.  


