module chn2

use parameters
use voper

contains

! TODO give x,y,z as arguments
! TODO numfaces as an argument
! otherwise DO NOT TOUCH

subroutine convhull(numedges,edge,nedge1,nedge2,x,y,z,numfaces)

!use constants_m,only:x,y,z,numfaces
implicit real(kind=dp) (a-h,o-z)
real(kind=dp),dimension(3)::stick,tvec,base,vectn,crvec,cvec 
real(kind=dp), dimension(0:nmax),intent(inout) :: x,y,z
integer::numedges(nedge1),edge(nedge1,nedge2),listch(nedge1),inlist(nedge1)

!print *,'convhull starts'

! origin is the zeroth point (if not inside, use cmass)
x(0)=0.
y(0)=0.
z(0)=0.
! First find the largest z-value (any min,max,coord would do) to establish
! one vertex of the CH
zmax=z(1)
mind=1
!mind2=0
do i=2,numfaces
        if (z(i).gt.zmax) then
          	zmax=z(i)
          	mind=i
        end if
end do

!print *,'finding largest z value'

! listch is the list of all CH vertices; inlist is either 1 or 0 depending
! on whether the point is in listch
do i=1,numfaces
        inlist(i)=0
        numedges(i)=0
end do 
listch(1)=mind
inlist(mind)=1

!print *,'inlist altered'

! Then find the point with the largest 'umbrella angle' from the first 
! point (umbrella 'stick' towards origin) to establish one edge of the CH
call mvector(mind,0,stick,x,y,z)
sticklen=sqrt(dot(stick,stick))
dmin=1.
do i=1,numfaces
        if (i.ne.mind) then
        	call mvector(mind,i,tvec,x,y,z)
          	tveclen=sqrt(dot(tvec,tvec))
          	dtest=dot(stick,tvec)/(sticklen*tveclen)
          	if (dtest.lt.dmin) then
           		dmin=dtest
           		mind2=i
          	end if
        end if
end do

!print *,'found point',mind2,mind
!print *,sticklen

listch(2)=mind2
inlist(mind2)=1
edge(mind,1)=mind2
edge(mind2,1)=mind
numedges(mind)=1
numedges(mind2)=1

!print *,'initial 2 points found'

! nvert counts the vertices for which all edges have been  determined; 
! nconn counts the points that have been connected to any other point (and
! are thus vertices of CH)
nvert=0
nconn=2

! Begin nvert-loop
10 inda=listch(nvert+1)
indb=edge(inda,1)
call mvector(inda,indb,base,x,y,z)
! Begin numedges-loop: find a new edge for the vertex inda.
! Go through all points; i is the trial point, indn the reference point
20 do i=0,numfaces
        if ((i.ne.inda).and.(i.ne.indb)) then
          	call mvector(inda,i,cvec,x,y,z)
          	clen=sqrt(dot(cvec,cvec))
          	if ((i.eq.0).or.(dot(cvec,crvec).lt.(clen*crlen*(-tiny)))) then
! the trial point becomes the new reference point (.lt.(-tiny) in the above
! test for counterclockwise ordering of edges, .gt.tiny for clockwise 
! ordering)
           		indn=i
           		do j=1,3
            			vectn(j)=cvec(j)
           		end do
           		call crossproduct(vectn,base,crvec)
           		crlen=sqrt(dot(crvec,crvec))
          	end if
        end if
end do

!print *,'finding new edge'

! check if we've not yet gone round inda (not found all its edges)
if ((numedges(inda).eq.1).or.((indn.ne.edge(inda,1)).and.(indn.ne.edge(inda,2)))) then
! check if a new vertex of CH has been found
        if (inlist(indn).eq.0) then
          	nconn=nconn+1
          	listch(nconn)=indn
          	inlist(indn)=1
          	numedges(indn)=1
          	edge(indn,1)=inda
        end if
        numedges(inda)=numedges(inda)+1
        if(numedges(inda) .gt. numfaces) then
          print *,'TROUBLE IN CH: hull failed to construct.'
          stop
        end if
        edge(inda,numedges(inda))=indn
        indb=indn
        do j=1,3
          	base(j)=vectn(j)
        end do
        goto 20
end if
nvert=nvert+1
! check if we haven't yet found edges for all connected points (otherwise 
! we have determined the CH)

!print *,'1 iteration complete'
if (nvert.lt.nconn) goto 10

end

end module chn2






