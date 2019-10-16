module chull

use parameters
use vectors


contains

! An interface for the below algorithm for GSPHERE.
! Uses column-major arrays instead of row major.
! (main reason for this interface to exist)
! 
! Also prints the results in a .vtk file.

subroutine GSPGW(XT,nnod,outfile,fflg)

implicit none

real(kind=dp),intent(in),allocatable :: XT(:,:)
integer,intent(in) :: nnod,fflg
character(32) :: outfile,matlabx,matlaby,matlabz

real(kind=dp),allocatable :: x(:),y(:),z(:)
integer,allocatable :: edges(:,:),facets(:,:)
integer :: j1,j2

allocate(x(nnod),y(nnod),z(nnod))
allocate(edges(2,0:nnod**2),facets(3,0:nnod**2))

! Assign XT-node coordinates to x,y,z.

do j1=1,nnod
    x(j1) = XT(j1,1)
    y(j1) = XT(j1,2)
    z(j1) = XT(j1,3)
end do

call giftwrap(x,y,z,nnod,edges,facets)

! Print output:

if (fflg .eq. 2) then

   open(unit=1, file=trim(outfile) // 'convex' // '.vtk')               ! VTK
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
   write (1,180) 'POLYGONS ',facets(1,0),4*facets(1,0)
180     format(a9,I7,I7)
   do j1=1,facets(1,0)
    write (1,*) 3,(facets(j2,j1)-1,j2=1,3)
   end do
   close(unit=1)

print *,'Wrote convex hull to '//trim(outfile) // 'convex.vtk'

elseif(fflg.eq.3) then

   open(unit=1, file=trim(outfile) // '.idf')               ! IDL
   write (1,*) nnod,facets(1,0)
   do j1=1,nnod
    write (1,*) (XT(j1,j2),j2=1,3)
   end do
   do j1=1,facets(1,0)
    write (1,*) 3
    write (1,*) (facets(j2,j1),j2=1,3)
   end do
  
   print *,'Wrote convex hull to '//trim(outfile) //'.idf'

   close(unit=1)

elseif(fflg.eq.1) then
   matlabx= trim(outfile) // 'xconvex.out'
   matlaby= trim(outfile) // 'yconvex.out'
   matlabz= trim(outfile) // 'zconvex.out' 
   open(unit=1, file=matlabx)           ! Matlab
   open(unit=2, file=matlaby)
   open(unit=3, file=matlabz)
   do j2=1,3
    write (1,125) (XT(facets(j2,j1),1),j1=1,facets(1,0))
    write (2,125) (XT(facets(j2,j1),2),j1=1,facets(1,0))
    write (3,125) (XT(facets(j2,j1),3),j1=1,facets(1,0))
   end do
125     format(130000(E14.8,1X))
   close(unit=3)
   close(unit=2)
   close(unit=1)

   print *,'Wrote convex hull x-coordinates in ' // matlabx
   print *,'Wrote convex hull y-coordinates in ' // matlaby
   print *,'Wrote convex hull z-coordinates in ' // matlabz

end if
end subroutine


! Computes the edges of the convex hull of x,y,z using the 3D gift wrapping algorithm.
! Runs in O(nF) where n is # of points, F is # of facets. Worst-case scenario 
! should be roughly O(n^2)
!
! Teo Korhonen, 2017 
!
! Input:
!
! real arrays x,y,z         coordinates of each point
! integer array edges       edges of convex hull
! integer n                 number of points
!

subroutine giftwrap(x,y,z,n,edges,facets)

implicit none

real(kind=dp),allocatable :: x(:),y(:),z(:) ! Store ALL points
real(kind=dp) :: tempr1, epsilon0 = 0.00000001

type(vector),allocatable :: points(:)
type(vector) :: tempv1, tempv2, tempv3, tempv4, &
                ez=vector(0.,0.,1.)

!integer :: point1,point2,point3 ! Temporary point storage


integer :: n,i0,i1,i2, & ! number of elements
        !zmax,zmin,xmax,xmin,ymax,ymin,& ! indices of extreme points
        init1,init2,&         ! indices of the initial points
        newpoint, currentedge1,currentedge2,	 &  ! index of a new point
        indexofedge, &                              ! index of currently processed edge.
        whichside,whichsidenow !, &
!        nfacets  ! Counts the number of facets in the hull.
integer,allocatable :: edges(:,:),unprocessededges(:)
       ! 1st dim is 2, 2nd dim stores the indices of the endpoints of each edge
       ! Edges stores the indices of the edges of each facet. unprocessededges stores
       ! each edge that is yet to be processed.
integer,allocatable :: facets(:,:)

logical :: added1=.false.,added2=.false.,containsallpoints=.true., &
            facetcontainspoints = .false.


allocate(points(n))
allocate(unprocessededges(0:n**2))

edges(1:2,0:n**2) = 0
unprocessededges = 0
! The 0th element in unprocessededges stores the number

open(unit=1,file='log.out')
! Write a log.

! Find two initial points are found for an initial edge.
! First one has the biggest z (bound to be in the convex hull):

init1=1

do i1=1,n

    if (z(i1) .gt. z(init1)) then
        init1=i1
    end if
    
    ! Assign the point into vector format for simpler operations afterwards
    points(i1) = vector(x=x(i1),y=y(i1),z=z(i1))
    
end do

write(1,*)  'largest z',z(init1),'with point',init1

! Second one is found by gift wrapping in the x,z plane:

tempr1 = 0.
do i1=1,n

    if (i1 .eq. init1) cycle
    
    tempv1 = points(init1) - points(i1)
    tempv1%y = 0.
    
    if(angle2(tempv1,ez) .gt. tempr1) then
        tempr1 = angle2(tempv1,ez)
        init2 = i1
    end if

end do

! Create the first edge between the two points.

edges(1,1) = init1; edges(2,1) = init2
edges(1,0) = 1

write(1,*) 'initial edge between points',init1,init2

! edges stores the endpoints of all the edges in the hull.
! edges(1,0) is the number of edges in the hull.

unprocessededges(0) = 1
unprocessededges(1) = 1

! unprocessededges stores the indices of the edges that haven't been checked for more connections.
! Each element in unprocessededges is turned into 0 after that edge has been processed.

! Subsequent points are found by picking one edge,
! picking one additional point, making sure that the point is not already connected to the edge, 
! and ensuring that it is a part of the convex hull (all of the other points are on one side of 
! the triangle formed by the edge and the point).
!
! When no more connections can be found, the algorithm is complete.

! BIG LOOP



 currentedge1 = 1
 currentedge2 = 2
 indexofedge = 1
 
 


write(1,*) n,'points. Begin convex sphere determination.'

facets(1,0)=0

do i0=1,n**2 ! At maximum, there can be n**2 runs
    
    write(1,*)edges(1,0),'edges in the hull'
    write(1,*)facets(1,0),'facets in the hull'
    
    ! If no more edges left to process, the hull is ready.
    
    write(1,*) 'unprocessed edges left:',unprocessededges(0)
    
    if (unprocessededges(0) .eq. 0) exit

    ! Pick an unprocessed edge:
    
    do i1 = 1,edges(1,0)
        if (unprocessededges(i1) .gt. 0) then
            currentedge1 = edges(1,unprocessededges(i1))
            currentedge2 = edges(2,unprocessededges(i1))
            !write(1,*)currentedge1,currentedge2
            indexofedge = i1
            exit
        end if
    end do
    
    ! Create a temporary vector between the points of the edge (useful later for calculations)
    
    tempv1 = points(currentedge1) - points(currentedge2)
    
    
    ! Loop through all points. Find one that can be added to the hull, or that can be connected to one of the
    ! endpoints of the current edge. 
    
    newpoint = 0
    
    do i1=1,n
        
        write(1,*) 'checking point',i1,'of',n
        
        added1 = .false.
        added2 = .false.
        
        if (i1 .eq. currentedge1 .or. i1.eq. currentedge2) then
            write(1,*) 'point already in edge'
            cycle
        end if
        !write(1,*)currentedge1,currentedge2,i1
        
        ! Loop through all the edges in the hull, see if the new point is connected to the
        ! current edge endpoints via an existing edge.
        
        do i2=1,edges(1,0)
        
            if ((i1 .eq. edges(1,i2) .and. currentedge1 .eq. edges(2,i2)) &
                  .or. (currentedge1 .eq. edges(1,i2) .and. i1 .eq. edges(2,i2))) then
                
                write(1,*) 'points ',currentedge1, i1, 'have already been connected'
                added1 = .true.
                
            end if
            
            if ((i1 .eq. edges(1,i2) .and. currentedge2.eq.edges(2,i2)) &
                  .or. (currentedge2 .eq. edges(1,i2) .and. i1 .eq. edges(2,i2))) then
                write(1,*) 'points ',currentedge2, i1, 'have already been connected'
                added2 = .true.
            end if
            
        end do
        
        ! If both edge endpoints are connected to the point with an existing edge, 
        ! check for a possible new facet, then move on to the next point.
        
        if (added1 .and. added2) then
            
            ! Loop through each facet, see if these 3 points are
            ! connected through some facet.    
            
            do i2=1,facets(1,0)
                
                facetcontainspoints = .false.
            
                if (facets(1,i2).eq.currentedge1.and.facets(2,i2).eq.currentedge2 &
                   .and.facets(3,i2).eq.i1) then
                   facetcontainspoints = .true.
                   exit
                   
                elseif (facets(1,i2).eq.currentedge1.and.facets(2,i2).eq.i1 &
                   .and.facets(3,i2).eq.currentedge2) then
                   facetcontainspoints = .true.
                   exit
                   
                elseif (facets(1,i2).eq.currentedge2.and.facets(2,i2).eq.i1 &
                   .and.facets(3,i2).eq.currentedge1) then
                   facetcontainspoints = .true.
                   exit
                   
                elseif (facets(1,i2).eq.currentedge2.and.facets(2,i2).eq.currentedge1 &
                   .and.facets(3,i2).eq.i1) then
                   facetcontainspoints = .true.
                   exit
                   
                elseif (facets(1,i2).eq.i1.and.facets(2,i2).eq.currentedge1 &
                   .and.facets(3,i2).eq.currentedge2) then
                   facetcontainspoints = .true.
                   exit
                   
                elseif (facets(1,i2).eq.i1.and.facets(2,i2).eq.currentedge2 &
                   .and.facets(3,i2).eq.currentedge1) then
                   facetcontainspoints = .true.
                   exit
                   
                end if
            end do
            
            ! If this is not the case, then add the new facet.
            
            if (.not. facetcontainspoints) then
                facets(1,0) = facets(1,0) + 1
                facets(1,facets(1,0)) = currentedge1
                facets(2,facets(1,0)) = currentedge2
                facets(3,facets(1,0)) = i1
                write(1,*) 'new facet:',currentedge1,currentedge2,i1
            end if
            
            cycle
        end if
        
        ! Calculate a temporary vector between point i1 and c.edge1
        ! and a normal vector for the triangle (c.edge1,c.edge2,i1)
        
        tempv2 = points(currentedge1) - points(i1)
        tempv3 = tempv1 .cross. tempv2
        
        ! Check if all other points are on one side of the triangle (c.edge1,c.edge2,i1)
        
        whichside = -2
        whichsidenow = -2
        
        do i2=1,n
        
            if (i2 .eq. currentedge1 .or. i2 .eq. currentedge2 .or. i2 .eq. i1) cycle
            
            ! This checks which side of the triangle the current point is
            
            tempv4 = points(currentedge1) - points(i2)
            tempr1 = tempv3 .dot. tempv4 
            
            if(tempr1 .gt. epsilon0) then
                whichsidenow = 1
            else if(tempr1 .lt. epsilon0) then
                whichsidenow = -1
            else
                whichsidenow = 0
            end if
            
            if (whichside .eq. -2 .or. whichside .eq. 0) then
                whichside = whichsidenow
            end if
            
            if ((whichside .eq. 1 .and. whichsidenow .eq. -1) .or. &
                 (whichside .eq. -1 .and. whichsidenow .eq. 1)) then
                 containsallpoints = .false.
                 exit
            end if
            
        end do
        
        ! If all other points are not on one side of the triangle,
        ! don't add the new point.
    
        if (.not. containsallpoints) then
            write(1,*) 'All points not on one side of triangle'
            containsallpoints = .true.
            cycle
        end if
        
        ! If all these tests were passed, select the new point.
        
        newpoint = i1
        write(1,*)'found new point',newpoint,'we will now connect it to points',currentedge1,currentedge2
        exit

    end do
    
    
    ! Adds a new facet to the hull
    ! TODO: Facets are not being counted properly. Find the conditions of a
    ! legitimate new face, or how to construct the faces from the hull.
    
    if (.not.(added1 .and. added2) .and. newpoint .ne. 0) then
        facets(1,0) = facets(1,0) + 1
        facets(1,facets(1,0)) = currentedge1
        facets(2,facets(1,0)) = currentedge2
        facets(3,facets(1,0)) = newpoint
        write(1,*)'new facet:',currentedge1,currentedge2,newpoint
    end if
    
    
    
    if (newpoint.eq.0) then
        ! Attempt to find a new facet, even if no new points were found:
    
        write(1,*) 'no new points found for edge',currentedge1,currentedge2
        unprocessededges(indexofedge) = 0
        unprocessededges(0) = unprocessededges(0)-1
        cycle
    end if
    
    ! If the current point is not already connected to the endpoint 1:
    
    if (.not. added1) then
        
        ! Add the edge between endpoint 1 and the new point to edges:
        
        edges(1,0) = edges(1,0) + 1
        edges(1,edges(1,0)) = currentedge1; edges(2,edges(1,0)) = newpoint
        
        unprocessededges(0) = unprocessededges(0) + 1
        unprocessededges(edges(1,0)) = edges(1,0)
        write(1,*)newpoint,currentedge1,'connected'
        
    end if
    
    ! If the current point is already connected to the endpoint 2,
    ! add that edge to the hull and mark it as unprocessed.
    
    if (.not. added2) then
    
        ! Add the edge between endpoint 1 and the new point to edges:
        
        edges(1,0) = edges(1,0) + 1
        edges(1,edges(1,0)) = currentedge2; edges(2,edges(1,0)) = newpoint
        
        ! Flag the new edge as unprocessed:
        
        unprocessededges(0) = unprocessededges(0) + 1
        unprocessededges(edges(1,0)) = edges(1,0)
        write(1,*)newpoint,currentedge2,'connected'
        
    end if
    
    ! Now, the current edge has been processed.
    
    !if(i0 .gt. 1) then
    unprocessededges(indexofedge) = 0
    unprocessededges(0) = unprocessededges(0) - 1
    !endif
end do

! END BIG LOOP


close(1)

deallocate(unprocessededges)
deallocate(points)

end subroutine giftwrap


end module
