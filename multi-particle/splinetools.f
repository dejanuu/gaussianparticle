!----------------------------------------------------
!Tools for handling splines
!Copyright (C) 2016 Karri Muinonen, Timo Väisänen, and University of Helsinki
!All rights reserved.
!The new BSD License is applied to this software, see LICENSE.txt
module splinetools
    use sirisconstants
    use mathroutines
    implicit none
    !private cubicspline,bisectionSearch
contains

    !---------------------------------------------
    !initPresentation(xp,coeffs,x,y)
    !Fills xp and coeffs usign data obtained from x and y
    !
    !IN:    x: data for x-axis
    !       y(a,b): the data for f, where a is the angle
    !              and b is the number of elements
    !INOUT  xp: data for x-axis
    !   coeffs: coefficients for multiple splines
    !---------------------------------------------
    subroutine initPresentation(xp,coeffs,P,np,dimy)
        integer, intent(in) :: np,dimy
        real(kind=dp), intent(inout) :: xp(:)            !points in x-axis
        real(kind=dp), intent(inout) :: coeffs(:,:,:)    !coefficients
        real(kind=dp), intent(in) :: P(dimy,np)                           !points in x-axis in
        real(kind=dp) :: val
        real(kind=dp) :: tmpP,tmpVal
        integer :: m,j1,j2,ind
        
        do j1=1,np
            xp(j1) = P(1,j1)*pi/180.0_dp
        enddo
        

        !valgrind requires initialization of these or otherwise
        !it detects false access
        coeffs(:,:,:) = 0.0_dp
        !Store d_k
        do j1=1,np
            do j2=2,dimy
                coeffs(4,j2-1,j1)=P(j2,j1)
            enddo
        enddo
        
        call cubicspline(xp,coeffs,0.0_dp,0.0_dp,np)
    end subroutine

    !---------------------------------------------
    !cubicspline(xp,coeffs,dfa,dfb)
    !Creates multiple clamped cubic spline presentations
    !to coeffs(1:3,:,1:np) from coeffs(4,:,1:np)
    !and xp(1:np).
    !
    !IN:    xp: x-axis data
    !      dfa: 1st derivative at f(x_0)
    !      dfb: 1st derivative at f(x_n)
    !INOUT: coeffs(): spline's coefficients
    !---------------------------------------------
    subroutine cubicspline(xp,coeffs,dfa,dfb,np)
        integer,intent(in) :: np
        real(kind=dp),intent(in) :: xp(np)            !x-axis data
        real(kind=dp),intent(inout) :: coeffs(4,6,np)   !coeffs, coeffs(4,:,:) should contain y-axis data
        real(kind=dp), intent(in) :: dfa,dfb      !derivatives at the end points                                
        real(kind=dp),allocatable :: h(:),a(:),aa(:),bb(:),cc(:),x(:)
        integer :: n,j1,j2
        
        n = np
        allocate(h(n),a(n),aa(n),bb(n),cc(n),x(n))
        do j1=1,n-1
            h(j1)=xp(j1+1)-xp(j1)
        enddo
       
        do j2=1,size(coeffs,dim=2)
            a(1)=3.0_dp*(coeffs(4,j2,2)-coeffs(4,j2,1))/h(1)-3.0_dp*dfa
            a(n)=3.0_dp*dfb-3.0_dp*(coeffs(4,j2,n)-coeffs(4,j2,n-1))/h(n-1)

            do j1=2,n-1       
                a(j1)=3.0_dp*(coeffs(4,j2,j1+1)-coeffs(4,j2,j1))/h(j1)&
                &   -3.0_dp*(coeffs(4,j2,j1)-coeffs(4,j2,j1-1))/h(j1-1)
            enddo

            do j1=1,n-1
                aa(j1)=h(j1)
                cc(j1)=h(j1)
            enddo

            bb(1) = 2.0_dp*h(1)
            do j1=2,n-1
                bb(j1)=2.0_dp*(h(j1)+h(j1-1))
            enddo
            bb(n) = 2.0_dp*h(n-1)
            call thomasAlgorithm(x,aa,bb,cc,a)
            
            do j1=1,n-1   
                coeffs(1,j2,j1)=(x(j1+1)-x(j1))/(3.0_dp*h(j1))
                coeffs(2,j2,j1)=x(j1)
                coeffs(3,j2,j1)=(coeffs(4,j2,j1+1)-coeffs(4,j2,j1))/h(j1)-&
                    & (2.0_dp*x(j1)+x(j1+1))*h(j1)/3.0_dp 
            enddo
        enddo
        deallocate(h,a,aa,bb,cc,x)
    end subroutine

    !---------------------------------------------
    !interpolate(tmpP,xp,coeffs,the)
    !Find tmpP=f(the) where f is expressed as a spline
    !All the splines are evaluated
    !
    !IN:    the:    =x in f(x)
    !    coeffs:    spline coefficients
    !        xp:    x-axis data
    !OUT:  tmpP:    multiple f(the) 
    !
    !---------------------------------------------
    pure subroutine interpolate(tmpP,xp,coeffs,the)
        real(kind=dp), intent(out) :: tmpP(:)         !interpolated values
        real(kind=dp), intent(in) :: xp(:)            !x-axis data
        real(kind=dp), intent(in) :: coeffs(:,:,:)    !spline coefficients
        real(kind=dp), intent(in) :: the              !point where spline is evaluated
        real(kind=dp) :: tmpX
        integer :: ind
        if(the<xp(size(xp))) then          
            !ind = bisectionSearch(xp,the)
            ind=1+int(the*(size(xp)-1)/pi,kind=dp)
            tmpX = the-xp(ind)      
            
        else
            ind=size(coeffs,dim=3)
            tmpP(:)=coeffs(4,:,ind)
            return
        endif
        tmpP(:)=coeffs(1,:,ind)*tmpX**3+coeffs(2,:,ind)*tmpX**2+   &
            & coeffs(3,:,ind)*tmpX+coeffs(4,:,ind)
    end subroutine
    
    !---------------------------------------------
    !interpolate2(tmpP,xp,coeffs,the)
    !Find tmpP=f(the) where f is expressed as a spline
    !Only the first spline is evaluated
    !
    !IN:    the:    =x in f(x)
    !    coeffs:    spline coefficients
    !        xp:    x-axis data
    !OUT:  tmpP:    =f(the)
    !---------------------------------------------
    pure subroutine interpolate2(tmpP,xp,coeffs,the)
        real(kind=dp), intent(out) :: tmpP            !interpolated values
        real(kind=dp), intent(in) :: xp(:)            !x-axis data
        real(kind=dp), intent(in) :: coeffs(:,:,:)    !spline coefficients
        real(kind=dp), intent(in) :: the              !point where spline is evaluated
        real(kind=dp) :: tmpX
        integer :: ind
        
        if(the<xp(size(xp))) then
            ind=1+int(the*(size(xp)-1)/pi,kind=dp)
            !ind = bisectionSearch(xp,the)
            tmpX = the-xp(ind) 
        else
            ind=size(coeffs,dim=3)
            tmpP=coeffs(4,1,ind)
            return
        endif
        tmpP=coeffs(1,1,ind)*tmpX**3+coeffs(2,1,ind)*tmpX**2+   &
            & coeffs(3,1,ind)*tmpX+coeffs(4,1,ind)
    end subroutine

    !---------------------------------------------
    !bisectionSearch(xp,the)
    !Search index for the value(lower) closest to
    !the value "the" using bisection search.
    !
    !IN: xp:    data   
    !   the:    the value to look for
    !---------------------------------------------
    pure function bisectionSearch(xp,the) result(retVal)
        real(kind=dp), intent(in) :: the         !point to be evaluated
        real(kind=dp), intent(in) :: xp(:)       !sorted array
        real(kind=dp) :: errCheck
        integer :: retVal
        integer :: x1,x2,xmid 
        x1=1
        x2=size(xp)
        do while(x2-x1>1)
            xmid=(x1+x2)/2
            if(xp(xmid) > the) then
                x2=xmid
            else
                x1=xmid
            endif
        enddo
        retVal=x1
    end function
        
    
    !---------------------------------------------
    !thomasAlgorithm(x,a0,b0,c0,d0)
    !solves tridiagonal matrix
    !
    !IN: a0,b0,c0:  diagonal vectors
    !          d0:  right side column vector
    !OUT:       x:  result
    !--------------------------------------------
    subroutine thomasAlgorithm(x,a0,b0,c0,d0)
        real(kind=dp), intent(in) :: a0(:),b0(:),c0(:),d0(:)
        real(kind=dp), intent(out) :: x(:)
        real(kind=dp), allocatable :: c(:),d(:)
        integer :: n,j1
        n = size(d0)
        allocate(c(n),d(n))
        c(1) = c0(1)/b0(1)
        do j1=2,n-1
            c(j1)=c0(j1)/(b0(j1)-a0(j1-1)*c(j1-1))
        enddo
        d(1)=d0(1)/b0(1)
        do j1=2,n
            d(j1)=(d0(j1)-a0(j1-1)*d(j1-1))/(b0(j1)-a0(j1-1)*c(j1-1))
        enddo
        x(n)=d(n)
        do j1=n-1,1,-1
            x(j1)=d(j1)-c(j1)*x(j1+1)
        enddo
        deallocate(c,d)
    end subroutine

    !------------------------------------------
    !Modified Bisection method which finds a root 
    !for f(x)-rn=0
    !
    !Input:    ncm: random points total
    !              rn: how much a function has 
    !                   cumulated at the point      
    !           x01,x02: interval
    !               
    !
    !Output:   retVal: solution to f(x)-rn=0
        !------------------------------------------
    function bisectionMod(rn,x1,x2,tol,ncm,xp,coeffs,np) result(retVal)
        integer, intent(in) :: np
        real(kind=dp), intent(in) :: xp(np),coeffs(4,6,np)
        
        real(kind=dp), intent(in) :: rn             !see manual
        real(kind=dp), intent(in) :: tol            !tolerance
        real(kind=dp), intent(in) :: x1,x2          !lower and upper bound
        integer, intent(in) :: ncm                  !number of integration points

        real(kind=dp) :: dx,fmid
        real(kind=dp) :: a,b,c,retVal,f,y1
        integer :: j1,nmax

        nmax = 60
        if(x1>=x2) then
            return
        endif
        a=x1
        b=x2
        y1=PCDFSPLI(a,ncm,xp,coeffs,np)-rn
        do j1=1,nmax
            c=(a+b)*0.5_dp
            f=PCDFSPLI(c,ncm,xp,coeffs,np)-rn
            if(f==0.0_dp .or. (b-a)*0.5_dp<tol) then
                retVal=c
                return            
            endif
            if(sign(1.0_dp,f)==sign(1.0_dp,y1)) then
                a=c
                y1=f        
            else
                b=c
            endif 
        enddo
        
        return
        !iee_denormal
    end function


    !---------------------------------------------
    !pcdfspli(xi,ncm)
    !
    !IN: ncm:   number of random points
    !     xi:   point f(x1) 
    !
    !Computes the cumulative distribution function for the spline
    !scattering phase function of the spline scattering phase matrix
    !---------------------------------------------
    function pcdfspli(xi,ncm,xp,coeffs,np) result(retVal)
        integer, intent(in) :: np
        real(kind=dp), intent(in) :: xp(np),coeffs(4,6,np)
        integer,intent(in) :: ncm       !number of integration points
        real(kind=dp), intent(in) :: xi !x-value
        integer :: j1,nnt
        real(kind=dp) :: retVal
        real(kind=dp) :: P
        real(kind=dp), allocatable :: X(:),WX(:)
        nnt=16+int(0.5_dp*(xi+1.0_dp)*ncm)
        allocate(X(nnt),WX(nnt))
        call gaussLQuads(-1.0_dp,xi,X,WX)
        retVal=0.0_dp
        do j1=1,nnt
            call interpolate2(P,xp,coeffs,acos(X(j1)))
            retVal=retVal+WX(j1)*0.5_dp*P
        enddo
        deallocate(X,WX)
    end function



end module
