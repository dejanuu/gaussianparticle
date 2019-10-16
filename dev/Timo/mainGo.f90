program mainGo
use objsupport
! Traces rays and Mueller matrices for gaussian random spheres.

  implicit none

    integer, parameter ::                              &
      sp = kind(1.0),                                  &
      dp = selected_real_kind(2*precision(1.0_sp)),    &
      qp = selected_real_kind(2*precision(1.0_dp))
  integer  :: nbin, j, jray, jsub, jbox, pout, nray, seed1, &
    ntri,j1,j2,totref,lmax, npar,pin,pinmax,poutmax,ntr,nis, savenum, &
        l1, flgin, np, nrn, pflg,ncm, nnod, nsub
  integer, dimension(50) :: PBOX, NISBOX
  integer :: run_num
  real (kind=dp), dimension(4,4) :: FOUT, FIN  
  real (kind=dp), dimension(2,50) :: MABOX  
  real (kind=dp), dimension(3,50) :: KEBOX, KFBOX, XBOX
    real (kind=dp), dimension(361) :: XP1
  real (kind=dp), dimension(361,4,4) :: S, YP1, YP21
    real (kind=dp), dimension(0:360,4,4) :: P1
  real (kind=dp), dimension(4,4,50) :: FBOX
  real (kind=dp):: wrayh, wraym, norm, ran2, mui, nui, phii, sphii, cphii, &
    rmax, r0, phi0, qsca,qabs,qin,qout,qis,qbox,abscf,len,lenabs,kekf,beta,nie, &
       nk,Fstop,bin,bin0,norml,normr,wray, wavelen, xa, m2real, &
       m2imag, qscaG, qabsG, qextG, renorm, lin, otin, lenrn, ghg1, ghg11, ghg21, &
        whg1, pmx1, dthe, pnorm1, radius, mmed, sig, nuc, rho, cs2d, cs4d, &
        ell,fintmp,volume0
    real(kind=dp), parameter :: pi = 4.0_dp*atan(1.0_dp)  
  real (kind=dp), dimension(2):: MA0, MAOUT, MAIN
    real (kind=dp), dimension(0:1000) :: CSRN1
    real (kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD
    real (kind=dp), dimension(0:256) :: CSCF
    real (kind=dp), dimension(130000) :: MUN, PHIN

    real (kind=dp), dimension(130000,3) :: XN, XN0
    real (kind=dp), dimension(260000,3) :: NT, NT0
    integer , dimension(260000,3) :: IT
    integer :: indA,indB
  complex (kind=dp):: m1,m2
    real (kind=dp), dimension(3) :: KEOUT, KFOUT, X, ELOUT, EROUT, Y, N, KEIN, KFIN, ELIN, ERIN
    complex (kind=dp), dimension(3) :: HLOUT, HROUT, HLIN, HRIN
    complex (kind=dp), dimension(3,50) :: HLBOX, HRBOX    
    character (len=16) :: fname 
    character (len=16) :: fname2     
    character (len=1) :: inhmg
    character (len=128) :: infile
    character(len=64) :: prefix,postfix
    character(64),allocatable :: fnames(:)
    integer :: indx
    real(kind=dp) :: newV

    newV = 4.0_dp/3.0_dp*pi*100.0_dp**3
    !allocate(fnames(100))
    !open(unit=1,file="list")
    !do j1=1,100
    !    read(1,'(A)') fnames(j1)
    !enddo
    !close(1)


    prefix="meshes/sirisg"
    postfix=".off"
! Read inputs from file:
    IF(IARGC() .EQ. 0) THEN
        infile = 'Input.in'
    ELSE
        CALL GETARG(1,infile)
    ENDIF
    open(unit=1, file=infile, status='old')


    read (1, *) nray     ! Number of rays
    read (1, *) npar        ! Number of sample particles.
    read (1, *) Fstop    ! Minimum relative flux
    read (1, *) pinmax    ! Max internal chord 
    read (1, *) poutmax    ! Max external chord
    read (1, *) nbin    ! Number of scattering angle bins
    read (1, *) mmed        ! Refractive index of the medium
    read (1, *) m2real     ! Particle refractive index (real)
    read (1, *) m2imag    ! Particle refractive index (imag)
    read (1, *) seed1     ! Seed for random number generation 
    read (1, *) wavelen   ! Wavelength
    read (1, *) radius      ! Mean radius (microns)
    read (1, *) run_num    ! Identification number for run -> output file is ice_go_'run'.out 
    read (1, *) inhmg    ! use inhomogeneous waves (y/n)  
    read (1, *) flgin       ! Internal medium: 1=yes
    read (1, *) pflg        ! Scattering matrix: input(1), Rayleigh(2)
    read (1, *) otin        ! Internal medium: single-scattering albedo
    read (1, *) lin         ! Internal medium: mean free path (in micrometers)
    read (1, *) ghg1        ! Diffuse-medium asymmetry parameter (internal)
    read (1, *) ghg11       ! Diffuse-medium forward asymmetry (internal)
    read (1, *) ghg21       ! Diffuse-medium backward parameter (internal)
    read (1, *) pmx1        ! Diffuse-medium single-scattering max. polarization
    read (1, *) np          ! Number of phase matrix angular points
    read (1, *) nrn         ! Number of points in random number array
    read (1, *) ncm         ! Number of G-L integration points for cosine map
    read (1, *) sig         ! Relative standard deviation of radius.
    read (1, *) nuc         ! Power law index for C_3 correlation.
    read (1, *) lmax        ! Maximum degree in C_1, C_2, C_3.
    read (1, *) ntr         ! Number of triangle rows in an octant.

    close (1)
    write(6,*) flgin
  write (*,*) 'Geometric optics approximation for a Gaussian random sphere...'
    write (*,'(A,F7.1,A)') 'Particle mean radius',radius,'micrometers'
  write (*,'(A,I12)') 'Number of rays: ', nray

  
    call init_random(seed1)


    nsub=nray/npar
    if (nsub*npar.ne.nray) stop 'Trouble in GEO: ray/particle number mismatch.'

!   Initialization of the Gaussian random sphere:

    beta=sqrt(log(sig**2+1.0_dp))
    call cs1cf(CSCF,nuc,2,lmax)  !power law
    call csini(CSCF,ell,cs2d,cs4d,lmax)

    call sgscfstd(SCFSTD,CSCF,beta,lmax)

    rho=beta/ell


!   Refractive indices:

  m2 = cmplx(m2real, m2imag,KIND=dp)
  m1 = cmplx(mmed, 0.0_dp,KIND=dp)   ! ref index of the medium

    lin=lin/radius


    whg1 = (ghg1-ghg21)/(ghg11-ghg21)
  
! Size parameter, absorption parameter, refractive indices:


    xa = 2.0_dp*pi*radius/wavelen
    abscf = -2.0_dp*xa*m2imag

! Discretization:

    call TRIDS(MUN,PHIN,IT,nnod,ntri,ntr)

! Initialize:

  qsca = 0.0_dp
  qabs = 0.0_dp
  qin = 0.0_dp
  qout = 0.0_dp
  qis = 0.0_dp
  qbox = 0.0_dp
  S = 0.0

! Initialize the scattering phase matrices:

    dthe=pi/np
    if (pflg .eq. 1) then
        call pmatrix1(P1,np)
        call psplivi(P1,CSRN1,XP1,YP1,YP21,pnorm1,dthe,nrn,np,ncm)
    end if

! Storage grid parameters:

    bin = pi/nbin
    bin0 = cos(0.5_dp*bin)

! Refractive indices:

    call refindapp(MA0,m1,1.0d0)





! Externally propagating rays. Initialization:

    jray = 0
    jsub = 0
    jbox = 0
    wrayh = 0.0_dp ! rays that hit
    wraym = 0.0_dp ! rays that miss

    !call sgscf(ACF,BCF,SCFSTD,lmax)
    call sgscf(ACF,BCF,SCFSTD,lmax)
    call rgstd(XN0,NT0,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)



    indx=1
    !call readMesh2(XN0,NT0,IT,ntri,nnod,rmax,fnames(indx),newV)
    call WriteMeshToFileOff(prefix,indx,postfix,XN0,IT)
    indx=indx+1



100  if (jbox .eq. 0) then
    pout = 0
    jray = jray + 1
   
    ! Normalizations:
    if (jray .gt. nray) then ! number of rays done --> quit program
           norm = nray/wrayh
           qsca = qsca*norm
           qabs = qabs*norm
           qin = qin*norm
           qout = qout*norm
           qis = qis*norm
           qbox = qbox*norm
           S = S*norm
           
           qscaG = qsca/nray           
         qabsG = qabs/nray
         qextG = qscaG + qabsG 
           
           renorm = 1.0_dp/qscaG
           norm = renorm*2.0_dp/(nray*(1.0_dp-bin0))
           
           do j1 = 1, 4
            do j2 = 1, 4
               S(1,j1,j2) = norm*S(1,j1,j2) ! forward scattering
               S(nbin+1,j1,j2) = norm*S(nbin+1,j1,j2) ! backscattering
        end do
      end do
           
           !do j = 1, nbin-1
           do j = 1, nbin-1           
            norm = renorm*2.0_dp/(nray*(cos((j-0.5_dp)*bin)-cos((j+0.5_dp)*bin)))      
            do j1 = 1, 4
               do j2 = 1, 4
                  S(j+1,j1,j2) = norm*S(j+1,j1,j2)
          end do
        end do
      end do    
           
           write (*,*) '----------------'
           write (*,*) 'Results:'
           write (*,'(A,E12.4)') 'Qsca = ', qscaG
           write (*,'(A,E12.4)') 'Qabs = ', qabsG
           write (*,'(A,F8.5)') 'SSA = ', qscaG/qextG
           write (*,'(A,F20.1)') 'Rays hit = ', wrayh
           write (*,'(A,F20.1)') 'Rays miss = ', wraym
  
! Hovenier's tests for scattering matrix: 

!      write (*,*) 'Hovenier tests for the scattering matrix:'
!      write (*,'(A,F15.7)') '(S22(0)-S33(0))/S11(0) = ', (S(1,2,2)-S(1,3,3))/S(1,1,1)
!      write (*,'(A,F15.7)') '(S22(pi)+S33(pi))/S11(pi) = ', (S(nbin+1,2,2)+S(nbin+1,3,3))/S(nbin+1,1,1)
!      write (*,'(A,F15.7)') 'S12(0)/S11(0) = ', S(1,1,2)/S(1,1,1)
!      write (*,'(A,F15.7)') 'S12(pi)/S11(pi) = ', S(nbin+1,1,2)/S(nbin+1,1,1)
!      write (*,'(A,F15.7)') 'S34(0)/S11(0) = ', S(1,3,4)/S(1,1,1)
!      write (*,'(A,F15.7)') 'S34(pi)/S11(pi) = ', S(nbin+1,3,4)/S(nbin+1,1,1)  
!      write (*,'(A,F15.7)') '(S11(pi)-2S22(pi)-S44(pi))/S11(pi) = ', &
!        (S(nbin+1,1,1)-2.0*S(nbin+1,2,2)-S(nbin+1,4,4))/S(nbin+1,1,1)                  

       write (fname,'(A,I3.3,A)') 'outputS_',run_num,'.out' ! output
       write (*,*) 'Scattering matrix elements written in file ', fname
       
        open (unit = 2, file = fname)
         savenum = 0
       do l1 = 1, nbin+1 ! vanha oli nolla...nbin  
              write (2,'((F5.1),1X,16(E16.7, 1X))') (l1-1)*180.0/nbin, S(l1,1,1), S(l1,1,2), S(l1,1,3), S(l1,1,4), &
                            S(l1,2,1), S(l1,2,2), S(l1,2,3), S(l1,2,4), &
                            S(l1,3,1), S(l1,3,2), S(l1,3,3), S(l1,3,4), &
                            S(l1,4,1), S(l1,4,2), S(l1,4,3), S(l1,4,4)     
        savenum = savenum + 1
        end do
        
        close (2)
        
        write (fname2,'(A,I3.3,A)') 'outputQ_',run_num,'.out' ! output
       write (*,*) 'Scattering efficiencies written in file ', fname2
       
        open (unit = 3, file = fname2)
        write (3,'(A,F15.8)') 'qsca = ',qscaG
        write (3,'(A,F15.8)') 'qabs = ',qabsG
        write (3,'(A,F15.8)') 'qext = ',qextG
        write (3,'(A,F15.8)') 'rhit = ',sqrt(wrayh/nray)
        write (3,'(A,F15.8)') 'qsca/qext = ',qscaG/qextG
            write (3,'(A,F15.8)') 'wavelen = ',wavelen
            write (3,'(A,F15.8)') 'meanRadius = ',radius
            write (3,'(A,F15.8)') 'mmed = ',mmed
            write (3,'(A,F15.8)') 'mreal = ',m2real
            write (3,'(A,F15.8)') 'mimg = ',m2imag
            write (3,'(A,F15.8)') 'meanFreePath = ',lin*radius
            write (3,'(A,F15.8)') 'qin = ',qin/nray
            write (3,'(A,F15.8)') 'qout = ',qout/nray
            write (3,'(A,F15.8)') 'qis = ',qis/nray
            write (3,'(A,F15.8)') 'qbox = ',qbox/nray

            close (3)
           
           stop
        end if

! rayinic initializes the Mueller matrix, wave vector, and 
! parallel/perpendicular coordinate axes.

        MAOUT = MA0 
    
        call rayinic(FOUT,KEOUT,KFOUT,HLOUT,HROUT)


    call random_number(ran2)
        mui = 1.0-2.0*ran2 ! -1...1 ! tama on cos(theta), tasaisesti jakautunut -1...1
    !mui = 0.0 
        nui = sqrt(1.0 - mui**2) ! taman on sin(theta)
    call random_number(ran2)          
        phii = 2.0*pi*ran2
    !phii = 0.0
        cphii = cos(phii)
        sphii = sin(phii)

        
! Rotation of KE, KF, HL, and HR to the particle coordinate system:
       
        call raypart0(KEOUT,mui,nui,cphii,sphii)
        call raypart0(KFOUT,mui,nui,cphii,sphii)
        call raypart0c(HLOUT,mui,nui,cphii,sphii)
        call raypart0c(HROUT,mui,nui,cphii,sphii)
        
! Generation of the spherical harmonics coefficients for the
! logradius. Discretization and random orientation of the Gaussian
! sample sphere:

        jsub=jsub+1
        if (jsub.gt.nsub) then

            call sgscf(ACF,BCF,SCFSTD,lmax)
            call rgstd(XN0,NT0,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)
            !call readMesh2(XN0,NT0,IT,ntri,nnod,rmax,fnames(indx),newV)
            call WriteMeshToFileOff(prefix,indx,postfix,XN0,IT)
            indx=indx+1
            jsub=1
        endif

        call pranor(XN,NT,XN0,NT0,nnod,ntri)





! Weight of incident ray.

        wray = rmax**2
        
! Offset of incident ray:

    call random_number(ran2)
        r0 = rmax*sqrt(ran2)
!    r0 = 0.01
    call random_number(ran2)
        phi0 = 2.0*pi*ran2
!    phi0 = 5.0
        X = (/r0*cos(phi0), r0*sin(phi0), -rmax/) ! - vai ei z-komponenttiin...?
        
        call raypart0(X,mui,nui,cphii,sphii) ! X is the location where the incident ray starts (particle coordinates)

! istri determines the possible interaction point on the particle surface.

    nie = 1.0
  
        call istri(X,KEOUT,N,XN,NT,IT,nk,len,ntri,nis,nie)
        
        if (nis .eq. 0) then
          wraym = wraym+wray
           goto 100
        endif
        
        wrayh = wrayh+wray
        do j1 = 1, 4
          do j2 = 1, 4
              FOUT(j1,j2) = wray*FOUT(j1,j2)
      end do
    end do
      
  else

! RAYGETD extracts one ray from the pile of external rays. ISTRI 
! determines whether that ray further interacts with the particle 
! surface. If it doesnt, SCATTER stores the ray among other scattered 
! rays.

      call raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT,HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    
      jbox = jbox - 1
      
      if (FOUT(1,1) .le. wray*Fstop .or. pout .gt. poutmax) then
          qout = qout+FOUT(1,1)
           goto 100
        end if
        N = NT(nis, 1:3)
        call prosca(nk,N,KEOUT)

  end if
  
! Scattering or external incidence:
 
    if (nk .ge. 0.0) then ! jos kyseessa sisaltapain ulostulo
        norml = 0.0
        do j1 = 1, 3
            KFOUT(j1)=KEOUT(j1)
            norml=norml+dreal(HLOUT(j1))**2
        end do

        norml = sqrt(norml)


        do j1 = 1, 3
            ELOUT = real(HLOUT(1:3))/norml
        end do
        call provec(EROUT,KEOUT,ELOUT)

        do j1=1,3
            HLOUT(j1)=dcmplx(ELOUT(j1),0.0_dp)
            HROUT(j1)=dcmplx(EROUT(j1),0.0_dp)
        end do
        MAOUT(1)=MA0(1)
        MAOUT(2)=MA0(2)

        nie = 1.0
        !write(6,*) X
        call istri(X,KEOUT,N,XN,NT,IT,nk,len,ntri,nis,nie) ! osuuko enaa uuteen kolmioon?

        !!!!!if (nis .ne. 0) write (*,*) 'Triangle ID at external->internal (should be 0!!!): ', nis
        if (nis .eq. 0) then ! ei osu. tallennetaan sade sirontamatriisiin.

           call partray0(KEOUT,mui,nui,cphii,sphii)
           call partray0(ELOUT,mui,nui,cphii,sphii)
           call partray0(EROUT,mui,nui,cphii,sphii)
           call scatter(S,FOUT,KEOUT,ELOUT,EROUT,qsca,bin,bin0,nbin)
           goto 100
        end if
    end if

! INCIDE determines the reflected and refracted Mueller matrices, 
! wave vectors, and parallel/perpendicular coordinate axes.

    pout = pout+1
    if (inhmg .eq. 'y') then

       
      call incide(FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
                  FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                  N,m1,m2,totref)

!    else if (inhmg .eq. 'n') then    
!      call incide_trad(FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
!                  FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
!                  N,m1,m2,totref)                
  end if
! Internally propagating rays. RAYPUTD stores the reflected ray into 
! the pile of external rays.

    pin = pout  
300 jbox = jbox+1
    if (jbox .gt. 49) then
      qbox = qbox+FIN(1,1)
        goto 100
    end if

    call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
           MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
           HLOUT,HROUT,X,MAOUT,pout,nis,jbox)

    if (pinmax .eq. 0) then
      qabs = qabs+FIN(1,1)
        goto 100
    end if


310 if (FIN(1,1) .le. wray*Fstop .or. pin .ge. pinmax) then
      qin = qin+FIN(1,1)
        goto 100
    end if

! ISTRI determines the next interaction point for the internal ray, 
! computes the inner normal at the interaction point, ABSORB determines the 
! amount of absorbed flux, and INCIDE determines the reflected and refracted 
! Mueller matrices, wave vectors, and parallel/perpendicular coordinate axes.

    Y = X
  
  nie = -1.0
    call istri(Y,KEIN,N,XN,NT,IT,nk,len,ntri,nis,nie)
    !!!!!write (*,*) 'Interaction number (pin): ', pin      
    !!!!write (*,*) 'Triangle ID at internal interaction (in->out): ', nis    
    !!!write (*,*) 'nis in->out = ', nis
    if (nis .eq. 0) then ! tahan ei kai pitais paatya mitaan? jos paatyy, niin kolmioinnissa vika. (tjs.)
      qis = qis+FIN(1,1)
        goto 100
    end if

!-------------------
!Internal diffuse medium
    if (flgin.eq.1) then
        call random_number(ran2)
        lenrn=-lin*log(ran2)
        if (lenrn .lt. len) then
        
            pin=pin+1
            pout=pout+1

            do j1 = 1, 3
                X(j1)=X(j1)+lenrn*KEIN(j1)
            end do


            call prosca(kekf,KEIN,KFIN)
            lenabs = lenrn*abs(kekf)

            call absorb(FIN,qabs,lenabs,abscf)

            fintmp=FIN(1,1)
            if (FIN(1,1) .le. wray*Fstop) then
                qin=qin+fintmp
                goto 100
            endif



            call absorbrt(FIN,qabs,otin)



            norml = 0.0
            do j1 = 1, 3
                norml=norml+dreal(HLIN(j1))**2
            end do

            norml = sqrt(norml)

            ELIN = real(HLIN(1:3))/norml

            call provec(ERIN,KEIN,ELIN)






!            norml = 0.0
!            normr = 0.0
!            do j1 = 1, 3
!                norml=norml+dreal(HLIN(j1))**2
!                normr=normr+dreal(HRIN(j1))**2
!            end do
!
!            norml = sqrt(norml)
!            normr = sqrt(normr)
!
!
!            ELIN = real(HLIN(1:3))/norml
!            ERIN = real(HRIN(1:3))/normr


            call scart(FIN,KEIN,ERIN,ELIN,CSRN1,XP1,YP1,YP21,pmx1,np,nrn,pflg)

            do j1 = 1, 3
                KFIN(j1)=KEIN(j1)
                HLIN(j1)=cmplx(ELIN(j1),0.0_dp,KIND=dp)
                HRIN(j1)=cmplx(ERIN(j1),0.0_dp,KIND=dp)
            end do



            goto 310
        endif


    endif

!------------------

    X = Y
    call prosca(kekf,KEIN,KFIN)

    lenabs = len*abs(kekf)

    call absorb(FIN,qabs,lenabs,abscf)



    call incide(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                  FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
                  N,m2,m1,totref)
    pin = pin+1
    pout = pout+1
  if (totref .eq. 1) goto 310 ! total internal reflection
  if (totref .eq. -1) goto 300 

contains

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
    integer :: j1,j2
    O=(/0.0_dp,0.0_dp,0.0_dp/)
    volume=0.0_dp
    do j1=1,ntri
        B=XN0(IT(j1,1),:)
        C=XN0(IT(j1,2),:)
        D=XN0(IT(j1,3),:)
        volume=volume+simplex3volume(B,D,C)
    end do
end subroutine

subroutine WriteMeshToFile(prefix,indx,postfix,XN0,IT)
    character(len=64), intent(in) :: prefix,postfix
    integer, intent(in) :: IT(:,:),indx
    real(kind=dp), intent(in) :: XN0(:,:)
    character(len=74) :: fname
    integer :: indA,indB,j1,j2
    !return
    indA=0
    indB=0
    do j1=1,30000
        if(dot_product(XN0(j1,:),XN0(j1,:))==0.0) cycle
        indA=indA+1
    enddo
    do j1=1,30000
        if(sum(IT(j1,:))==0) cycle
        indB=indB+1
    enddo

    write(fname,'(A,I4.4,A)') adjustl(trim(prefix)),indx,adjustl(trim(postfix))
     
    open(unit=123,file=fname)
    write(123,*) "# vtk DataFile Version 2.0"
    write(123,*) "gsphere output"   
    write(123,*) "ASCII"         
    write(123,*) "DATASET POLYDATA"   
    write(123,*) "POINTS    ",indA,"  float"
    do j1=1,indA
        if(dot_product(XN0(j1,:),XN0(j1,:))==0.0) cycle
        write(123,*) XN0(j1,:)
    enddo
    write(123,*) "POLYGONS ",indB,indB*4
    do j1=1,indB
        if(sum(IT(j1,:))==0) cycle
        write(123,*) 3,IT(j1,:)-1
    enddo

    close(123)
    write(6,*) "mesh written to:", fname
end subroutine


subroutine WriteMeshToFileOff(prefix,indx,postfix,XN0,IT)
    character(len=64), intent(in) :: prefix,postfix
    integer, intent(in) :: IT(:,:),indx
    real(kind=dp), intent(in) :: XN0(:,:)
    character(len=74) :: fname
    integer :: indA,indB,j1,j2
    !return
    indA=0
    indB=0
    do j1=1,30000
        if(dot_product(XN0(j1,:),XN0(j1,:))==0.0) cycle
        indA=indA+1
    enddo
    do j1=1,30000
        if(sum(IT(j1,:))==0) cycle
        indB=indB+1
    enddo

    write(fname,'(A,I4.4,A)') adjustl(trim(prefix)),indx,adjustl(trim(postfix))
    
    open(unit=123,file=fname)
    
    write(123,'(A)') "OFF"
    write(123,'(A)') "# made by siris"  
    write(123,'(A)') "# A cube"         
    write(123,*)
    write(123,*) indA,indB,1
    do j1=1,indA
        if(dot_product(XN0(j1,:),XN0(j1,:))==0.0) cycle
        write(123,*) XN0(j1,:)
    enddo
    do j1=1,indB
        if(sum(IT(j1,:))==0) cycle
        write(123,*) 3,IT(j1,:)-1
    enddo

    close(123)
    write(6,*) "mesh written to:", fname
end subroutine


end program mainGo
