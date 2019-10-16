program mainGo

! Traces rays and Mueller matrices for gaussian random spheres.

  use common
  implicit none

  integer  :: nbin, j, jray, jsub, jbox, jboxin, pout, nray, seed1,ntri,&
  j0,j1,j2,j3,totref,lmax, npar,pin,pinmax,poutmax,ntr,&
  nis, nis2, savenum, l1, rtflg1, rtflg2, np, nrn,ncm, nnod, nsub, nie, nie2
  integer, dimension(500) :: PBOX, NISBOX, PINBOX, NISINBOX
  integer :: run_num
  real (kind=dp), dimension(4,4) :: FOUT, FIN, FCORE  
  real (kind=dp), dimension(2,500) :: MABOX, MAINBOX
  real (kind=dp), dimension(3,500) :: KEBOX, KFBOX, XBOX, KEINBOX, KFINBOX, XINBOX
  real (kind=dp), dimension(361) :: XP1, XP2
  real (kind=dp), dimension(361,4,4) :: S, YP1,YP2, YP21, YP22
  real (kind=dp), dimension(0:360,4,4) :: P,P1,P2
  real (kind=dp), dimension(4,4,500) :: FBOX, FINBOX

  real (kind=dp):: wrayh, wraym, norm, ran2, mui, nui, phii, sphii,&
  cphii,rmax, r0,phi0,qsca,qabs,qin,qout,qis,qbox,abscf1,abscf2,len,len1,len2,len3,lenabs,&
  kekf,beta,nk1,nk2,Fstop,bin,bin0,norml,wray, wavelen, xa, m1real,m1imag,&
  m2real, m2imag, qscaG, qabsG,qextG, renorm, lenrt1, lenrt2, omg1,omg2, lenrn, dthe,&
  pnorm1, pnorm2, radius1, radius2, mmed, sig, nuc, rho, cs2d, cs4d, ell, rsc

  real (kind=dp), dimension(2):: MA0, MAOUT, MAIN, MACORE
  real (kind=dp), dimension(0:1000) :: CSRN1, CSRN2
  real (kind=dp), dimension(0:256,0:256) :: ACF, BCF, SCFSTD
  real (kind=dp), dimension(0:256) :: CSCF
  real (kind=dp), dimension(130000) :: MUN, PHIN

  real (kind=dp), dimension(130000,3) :: XN1, XN2
  real (kind=dp), dimension(260000,3) :: NT
  integer , dimension(260000,3) :: IT

  complex (kind=dp):: m0,m1,m2
  real (kind=dp), dimension(3) :: KEOUT,KFOUT,X,Y1,Y2, ELOUT,EROUT, N1, N2, KEIN,&
  KFIN, KECORE, KFCORE
  complex (kind=dp), dimension(3) :: HLOUT, HROUT, HLIN, HRIN, HLCORE, HRCORE
  complex (kind=dp), dimension(3,500) :: HLBOX, HRBOX, HLINBOX, HRINBOX
  character (len=32) :: fname 
  character (len=32) :: fname2     
  character (len=128) :: infile

  ! Read inputs from file:
  IF(IARGC() .EQ. 0) THEN
    infile = 'Input_1.in'
  ELSE
    CALL GETARG(1,infile)
  ENDIF
  open(unit=1, file=infile, status='old')


  read (1, *) nray        ! Number of rays
  read (1, *) npar        ! Number of sample particles.
  read (1, *) Fstop       ! Minimum relative flux
  read (1, *) pinmax      ! Max internal chord 
  read (1, *) poutmax     ! Max external chord
  read (1, *) nbin        ! Number of scattering angle bins
  read (1, *) mmed        ! Refractive index of the medium
  read (1, *) m1real      ! Mantle: Particle refractive index (real)
  read (1, *) m1imag      ! Mantle: Particle refractive index (imag)
  read (1, *) m2real      ! Core: Particle refractive index (real)
  read (1, *) m2imag      ! Core: Particle refractive index (imag)
  read (1, *) seed1       ! Seed for random number generation 
  read (1, *) wavelen     ! Wavelength
  read (1, *) radius1     ! Outer radius (microns)
  Read (1, *) radius2     ! Inner radius (microns)
  read (1, *) run_num     ! Identification number for run -> output file is    outputS_'run'.out 
  read (1, *) rtflg1      ! External medium: 1=yes
  read (1, *) rtflg2      ! Internal medium: 1=yes
  read (1, *) omg1        ! External medium: single-scattering albedo
  read (1, *) omg2        ! Internal medium: single-scattering albedo
  read (1, *) lenrt1      ! External medium: mean free path (in micrometers)
  read (1, *) lenrt2      ! Internal medium: mean free path (in micrometers)
  read (1, *) np          ! Number of phase matrix angular points
  read (1, *) nrn         ! Number of points in random number array
  read (1, *) ncm         ! Number of G-L integration points for cosine map
  read (1, *) sig         ! Relative standard deviation of radius.
  read (1, *) nuc         ! Power law index for C_3 correlation.
  read (1, *) lmax        ! Maximum degree in C_1, C_2, C_3.
  read (1, *) ntr         ! Number of triangle rows in an octant.

  close (1)


  write (*,*) 'Geometric optics approximation for a Gaussian random sphere...'
  write (*,'(A,F7.1,A)') 'Particle outer radius',radius1,' micrometers'
  write (*,'(A,F7.1,A)') 'Particle inner radius',radius2,' micrometers'
  write (*,'(A,I0)') 'Number of rays: ', nray

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

  m0 = cmplx(mmed, 0.0_dp,KIND=dp)    ! ref index of the medium
  m1 = cmplx(m1real, m1imag,KIND=dp) ! Mantle
  m2 = cmplx(m2real, m2imag,KIND=dp)  ! Core


  ! Reality checks:

  if (aimag(m0) .ne. 0.0d0) stop 'Trouble in GRO: Im(m0).ne.0.'


  
  ! Size parameter, absorption parameter, refractive indices:

  rsc = radius2/radius1
  xa = 2.0_dp*pi*radius1/(wavelen/real(m0,dp))
  abscf1 = -2.0_dp*xa*aimag(m1/m0)
  abscf2 = -2.0_dp*xa*aimag(m2/m0)

  lenrt1 = lenrt1/radius1
  lenrt2 = lenrt2/radius1

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

  ! Mantle:

  dthe = pi/np    
  if (rtflg1 .eq. 1) then
      call pmatrix1(P1,np)
      call PSPLIVI(P1,CSRN1,XP1,YP1,YP21,pnorm1,dthe,nrn,np,ncm)

!        open (unit=1,file='p1.dat')
!        dthe=pi/np
!        do  j2=0,np
!         read(1,*) XP1(j2+1),(P(j2,j3,1),j3=1,4),&
!                    (P(j2,j3,2),j3=1,4),&
!                    (P(j2,j3,3),j3=1,4),&
!                    (P(j2,j3,4),j3=1,4)
!        end do
!        call PSPLIVI(P,CSRN1,XP1,YP1,YP21,pnorm1,dthe,nrn,np,ncm)
!        close(unit=1)
  end if


  ! Core:

  dthe = pi/np

  if (rtflg2 .eq. 1) then
    call pmatrix2(P2,np)
    call PSPLIVI(P2,CSRN2,XP2,YP2,YP22,pnorm2,dthe,nrn,np,ncm)


!        open (unit=1,file='p2.dat')
!        dthe=pi/np
!        do  j2=0,np
!        read(1,*) XP2(j2+1),(P(j2,j3,1),j3=1,4),&
!                    (P(j2,j3,2),j3=1,4),&
!                    (P(j2,j3,3),j3=1,4),&
!                    (P(j2,j3,4),j3=1,4)
!        end do
!        call PSPLIVI(P,CSRN2,XP2,YP2,YP22,pnorm2,dthe,nrn,np,ncm)
!        close(unit=1)
   end if




  ! Storage grid parameters:

  bin = pi/nbin
  bin0 = cos(0.5_dp*bin)

  ! Refractive indices:

  call refindapp(MA0,m0,1.0_dp)





  ! Externally propagating rays. Initialization:

  jray = 0
  jsub = 0
  jbox = 0
  jboxin = 0
  wrayh = 0.0_dp ! rays that hit
  wraym = 0.0_dp ! rays that miss

  call sgscf(ACF,BCF,SCFSTD,lmax)
  call rgstd(XN1,NT,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)

  do j0 = 1, nnod
    do j1 = 1, 3
      XN2(j0,j1) = rsc*XN1(j0,j1)
    end do
  end do


100 if (jbox .eq. 0) then
    pout = 0
    jray = jray + 1
    
    ! Print status
    IF(MOD(jray,5000) == 0) THEN
      WRITE(*,'(A19,I0,A1,I0)') "  computing at ray ", jray, "/", nray
    END IF


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
      write (*,*) 'Scattering matrix elements written in file ', TRIM(fname)

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



      write (fname,'(A,I3.3,A)') 'pmatrix_new',run_num,'.out' ! output
        write(6,*) fname
      open (unit = 2, file = fname)
      savenum = 0
      do l1 = 1, nbin+1 ! vanha oli nolla...nbin  
        write (2,'((F5.1),1X,16(E16.7, 1X))') (l1-1)*180.0/nbin, S(l1,1,1), S(l1,1,2), S(l1,2,2), S(l1,3,3), S(l1,3,4), S(l1,4,4)     
        savenum = savenum + 1
      end do
      close (2)




      write (fname2,'(A,I3.3,A)') 'outputQ_',run_num,'.out' ! output
      write (*,*) 'Scattering efficiencies written in file ', TRIM(fname2)

      open (unit = 3, file = fname2)
      write (3,'(A,F15.8)') 'qsca = ',qscaG
      write (3,'(A,F15.8)') 'qabs = ',qabsG
      write (3,'(A,F15.8)') 'qext = ',qextG
      write (3,'(A,F15.8)') 'rhit = ',sqrt(wrayh/nray)
      write (3,'(A,F15.8)') 'qsca/qext = ',qscaG/qextG
      write (3,'(A,F15.8)') 'wavelen = ',wavelen
      write (3,'(A,F15.8)') 'OuterRadius = ',radius1
      write (3,'(A,F15.8)') 'InnerRadius = ',radius2
      write (3,'(A,F15.8)') 'mmed = ',mmed
      write (3,'(A,F15.8)') 'm1real = ',m1real
      write (3,'(A,F15.8)') 'm1imag = ',m1imag
      write (3,'(A,F15.8)') 'm2real = ',m2real
      write (3,'(A,F15.8)') 'm2img = ',m2imag
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
    !    mui = 0.0 
    nui = sqrt(1.0 - mui**2) ! taman on sin(theta)
    call random_number(ran2)          
    phii = 2.0*pi*ran2
    !    phii = 0.0
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
      call rgstd(XN1,NT,MUN,PHIN,ACF,BCF,rmax,beta,IT,nnod,ntri,lmax)
      do j0 = 1, nnod
        do j1 = 1, 3
          XN2(j0,j1) = rsc*XN1(j0,j1)
        end do
      end do
      jsub=1
    end if




    ! Weight of incident ray.

    wray = rmax**2

    ! Offset of incident ray:

    call random_number(ran2)
    r0 = rmax*sqrt(ran2)
    !    r0 = 0.01
    call random_number(ran2)
    phi0 = 2.0*pi*ran2
    !    phi0 = 5.0
    X = (/r0*cos(phi0), r0*sin(phi0), -sqrt(rmax**2-r0**2)/) 
    call raypart0(X,mui,nui,cphii,sphii) ! X is the location where the incident ray starts (particle coordinates)

    ! istri determines the possible interaction point on the particle surface.


    nie = 1
    call istri(X,KEOUT,N1,XN1,NT,IT,nk1,len,ntri,nis,nie)

    if (nis .eq. 0) then
      wraym = wraym+wray
      goto 100
    end if
    
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
    
  
    call raygetd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX,MABOX,PBOX,&
                   NISBOX,FOUT,KEOUT,KFOUT,HLOUT,HROUT,X,MAOUT,pout,nis,jbox)
    
    jbox = jbox - 1
      
    if (FOUT(1,1) .le. wray*Fstop .or. pout .gt. poutmax) then
      qout = qout+FOUT(1,1)
       goto 100
    end if
    N1 = NT(nis, 1:3)
    call prosca(nk1,N1,KEOUT)

  end if
  
  ! Scattering or external incidence:
 
  if (nk1 .ge. 0.0_dp) then ! jos kyseessa sisaltapain ulostulo
    call istri(X,KEOUT,N1,XN1,NT,IT,nk1,len,ntri,nis,1)

    if (nis .eq. 0) then
      call ehk(ELOUT,EROUT,HLOUT,HROUT,KEOUT)
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
  pin = pout
  nis2 = nis

  call incide(FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
              FIN,KEIN,KFIN,HLIN,HRIN,MAIN,&
              N1,m0,m1,totref)

  if (totref .eq. 0) then

    jbox = jbox+1
    if (jbox .gt. 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    end if

    call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                   MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jbox)

    jboxin = jboxin+1
    if (jboxin .gt. 499) then
      qbox=qbox+FIN(1,1)
      goto 320
    end if

    call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                   MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                   HLIN,HRIN,X,MAIN,pin,nis,jboxin)

    goto 320

  else if (totref .eq. 1) then

    jbox = jbox+1
    if (jbox .gt. 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    endif

    call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                   MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jbox)

    goto 100

  else

    jboxin = jboxin+1
    if (jboxin .gt. 499) then
      qbox=qbox+FOUT(1,1)
      goto 100
    end if

    call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                   MAINBOX,PINBOX,NISINBOX,FOUT,KEOUT,KFOUT, &
                   HLOUT,HROUT,X,MAOUT,pout,nis,jboxin)

    goto 320

  end if




  ! Internally propagating rays. Ray tracing in the internal medium:


       
320 if (jboxin .eq. 0) then
    goto 100
  end if

  call raygetd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX,&
                 MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN,&
                 HLIN,HRIN,X,MAIN,pin,nis2,jboxin)
       
  jboxin=jboxin-1


  if (FIN(1,1) .le. wray*Fstop .or. pin .gt. pinmax) then
    qin=qin+FIN(1,1)
    goto 320
  end if
     

  ! ISTRI determines the next interaction point for the internal ray, 
  ! computes the inner normal at the interaction point, ABSORB determines the 
  ! amount of absorbed flux, and INCIDE determines the reflected and refracted 
  ! Mueller matrices, wave vectors, and parallel/perpendicular coordinate axes.

  Y1 = X
  Y2 = X

  nie = -1
  nie2 = 1
  call istri(Y1,KEIN,N1,XN1,NT,IT,nk1,len1,ntri,nis,nie)
  call istri(Y2,KEIN,N2,XN2,NT,IT,nk2,len2,ntri,nis2,nie2)

  !!!!!write (*,*) 'Interaction number (pin): ', pin      
  !!!!write (*,*) 'Triangle ID at internal interaction (in->out): ', nis    
  !!!write (*,*) 'nis in->out = ', nis
  if (nis .eq. 0) then ! tahan ei kai pitais paatya mitaan? jos paatyy, niin kolmioinnissa vika. (tjs.)
    qis = qis+FIN(1,1)
    goto 320
  end if

  ! MANTLE INTERSECTION:

  ! Radiative transfer in the mantle:

  if (rtflg1 .eq. 1) then
    call random_number(ran2)
    lenrn = -lenrt1*log(ran2)
    if (lenrn .lt. min(len1,len2)) then
      pin = pin + 1
      pout = pout + 1

      call rtmc(FIN,KEIN,KFIN,HLIN,HRIN,X,XP1,YP1,YP21,CSRN1,&
      qabs,lenrn,abscf1,omg1,np,nrn)


      jboxin=jboxin+1
      if (jboxin .gt. 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

      goto 320
    end if
  end if

  ! Ray tracing in the mantle:

  if (len1 .lt. len2) then

    X = Y1

    call prosca(kekf,KEIN,KFIN)
    lenabs=len1*abs(kekf)
    call absorb(FIN,qabs,lenabs,abscf1)

    if (FIN(1,1).le.wray*Fstop) then
      qin=qin+FIN(1,1)
      goto 320
    endif

    call incide(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                FOUT,KEOUT,KFOUT,HLOUT,HROUT,MAOUT, &
                N1,m1,m0,totref)

    pin=pin+1
    pout=pout+1

    if (totref .eq. 0) then

      jbox=jbox+1
      if (jbox.gt.499) then
        qbox=qbox+FOUT(1,1)
        goto 320
      end if

      call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                     MABOX,PBOX,NISBOX,FOUT,KEOUT,KFOUT, &
                     HLOUT,HROUT,X,MAOUT,pout,nis,jbox)

      jboxin=jboxin+1
      if (jboxin.gt.499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jboxin)

    else if (totref.eq.1) then

      jboxin=jboxin+1
      if (jboxin.gt.499) then
        qbox=qbox+FIN(1,1)
        goto 320
      endif

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jboxin)


    else

      jbox=jbox+1
      if (jbox .gt. 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      endif

      call rayputd1c(FBOX,KEBOX,KFBOX,HLBOX,HRBOX,XBOX, &
                     MABOX,PBOX,NISBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis,jbox)

    end if

    goto 320



  ! CORE INTERSECTION:

  else

    ! Ray tracing in the mantle:

    X = Y2

    call prosca(kekf,KEIN,KFIN)
    lenabs=len2*abs(kekf)
    call absorb(FIN,qabs,lenabs,abscf1)

    if (FIN(1,1) .le. wray*Fstop) then
      qin=qin+FIN(1,1)
      goto 320
    end if

    call incide(FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                FCORE,KECORE,KFCORE,HLCORE,HRCORE,MACORE, &
                N2,m1,m2,totref)

    pin=pin+1
    pout=pout+1

    if (totref .eq. 0) then

      jboxin=jboxin+1
      if (jboxin .gt. 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

    else if (totref .eq. 1) then

      jboxin=jboxin+1
      if (jboxin .gt. 499) then
        qbox=qbox+FIN(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

      goto 320

    else

      goto 410

    end if

    ! Ray tracing in the core:

410 if (FCORE(1,1).le.wray*Fstop .or. pin.ge.pinmax) then
      qin=qin+FCORE(1,1)
      goto 320
    end if

    Y2 = X
    call istri(Y2,KECORE,N2,XN2,NT,IT,nk2,len3,ntri,nis2,-1)

    if (nis2 .eq. 0) then
      qis=qis+FCORE(1,1)
      goto 320
    end if

    ! Radiative transfer in the core:

    if (rtflg2 .eq. 1) then

      call random_number(ran2)
      lenrn=-lenrt2*log(ran2)
      if (lenrn.lt.len3) then
        pin=pin+1
        pout=pout+1
        call rtmc(FCORE,KECORE,KFCORE,HLCORE,HRCORE, &
        X,XP2,YP2,YP22,CSRN2, &
        qabs,lenrn,abscf2,omg2,np,nrn)
        goto 410
      end if
    end if

    ! Ray tracing in the core continued.

    X = Y2

    call prosca(kekf,KECORE,KFCORE)
    lenabs=len3*abs(kekf)
    call absorb(FCORE,qabs,lenabs,abscf2)

    if (FCORE(1,1).le.wray*Fstop) then
      qin=qin+FCORE(1,1)
      goto 320
    end if

    call incide(FCORE,KECORE,KFCORE,HLCORE,HRCORE,MACORE, &
                FIN,KEIN,KFIN,HLIN,HRIN,MAIN, &
                N2,m2,m1,totref)

    pin=pin+1
    pout=pout+1

    if (totref .eq. 0) then

      jboxin=jboxin+1
      if (jboxin .gt. 499) then
        qbox=qbox+FCORE(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FIN,KEIN,KFIN, &
                     HLIN,HRIN,X,MAIN,pin,nis2,jboxin)

      goto 410

    else if (totref .eq. 1) then

      goto 410

    else

      jboxin=jboxin+1
      if (jboxin .gt. 499) then
        qbox=qbox+FCORE(1,1)
        goto 320
      end if

      call rayputd1c(FINBOX,KEINBOX,KFINBOX,HLINBOX,HRINBOX,XINBOX, &
                     MAINBOX,PINBOX,NISINBOX,FCORE,KECORE,KFCORE, &
                     HLCORE,HRCORE,X,MAIN,pin,nis2,jboxin)

      goto 320

    end if
  end if


end program mainGo
