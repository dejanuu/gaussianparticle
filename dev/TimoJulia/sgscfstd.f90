
subroutine sgscfstd(SCFSTD,CSCF,beta,lmax)

    ! Generates the standard deviations for the spherical harmonics
    ! coefficients of the logarithmic radial distance. Version 2002-12-16.


    implicit none

    integer, parameter ::                              &
    sp = kind(1.0),                                  &
    dp = selected_real_kind(2*precision(1.0_sp)),    &
    qp = selected_real_kind(2*precision(1.0_dp))
    integer  :: l,lmax,m
    real (kind=dp) :: FACTRL,beta
    real (kind=dp), dimension(0:256,0:256) :: SCFSTD
    real (kind=dp), dimension(0:256) :: CSCF



    SCFSTD(0,0)=beta*sqrt(CSCF(0))

    do l=1,lmax
        SCFSTD(l,0)=beta*sqrt(CSCF(l))
        do m=1,l
            SCFSTD(l,m)=SCFSTD(l,0)*sqrt(2.0_dp*FACTRL(l-m)/FACTRL(l+m))
        end do
    end do


!    nu=0
!    if (lmin.eq.0) then
!        SCFSTD(0,0)=beta*sqrt(CSCF(0))
!        nu=nu+1
!        do l=1,lmax
!            SCFSTD(l,0)=beta*sqrt(CSCF(l))
!            nu=nu+1
!            do m=1,l
!                SCFSTD(l,m)=SCFSTD(l,0)*sqrt(2.0_dp*FACTRL(l-m)/FACTRL(l+m))
!                nu=nu+2
!            end do
!        end do
!    else
!        do l=lmin,lmax
!            SCFSTD(l,0)=beta*sqrt(CSCF(l))
!            nu=nu+1
!            do m=1,l
!                SCFSTD(l,m)=SCFSTD(l,0)*sqrt(2.0_dp*FACTRL(l-m)/FACTRL(l+m))
!                nu=nu+2
!            end do
!        end do
!    endif

end subroutine sgscfstd