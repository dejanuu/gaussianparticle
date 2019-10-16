program test_bench
use pcgrng
use ginterface
use parameters 
implicit none

real(kind=dp) :: val
integer :: seed, i_stream,nstates
integer :: j1
real(kind=dp) :: rmax,volume


seed = 2
i_stream = 1
nstates = 10**8


call init_rng2(seed,i_stream,nstates)


call init_ginterface("ttest1.in")

do j1=1,5
    call generate_gsphere(rmax)
    call compute_volume(volume)
    write(6,*) "max radius",rmax,"volume",volume
    call save_gsphere()
enddo






end program
