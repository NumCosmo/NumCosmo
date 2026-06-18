program validate_comm_lowl_mod
   use comm_lowl_mod_dist

   integer,      parameter :: i4b = selected_int_kind(9)
   integer,      parameter :: sp  = selected_real_kind(5,30)
   integer,      parameter :: dp  = selected_real_kind(12,200)
   integer,      parameter :: lgt = kind(.true.)

   character(len=512)   :: paramfile, clfile1, clfile2
   real(dp), allocatable, dimension(:, :) :: cl1, cl2
   real(dp)     :: lnL1, lnL2, lnL1ex, lnL2ex
   integer(i4b) :: handle

   if (iargc() /= 3) then
      stop "Usage: validate_comm_lowl_mod paramfile clfile1 clfile2"
   end if
   call getarg(1, paramfile)

   call getarg(2, clfile1)
   call getarg(3, clfile2)

   call comm_lowl_initialize_object(paramfile, handle)

   call read_fiducial_spectrum(clfile1, cl1)
   call read_fiducial_spectrum(clfile2, cl2)

   lnL1ex = -4915.03736417d0
   lnL2ex = -4914.99669525d0
   lnL1 = comm_lowl_compute_lnL(cl1)
   lnL2 = comm_lowl_compute_lnL(cl2)

   write(*, *) '*********************************************************************'
   write(*, *) 'Test no.    Computed lnL       Expected lnL       Relative difference'
   write(*, fmt='(a, I1,f20.8, f20.8, f20.8)') '    ', 1, lnL1, lnL1ex, abs((lnL1-lnL1ex)/lnL1ex)
   write(*, fmt='(a, I1,f20.8, f20.8, f20.8)') '    ', 2, lnL2, lnL2ex, abs((lnL2-lnL2ex)/lnL2ex)
   write(*, *) '*********************************************************************'

end program
