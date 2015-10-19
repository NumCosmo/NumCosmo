!-----------------------------------------------------------------------------
!
!  Copyright (C) 1997-2013 Krzysztof M. Gorski, Eric Hivon,
!                          Benjamin D. Wandelt, Anthony J. Banday, 
!                          Matthias Bartelmann, Hans K. Eriksen, 
!                          Frode K. Hansen, Martin Reinecke
!
!
!  This file is part of HEALPix.
!
!  HEALPix is free software; you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation; either version 2 of the License, or
!  (at your option) any later version.
!
!  HEALPix is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with HEALPix; if not, write to the Free Software
!  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
!
!  For more information about HEALPix see http://healpix.jpl.nasa.gov
!
!-----------------------------------------------------------------------------

!module fitstools
module fitstools_smw
!lpl renamed to avoid possible conflicts

 
  use healpix_types
  implicit none

  real(kind=SP),     private, parameter :: s_bad_value = HPX_SBADVAL
  real(kind=DP),     private, parameter :: d_bad_value = HPX_DBADVAL
  integer(kind=I4B), private, parameter :: i_bad_value = -1637500000
  integer(I4B) ,     private, parameter :: nchunk_max  = 12000


  interface fatal_error
        module procedure fatal_error_womsg, fatal_error_msg
  end interface

!  interface read_bintab
!     module procedure read_bintab8_s, read_bintab8_d, read_bintab4_s, read_bintab4_d
!  end interface

  interface fits2cl
     module procedure fits2cl_s, fits2cl_d
  end interface

  interface map_bad_pixels
     module procedure map_bad_pixels_s, map_bad_pixels_d
  end interface

  private

  public :: fits2cl ,read_dbintab ,getsize_fits ,printerror

contains


  ! define routine with SP I/O
!#include "fits_s_inc.f90"

  ! define routine with DP I/O
!#include "fits_d_inc.f90"

  !=======================================================================
  ! FITS2CL(filename, clin, lmax, ncl, header, units, fmissval)
  !     Read C_l from a FITS file
  !   Aug 2000 : modification by EH of the number of columns actually read
  !
  !     This routine is used for reading input power spectra for synfast
  !
  !  Dec 2004: overloading for single and double precision output array (clin)
  ! Feb 2013: added fmissval (value to be given to NaN or Inf C(l))
  !      if absent, those C(l) are left unchanged
  !   
  !=======================================================================
  subroutine fits2cl_s(filename, clin, lmax, ncl, header, units, fmissval)
    !=======================================================================
    !        single precision
    !=======================================================================
    integer(I4b), parameter :: KCL = SP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),        DIMENSION(0:lmax,1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),            INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    real(KCL),                       optional,  intent(IN):: fmissval

    real(DP), dimension(:,:), allocatable :: cl_dp
    real(DP) :: fmvd

    ! since the arrays involved are small
    ! read in double precision, and then cast in single
    allocate(cl_dp(0:lmax, 1:ncl))
    if (present(fmissval)) then
       fmvd = fmissval * 1.0_DP
       call fits2cl_d(filename, cl_dp, lmax, ncl, header, units, fmvd)
    else
       call fits2cl_d(filename, cl_dp, lmax, ncl, header, units)
    endif
    clin(0:lmax, 1:ncl) = cl_dp(0:lmax, 1:ncl)

    return
  end subroutine fits2cl_s

  subroutine fits2cl_d(filename, clin, lmax, ncl, header, units, fmissval)
    !=======================================================================
    !        double precision
    !=======================================================================
    integer(I4b), parameter :: KCL = DP
    !
    CHARACTER(LEN=*),                          INTENT(IN) :: filename
    INTEGER(I4B),                              INTENT(IN) :: lmax, ncl
    REAL(KCL),         DIMENSION(0:lmax, 1:ncl), INTENT(OUT) :: clin
    CHARACTER(LEN=*), DIMENSION(1:),   INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units
    real(KCL),                       optional,  intent(IN):: fmissval

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),ncl_file, naxis
    INTEGER(I4B) :: firstpix,lmax_file,lmax_min,nelems,datacode,repeat,width
    CHARACTER(LEN=80) :: comment, pdmstr
    LOGICAL :: extend
    INTEGER(I4B) :: nmove, hdutype, hdunum
    INTEGER(I4B) :: column, frow
    REAL(KCL), DIMENSION(:), ALLOCATABLE :: clin_file

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: rowlen, nrows, varidat
    INTEGER(I4B),      dimension(1:maxdim) :: tbcol
    CHARACTER(LEN=20), dimension(1:maxdim) :: ttype, tform, tunit
    CHARACTER(LEN=20)                      :: extname
    LOGICAL :: anynull
    REAL(KCL) ::  nullval
    logical :: planck_format
    integer(i8b) :: nbads


    !-----------------------------------------------------------------------
    status=0
    anynull = .false.
    ttype=''
    tform=''
    tunit=''
    extname=''
    comment=''

    unit = 110
    naxes(1) = 1
    naxes(2) = 1

    readwrite=0
    call ftnopn(unit,filename,readwrite, status) !open primary HDU or specified HDU
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftghdn(unit, hdunum)
    if (hdunum == 1) then  ! in primary HDU: move to next HDU
       !     determines the presence of image
       call ftgkyj(unit,'NAXIS', naxis, comment, status)
       if (status > 0) call printerror(status)

       !     determines the presence of an extension
       call ftgkyl(unit,'EXTEND', extend, comment, status)
       if (status > 0) call printerror(status)

       nmove = +1
       call ftmrhd(unit, nmove, hdutype, status)
    else ! already in non primary HDU
       extend = .true.
       call ftghdt(unit, hdutype, status)
    endif

    if (extend) then ! there is an extension

       call assert ((hdutype==1).or.(hdutype==2), 'this is not a table')

       !        reads keywords related to table layout
       if (hdutype==1) then ! ASCII table
         call ftghtb(unit, maxdim, rowlen, &
              &      nrows, ncl_file, ttype, tbcol, &
              &      tform, tunit, extname, status)
         repeat = 1
       else ! binary table
         call ftghbn(unit, maxdim, &
              &      nrows, ncl_file, ttype,        &
              &      tform, tunit,extname, varidat, status)
          call ftbnfm(tform(1), datacode, repeat, width, status)
       endif

       status = 0
       header = ""
       call get_clean_header(unit, header, filename, status)

       if (present(units)) then 
          do column = 1, min(ncl, ncl_file)
             units(column) = adjustl(tunit(column))
          enddo
       endif

       !        reads the columns
       column = 1
       frow = 1
       firstpix = 1
       lmax_file = nrows*repeat - 1
       lmax_min = MIN(lmax,lmax_file)
       nelems = lmax_min + 1
       nullval = 0.0_KCL ! CFITSIO will leave bad pixels unchanged, and so will this routine
       if (present(fmissval)) then
          if (fmissval == 0.0_KCL) then
             nullval = HPX_DBADVAL ! sentinel value to be given to bad pixels by CFITSIO, will later be mapped to user defined 0
          else
             nullval = fmissval ! user defined value to be given to bad pixels by CFITSIO, and by this routine
          endif
       endif

! check for the special Planck format (i.e. one additional column)
       planck_format=.true.
       call ftgkys(unit,'PDMTYPE',pdmstr,comment,status)
       if (status==202) then
         planck_format=.false.
         status=0
       endif
       allocate(clin_file(0:lmax_min),stat=status)
       clin = 0.0_KCL                         ! modification by EH
       if (planck_format) then
         do column = 1, MIN(ncl, ncl_file-1) ! modification by EH
            clin_file(:) = 0.0_KCL
            call ftgcvd(unit, column+1_i4b, frow, firstpix, nelems, nullval, &
                 &        clin_file(0), anynull, status)
            clin(0:lmax_min, column) = clin_file(0:lmax_min)
         enddo
       else
         do column = 1, MIN(ncl, ncl_file) ! modification by EH
            clin_file(:) = 0.0_KCL
            call ftgcvd(unit, column, frow, firstpix, nelems, nullval, &
                 &        clin_file(0), anynull, status)
            clin(0:lmax_min, column) = clin_file(0:lmax_min)
         enddo
       endif
       deallocate(clin_file)
       if (present(fmissval)) then
          if (nullval /= fmissval) call map_bad_pixels(clin, nullval, fmissval, nbads)
       endif
    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension in '//trim(filename))
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine fits2cl_d

  !=======================================================================
  !=======================================================================
  subroutine read_dbintab(filename,map,npixtot,nmaps,nullval,anynull, units)
    !=======================================================================
    !     Read a FITS file
    !
    !     slightly modified to deal with vector column 
    !     in binary table       EH/IAP/Jan-98
    !
    !     Reads a double-precision array with precomputed plms, used by syn/anafast
    !                FKH/Apr-99
    !=======================================================================
    CHARACTER(LEN=*),                           INTENT(IN) :: filename
    INTEGER(I4B),                               INTENT(IN) :: npixtot, nmaps
    !       REAL(DP), DIMENSION(0:npixtot-1,1:nmaps),   INTENT(OUT) :: map
    REAL(DP), DIMENSION(0:,1:),                 INTENT(OUT) :: map
    REAL(DP),                                   INTENT(OUT) :: nullval
    LOGICAL(LGT),                               INTENT(OUT) ::  anynull
    !       CHARACTER(LEN=*), dimension(1:), optional,           INTENT(OUT) :: header
    CHARACTER(LEN=*), dimension(1:), optional,  INTENT(OUT) :: units

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B) :: group,firstpix,npix,i
    REAL(SP) :: blank, testval
    REAL(DP) ::  bscale,bzero
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT) :: extend
    INTEGER(I4B) :: nmove, hdutype
    INTEGER(I4B) :: column, frow, imap
    INTEGER(I4B) :: datacode, repeat, width

    INTEGER(I4B), PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B) :: nrows, tfields, varidat
    CHARACTER(LEN=20), dimension(1:maxdim) :: ttype, tform, tunit
    CHARACTER(LEN=20)                      :: extname
    !-----------------------------------------------------------------------
    status=0

    unit = 147
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    anynull = .false.
    bscale = 1.0d0
    bzero = 0.0d0
    blank = -2.e25
    nullval = bscale*blank + bzero

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------

    !     determines the presence of image
    call ftgkyj(unit,'NAXIS', naxis, comment, status)
    if (status > 0) call printerror(status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    if (status > 0) status = 0 ! no extension : 
    !     to be compatible with first version of the code

    if (naxis > 0) then ! there is an image
       !        determine the size of the image (look naxis1 and naxis2)
       call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)

       !        check that it found only NAXIS1
       if (nfound == 2 .and. naxes(2) > 1) then
          print *,'multi-dimensional image'
          print *,'expected 1-D data.'
          call fatal_error
       end if

       if (nfound < 1) then
          call printerror(status)
          print *,'can not find NAXIS1.'
          call fatal_error
       endif

       npix=naxes(1)
       if (npix /= npixtot) then
          print *,'found ',npix,' plms'
          print *,'expected ',npixtot
          call fatal_error
       endif

       call ftgkyd(unit,'BSCALE',bscale,comment,status)
       if (status == 202) then ! BSCALE not found
          bscale = 1.0d0
          status = 0
       endif
       call ftgkyd(unit,'BZERO', bzero, comment,status)
       if (status == 202) then ! BZERO not found
          bzero = 0.0d0
          status = 0
       endif
       call ftgkye(unit,'BLANK', blank, comment,status)
       if (status == 202) then ! BLANK not found 
          ! (according to fitsio BLANK is integer)
          blank = -2.e25
          status = 0
       endif
       nullval = bscale*blank + bzero

       !        -----------------------------------------

       group=1
       firstpix = 1
       call ftgpvd(unit,group,firstpix,npix,nullval,map,anynull,status)
       ! if there are any NaN pixels, (real data)
       ! or BLANK pixels (integer data) they will take nullval value
       ! and anynull will switch to .true.
       ! otherwise, switch it by hand if necessary
       testval = 1.e-6 * ABS(nullval)
       do i=0, npix-1
          if (ABS(map(i,1)-nullval) < testval) then
             anynull = .true.
             goto 111
          endif
       enddo
111    continue

    else if (extend) then ! there is an extension
       nmove = +1
       call ftmrhd(unit, nmove, hdutype, status)
       !cc         write(*,*) hdutype

       call assert (hdutype==2, 'this is not a binary table')

       !        reads all the keywords
       call ftghbn(unit, maxdim, &
            &        nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &        status)

       if (tfields < nmaps) then
          print *,'found ',tfields,' maps in the file'
          print *,'expected ',nmaps
          call fatal_error
       endif

       !        finds the bad data value
       call ftgkyd(unit,'BAD_DATA',nullval,comment,status)
       if (status == 202) then ! bad_data not found
          nullval = d_bad_value
          status = 0
       endif

       !          if (nlheader > 0) then
       !             call get_clean_header(unit, header, filename, status)
       !          endif

       if (present(units)) then
          units = 'unknown'
          do imap = 1, nmaps
             units(imap) = adjustl(tunit(imap))
          enddo
       endif

       do imap = 1, nmaps
          !parse TFORM keyword to find out the length of the column vector
          call ftbnfm(tform(imap), datacode, repeat, width, status)

          !reads the columns
          column = imap
          frow = 1
          firstpix = 1
          npix = nrows * repeat
          if (npix /= npixtot) then
             print *,'found ',npix,' plms'
             print *,'expected ',npixtot
             call fatal_error("read_dbintab "//trim(filename))
          endif
          call ftgcvd(unit, column, frow, firstpix, npix, nullval, &
               &        map(0:npix-1,imap), anynull, status)
       enddo

    else ! no image no extension, you are dead, man
       call fatal_error(' No image, no extension')
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    return
  end subroutine read_dbintab


  !=======================================================================
  subroutine printerror(status)
    !=======================================================================
    !     Print out the FITSIO error messages to the user
    !=======================================================================
    INTEGER(I4B), INTENT(IN) :: status
    CHARACTER ::  errtext*30,errmessage*80
    !-----------------------------------------------------------------------
    !     check if status is OK (no error); if so, simply return
    if (status .le. 0)return

    !     get the text string which describes the error
    call ftgerr(status,errtext)
    print *,'FITSIO Error Status =',status,': ',errtext

    !     read and print out all the error messages on the FITSIO stack
    call ftgmsg(errmessage)
    do while (errmessage /= ' ')
       print *,errmessage
       call ftgmsg(errmessage)
    end do

    return
  end subroutine printerror


  !=======================================================================
  subroutine read_par(filename,nside,lmax,tfields,mmax)
    !=======================================================================
    !        Read nside, lmax, tfields and mmax from a FITS file          
    !    parameters not found take a value of -1
    !
    !         Frode K. Hansen, April 1999
    !         EH, Dec 2004
    !
    !=======================================================================
    CHARACTER(LEN=*), INTENT(IN)            :: filename

    INTEGER(I4B), INTENT(OUT)   :: nside, lmax, tfields
    integer(i4b), intent(out), optional :: mmax

    INTEGER(I4B) :: status,unit,readwrite,blocksize,naxis
    CHARACTER(LEN=80) :: comment, ttype1
    LOGICAL(LGT) ::  extend, anyf
    INTEGER(I4B)::  nmove, hdutype, idmax, nrows

    !-----------------------------------------------------------------------
    status=0
    unit = 120
    comment=''
    ttype1=''

    readwrite=0
    call ftopen(unit,filename,readwrite,blocksize,status)
    if (status > 0) call printerror(status)
    !     -----------------------------------------
    call ftgkyj(unit,'NAXIS', naxis, comment, status)

    !     determines the presence of an extension
    call ftgkyl(unit,'EXTEND', extend, comment, status)
    call assert (status<=0, 'No Extension in FITS file!')

    nmove = +1
    call ftmrhd(unit, nmove, hdutype, status)

    call assert (hdutype==2, 'This is not a FITS binary-table')

    call ftgkyj(unit,'NSIDE',nside,comment,status)
    if (status == 202) then
       print*,'WARNING: NSIDE keyword not found!'
       nside = -1
       status = 0
    endif

    call ftgkyj(unit,'TFIELDS',tfields,comment,status)
    if (status == 202) then
       print*,'WARNING: TFIELDS keyword not found!'
       tfields = -1
       status = 0
    endif
    
    call ftgkyj(unit,'MAX-LPOL',lmax,comment,status)
    if (status == 202) then
       status = 0
       ! if not found, determines if file contains indexed list of alm
       ! if so, find out largest l
       if (tfields >= 3 .and. hdutype==2) then ! 3 column binary table
          call ftgkys(unit,'TTYPE1',ttype1,comment,status)
          ttype1 = trim(strupcase(adjustl(ttype1)))
          if (trim(ttype1(1:5)) == 'INDEX') then
             call ftgkyj(unit, 'NAXIS2', nrows, comment, status) ! find number of rows
             call ftgcvj(unit, 1_i4b, nrows, 1_i4b, 1_i4b, 0_i4b, idmax, anyf, status) ! read element on last row of first column
             if (status == 0) then
                lmax = int(sqrt(   real(idmax-1, kind = DP)  ) )
                if (lmax > 0) goto 1000
             endif
          endif
       endif
       print*,'WARNING: MAX-LPOL keyword not found!'
       lmax = -1
       status = 0
    endif
1000 continue

    if (present(mmax)) then
       call ftgkyj(unit,'MAX-MPOL',mmax,comment,status)
       if (status == 202) then
          print*,'WARNING: MAX-MPOL keyword not found!'
          mmax = -1
          status = 0
       endif
    endif
    call ftclos(unit, status)

  end subroutine read_par

  !=======================================================================
  function getnumext_fits(filename)
    !=======================================================================
    !  result = getnumext_fits(filename)
    !    returns the number of extensions present in FITS file 'filename'
    !
    ! EH, Nov 2004
    ! April 2007: close file on exit
    !=======================================================================
    character(LEN=*), intent(IN)             :: filename
    integer(i4b)                             :: getnumext_fits
    !
    integer(i4b) :: status, unit, readwrite, blocksize, nhdu
    !-----------------------------------------------------------------------
    status         = 0
    unit           = 149
    getnumext_fits = 0
    readwrite      = 0 ! Read only
    call ftopen(unit, filename, readwrite, blocksize, status)
    if (status > 0) then
       call printerror(status)
       call ftclos(unit, status)
       return
    endif

    call ftthdu(unit, nhdu, status)
    getnumext_fits = nhdu - 1

    call ftclos(unit, status)
    return
  end function getnumext_fits
  !=======================================================================
  function getsize_fits(filename, nmaps, ordering, obs_npix, &
       &    nside, mlpol, type, polarisation, fwhm_arcmin, beam_leg, & 
       &    coordsys, polcconv, extno)
    !=======================================================================
    !  result = getsize_fits(filename, nmaps, ordering, nside, mlpol, type)
    !     Get the number of pixels stored in a map FITS file.
    !     Each pixel is counted only once 
    !     (even if several information is stored on each of them, see nmaps).
    !     Depending on the data storage format, this may be :
    !       - equal or smaller to the number Npix of Healpix pixels available 
    !          over the sky for the given resolution (Npix = 12*nside*nside)
    !       - equal or larger to the number of non blank pixels (obs_npix)
    !     
    !     filename = (IN) name of the FITS file containing the map
    !     nmaps = (OPTIONAL, OUT) number of maps in the file
    !     ordering = (OPTIONAL, OUT) Healpix ordering scheme used
    !                  0 = unknown
    !                  1 = RING
    !                  2 = NESTED
    !     obs_npix =  (OPTIONAL, OUT) number of non blanck pixels.    
    !                It is set to -1 if it can not be determined from header
    !                information alone
    !     nside = (OPTIONAL, OUT) Healpix parameter Nside
    !                 returns a negative value if not found
    !     mlpol = (OPTIONAL, OUT) maximum multipole used to generate the map (for simulated map)
    !                 returns a negative value if not found
    !     type = (OPTIONAL, OUT) Healpix/FITS file type
    !                  <0 : file not found, or not valid
    !                  0  : image only fits file, deprecated Healpix format
    !                        (result = 12 * nside * nside)
    !                  1  : ascii table, generally used for C(l) storage
    !                  2  : binary table : with implicit pixel indexing (full sky)
    !                        (result = 12 * nside * nside)
    !                  3  : binary table : with explicit pixel indexing (generally cut sky)
    !                        (result <= 12 * nside * nside)
    !                999  : unable to determine the type
    !     polarisation = (OPTIONAL, OUT) presence of polarisation data in the file
    !                  <0 : can not find out
    !                   0 : no polarisation
    !                   1 : contains polarisation (Q,U or G,C)
    !     fwhm_arcmin     = (OPTIONAL, DP, OUT) returns the beam FWHM read from FITS header, 
    !                        translated from Deg (hopefully) to arcmin
    !                     returns a negative value if not found
    !
    !     beam_leg     = (OPTIONAL, CHR, OUT) filename of beam or filtering window function applied to data
    !                     returns a empty string if not found
    !     coordsys     = (OPTIONAL, CHR, OUT) string describing coordinate system,  
    !                     G = Galactic, E = ecliptic, C = celestial = equatorial.
    !                     empty if not found.
    !
    !     polcconv     = (OPTIONAL, I4B, OUT) coordinate convention for polarisation
    !                   0: unknown
    !                   1: COSMO (default for Healpix)
    !                   2: IAU
    !
    !     extno = (OPTIONAL, IN) specify FITS extension to look at (0 based)
    !                  default = 0 (first extension)
    !
    !     Benjamin D. Wandelt January 1998
    !     includes Eric Hivon's modification of FITS columns
    !     and draws heavily on the read_bintab routine.
    !
    !     addition of optional Ordering by E.H. (July 98)
    !     addition of optional Nside and Mlpol by ?? (??)
    !     addition of optional type by E.H. (Sept 00)
    !     improved for ASCII table E.H (Feb 03)
    !     addition of extno, E.H (Nov 04)
    !     addition of fwhm_arcmin and beam_leg, EH (Jan 05).
    !     addition of polcconv, EH (June 05).
    !=======================================================================
    character(LEN=*), intent(IN)             :: filename
    integer(I4B),     intent(out),  optional :: nmaps
    integer(I4B),     intent(out),  optional :: ordering
    integer(I4B),     intent(out),  optional :: nside
    integer(I4B),     intent(out),  optional :: mlpol
    integer(I4B),     intent(out),  optional :: obs_npix
    integer(I4B),     intent(out),  optional :: type
    integer(I4B),     intent(out),  optional :: polarisation
    real(DP),         intent(out),  optional :: fwhm_arcmin
    character(LEN=*), intent(out),  optional :: beam_leg
    character(LEN=*), intent(out),  optional :: coordsys
    integer(I4B),     intent(out),  optional :: polcconv
    integer(I4B),     intent(in),   optional :: extno

    INTEGER(I4B)           :: nmaps_in, ordering_in, nside_in
    INTEGER(I4B)           :: mlpol_in, obs_npix_in, ftype_in
    INTEGER(I4B)           :: extno_in, polcconv_in
    real(DP)               :: fwhm_arcmin_in
    character(len=FILENAMELEN)        :: beam_leg_in
    character(len=20)        :: coordsys_in
    INTEGER(I4B)           :: grain
    CHARACTER(len=20)      :: indxschm
    CHARACTER(LEN=20)      :: order_val, object_val, ttype_val, polcconv_val
!     INTEGER(I4B)           :: getsize_fits
    INTEGER(I8B)           :: getsize_fits
    LOGICAL(LGT)           ::  polar_in
    character(len=3),  dimension(1:10,1:2)  :: defpol
    logical(kind=LGT), dimension(1:2)       :: pf
    integer(kind=I4B)                       :: ndp, j, k

    INTEGER(I4B)      :: status,unit,readwrite,blocksize,naxes(2),nfound, naxis
    INTEGER(I4B)      :: i
!     INTEGER(I4B) :: nsize
    INTEGER(I8B)      :: nsize
    CHARACTER(LEN=80) :: comment
    LOGICAL(LGT)      ::  extend
    INTEGER(I4B)      ::  nmove, hdutype, hdunum
    INTEGER(I4B)      :: datacode, repeat1, repeat2, width

    INTEGER(I4B),           PARAMETER :: maxdim = 40 !number of columns in the extension
    INTEGER(I4B)                      :: nrows, tfields, varidat, rowlen
    CHARACTER(LEN=20), dimension(1:maxdim) :: ttype, tform, tunit
    INTEGER(I4B),      dimension(1:maxdim) :: tbcol
    CHARACTER(LEN=20)                      :: extname
    !-----------------------------------------------------------------------
    status=0
    order_val = ''
    nmaps_in = 0
    ordering_in = 0
    mlpol_in = -1
    nside_in = -1
    ftype_in = 999 !   
    obs_npix_in = -1
    unit = 119
    naxes(1) = 1
    naxes(2) = 1
    nfound = -1
    extno_in = 0
    polcconv_in = 0
    comment=''
    ttype=''
    tform=''
    tunit=''
    extname=''
    if (present(extno)) extno_in = extno
    naxis = 0
    extend = .false.

    readwrite=0
    !call ftopen(unit,filename,readwrite,blocksize,status)
    call ftnopn(unit,filename,readwrite,status)
    if (status > 0) then
       ftype_in = -1
       getsize_fits = -1
       call printerror(status)
       call ftclos(unit, status)
       return
    endif
    !     -----------------------------------------
    call ftghdn(unit, hdunum)
    if (hdunum == 1) then  ! in primary HDU: move to next HDU
       !     determines the presence of image
       call ftgkyj(unit,'NAXIS', naxis, comment, status)

       !     determines the presence of an extension
       call ftgkyl(unit,'EXTEND', extend, comment, status)
       if (status > 0) then
          ftype_in = 0
          status = 0 ! no extension : 
       !     to be compatible with first version of the code
       endif

    else ! already in non primary HDU
       extend = .true.
       call ftghdt(unit, hdutype, status)
    endif

    if (naxis > 0) then 
       !---------------------------------
       ! there is an image
       !---------------------------------
       !        determine the size of the image (look naxis1 and naxis2)
       call ftgknj(unit,'NAXIS',1_i4b,2_i4b,naxes,nfound,status)

       !        check that it found only NAXIS1
       if (nfound == 2 .and. naxes(2) > 1) then
          print *,'multi-dimensional image'
          print *,'expected 1-D data.'
          call ftclos(unit, status)
          call fatal_error
       end if

       if (nfound < 1) then
          call printerror(status)
          print *,'can not find NAXIS1.'
          call ftclos(unit, status)
          call fatal_error
       endif

       nsize=naxes(1)

       call ftgkys(unit,'ORDERING',order_val,comment,status)
       if (status == 202) then ! Ordering not found
          ordering_in = 0
          order_val = ''
          status = 0
       endif

       call ftgkyj(unit,'NSIDE',nside_in,comment,status)
       if (status == 202) then ! Nside not found
          nside_in = -1
          status = 0
       endif

       if (present(mlpol)) then
          call ftgkyj(unit,'MAX-LPOL',mlpol_in,comment,status)
          if (status == 202) then ! max-lpol not found
             mlpol_in = -1
             status = 0
          endif
       endif

       if (present(polarisation)) then
          polarisation = 0
       endif
    else if (extend) then 
       !-----------------------------------------------------------------
       ! there is an extension
       !-----------------------------------------------------------------

       nmove =  extno_in
       if (hdunum == 1) nmove = extno_in + 1
       call ftmrhd(unit, nmove, hdutype, status)
       if (status > 0) then ! extension not found
          print*,'Extension #',extno_in,' not found in '//trim(filename)
          call printerror(status)
          call ftclos(unit, status)
          call fatal_error
       endif
       !c         write(*,*) hdutype

       !        reads all the keywords
       if (hdutype == 2) then ! binary table
          call ftghbn(unit, maxdim, &
            &         nrows, tfields, ttype, tform, tunit, extname, varidat, &
            &         status)
       else ! ASCII table (hdutype = 1)
          ftype_in = 1
          call ftghtb(unit, maxdim, &
            &         rowlen, nrows, tfields, ttype, tbcol, tform, tunit, &
            &         extname, status)
       endif

       !        parse TFORM keyword to find out the length of the column vector
       repeat1 = 1
       repeat2 = 1
       call ftbnfm(tform(1), datacode, repeat1, width, status)
       if (tfields > 1) call ftbnfm(tform(2), datacode, repeat2, width, status)

       nsize = int(nrows, kind=i8b) * max(repeat1,repeat2) ! corrected Oct-03

       nmaps_in = tfields

       call ftgkys(unit,'ORDERING',order_val,comment,status)
       if (status == 202) then ! Ordering not found
          ordering_in = 0
          order_val = ''
          status = 0
       endif

       call ftgkyj(unit,'NSIDE',nside_in,comment,status)
       if (status == 202) then ! Nside not found
          nside_in = -1
          status = 0
       endif

       call ftgkyj(unit,'OBS_NPIX',obs_npix_in,comment,status)
       if (status == 202) then ! obs_npix not found
          obs_npix_in = -1
          status = 0
       endif

       if (present(mlpol)) then
          call ftgkyj(unit,'MAX-LPOL',mlpol_in,comment,status)
          if (status == 202) then ! max-lpol not found
             status = 0
             call ftgkyj(unit,'LMAX',mlpol_in,comment,status)
             if (status == 202) then
                mlpol_in = -1
                status = 0
             endif
          endif
       endif

       if (present(fwhm_arcmin)) then
          call ftgkyd(unit,'FWHM',fwhm_arcmin_in, comment, status)
          if (status == 202) then ! fwhm not found
             fwhm_arcmin_in = -1.
             status = 0
          else
             fwhm_arcmin_in = 60.0_dp * fwhm_arcmin_in
          endif
       endif

       if (present(beam_leg)) then
          call ftgkys(unit,'BEAM_LEG',beam_leg_in, comment, status)
          if (status == 202) then ! beam_leg not found
             beam_leg_in = ' '
             status = 0
          endif
       endif

       if (present(coordsys)) then
          call ftgkys(unit,'COORDSYS',coordsys_in, comment, status)
          if (status == 202) then ! coordsys not found
             coordsys_in = ' '
             status = 0
          endif
       endif

       ! determines pixel indexing (for binary tables)
       if (present(type) .and. ftype_in == 999) then
          ! most stringent test
          call ftgkys(unit,'INDXSCHM',indxschm,comment,status)
          if (status == 0) then ! found
             ftype_in = 3
             if (trim(indxschm) == 'IMPLICIT') ftype_in = 2
             goto 1000
          else
             status = 0
          endif
          ! 2nd most stringent test
          call ftgkyj(unit,'GRAIN',grain,comment,status)
          if (status == 0) then ! found
             ftype_in = 3
             if (grain == 0) ftype_in = 2
             goto 1000
          else
             status = 0
          endif
          ! 3rd most stringent test
          if (trim(ttype(1)) /= '') then
             if (trim(ttype(1)) == 'PIXEL') then
                ftype_in = 3
             else
                ftype_in = 2
             endif
             goto 1000
          endif
          ! lousy test
          call ftgkys(unit,'OBJECT',object_val,comment,status)
          if (status == 0) then
             if (trim(object_val) == 'PARTIAL') ftype_in = 3
             if (trim(object_val) == 'FULLSKY') ftype_in = 2
             if (ftype_in /= 999) goto 1000
          else
             status = 0
          endif
          ! very lousy test
          if (nside_in > 0) then
             ftype_in = 3
             if (12*nside_in*nside_in == nsize) ftype_in = 2
             goto 1000
          endif
       endif
1000   continue

       ! find out if polarisation data is present
       if (present(polarisation)) then
          if (tfields < 3) then
             polarisation = 0 ! no polarisation
             goto 2000
          endif
          call ftgkyl(unit,'POLAR',polar_in,comment,status)
          if (status == 0) then ! polar found
             polarisation = 0
             if (polar_in) polarisation = 1
             goto 2000
          else ! polar not found
             status = 0
             polarisation = -1
             if (hdutype <= 0) goto 2000
             if (hdutype == 1) then ! ascii table -> power spectra
                ndp = 4
                defpol(1:ndp,1) = (/ "GRA","E-M","POW","EE " /)
                defpol(1:ndp,2) = (/ "CUR","B-M","POW","BB " /)    
             endif
             if (hdutype == 2) then ! binary table -> maps
                ndp = 4
                defpol(1:ndp,1) = (/ "Q-P","Q_P","Q P","Q  " /)
                defpol(1:ndp,2) = (/ "U-P","U_P","U P","U  " /)
             endif
             pf(:) = .false.
             do i = 2, tfields ! do not consider first field (generally temperature)
                ttype_val = adjustl(ttype(i))
                call ftupch(ttype_val) ! upper case
                do k=1,2
                   do j=1, ndp
                      if (index( ttype_val, trim(defpol(j,k)) ) == 1) then
                         pf(k) = .true.
                         goto 1500 ! move to next field
                      endif
                   enddo
                enddo ! loop on k
1500            continue
             enddo ! loop on i
             if (pf(1) .and. pf(2)) polarisation = 1
             if (.not. (pf(1) .or. pf(2))) polarisation = 0
          endif ! polar not found
       endif ! present(polarisation)
2000   continue

       call ftgkys(unit,'POLCCONV',polcconv_val,comment,status)
       if (status == 0) then
          if (trim(polcconv_val) == 'COSMO') polcconv_in = 1
          if (trim(polcconv_val) == 'IAU')   polcconv_in = 2
       else
          status = 0
       endif

    else ! no image no extension, you are dead, man
       ftype_in = -1
       call fatal_error(' No image, no extension')
    endif

    !     close the file
    call ftclos(unit, status)

    !     check for any error, and if so print out error messages
    if (status > 0) call printerror(status)

    getsize_fits=nsize

    call ftupch(order_val) ! convert order_val to upper case
    if (order_val(1:4) == 'RING') ordering_in = 1
    if (order_val(1:4) == 'NEST') ordering_in = 2

    if (present(nmaps)) nmaps = nmaps_in
    if (present(mlpol)) mlpol = mlpol_in
    if (present(obs_npix)) obs_npix = obs_npix_in
    if (present(ordering)) ordering = ordering_in
    if (present(nside)) nside = nside_in
    if (present(type)) type = ftype_in
    if (present(fwhm_arcmin)) fwhm_arcmin = fwhm_arcmin_in
    if (present(beam_leg)) beam_leg = adjustl(beam_leg_in)
    if (present(coordsys)) coordsys = adjustl(coordsys_in)
    if (present(polcconv)) polcconv = polcconv_in

    return
  end function getsize_fits
  !=======================================================================
 
  !======================================================================
  subroutine putrec(unit, card, status)
    !======================================================================
    ! append, delete or update a card from the current header
    ! (deletion of keyword KWD is described by : - KWD)
    ! EH, version 1.0, Dec 2001
    !======================================================================
    integer(kind=I4B), intent(IN)  :: unit
    character(len=*),  intent(IN)  :: card
    integer(kind=I4B), intent(OUT) :: status

    character(len=80) :: cardfits,record
    character(len=8)  :: kwd
    integer(kind=I4B) :: hdtype
    character(len=80) :: nullarr(0)
    !======================================================================

    status = 0
    cardfits=''
    record=''
    call ftgthd(card, cardfits, hdtype, status)
    kwd = cardfits(1:8)
    status = 0

    select case (hdtype)
    case (-1) ! delete keyword (starting from the top)
       call ftgrec(unit,0_i4b,record,status)
       !             call ftdkey(unit, kwd, status)
       ! patch around cfitsio bug 
       ! (ftdkey does not recognize wild cards while ftgnxk does)
       do
          call ftgnxk(unit, kwd, 1_i4b, nullarr, 0_i4b, record, status)
          if (status /= 0) exit
          call ftdkey(unit, record(1:8), status)
       enddo
    case (0) ! append or update
       if (kwd == 'CONTINUE') then
          call ftprec(unit, trim(card), status)
          call ftplsw(unit, status)
       else
          ! delete keyword in its current location (if any)
          call ftdkey(unit, kwd, status)
          status = 0
          ! append
          call ftprec(unit, cardfits, status)
       endif
    case (1) ! append (for HISTORY and COMMENT)
       call ftprec(unit, cardfits, status)
    case default
       write(unit=*,fmt=*)" Unexpected card format in fits header :"
       write(unit=*,fmt="(a80)") card
       write(unit=*,fmt=*)" card not written."
    end select
    status = 0

    return
  end subroutine putrec

  !====================================================================
  subroutine get_clean_header(unit, header, filename, error, xalso, xonly)
    !====================================================================
    ! get_clean_header(unit, header, error [, xalso, xonly])
    ! gets the FITS header from unit, after excluding some default keywords 
    !  defined in def_excl
    ! if header in non void on input, its content will be concatenated with that
    !  of the FITS header
    ! if xalso is defined as a list of keywords, they are also excluded from the header
    ! if xonly is defined as a list of keywords, only those keywords are excluded from
    ! the header.
    ! xonly and xalso are exclusive
    !====================================================================
    INTEGER(I4B),                    intent(IN)           :: unit
    CHARACTER(LEN=*), DIMENSION(1:), INTENT(IN OUT)       :: header
    CHARACTER(LEN=*),                INTENT(IN)           :: filename
    INTEGER(I4B),                    intent(OUT)          :: error
    character(len=8), dimension(1:), intent(IN), optional :: xalso
    character(len=8), dimension(1:), intent(IN), optional :: xonly

    INTEGER(I4B) :: nlheader, status, i, n_excl
    CHARACTER(LEN=80) :: record
    CHARACTER(len=8), dimension(:), allocatable :: to_excl

    CHARACTER(len=8), dimension(1:21) :: def_excl
    !====================================================================

    ! keywords to be excluded by default from output header
    ! Note that TTYPE# and TUNIT# keywords are not excluded
    ! from the output header as they might be useful to the 
    ! calling routines
    def_excl=(/&
         & "SIMPLE  ","BITPIX  ","NAXIS   ",&
         & "NAXIS#  ","PCOUNT  ","GCOUNT  ",&
         & "EXTEND  ","ORIGIN  ","DATE*   ",&
         & "TFIELDS ","TFORM#  ",           & 
         & "TBCOL#  ","EXTNAME ","CTYPE#  ",&
         & "CRVAL#  ","CRPIX#  ","CDELT#  ",&
         & "XTENSION","INSTRUME","TELESCOP",&
         & "PDMTYPE "/)

    error = 0
    record=''

    if (present(xonly)) then 
       n_excl = size(xonly)
       allocate(to_excl(1:n_excl))
       to_excl = xonly

    else if (present(xalso)) then
       n_excl = size(xalso) + size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl(1:size(def_excl)) = def_excl
       to_excl(size(def_excl)+1:n_excl) = xalso

    else
       n_excl = size(def_excl)
       allocate(to_excl(1:n_excl))
       to_excl = def_excl
    endif

    nlheader=size(header)
    ! go to end of fortran header
    do i = 1, nlheader
       if (trim(header(i)) == "") exit
    enddo
    ! go to top of fits file header
    status=0
    call ftgrec(unit,0_i4b,record,status)
    ! read in all header lines except those excluded
    do
       call ftgnxk(unit,'*',1_i4b,to_excl,n_excl,record,status)
       if (status > 0) exit ! end of header
       if (i > nlheader) then
          write(unit=*,fmt="(a,i5,a)") &
               & " WARNING : The header in "//  &
               &    trim(filename)//" has more than ", &
               &  nlheader," lines."
          print*," It will be truncated."
          error = 1
          exit
       endif
       header(i)=record
       i=i+1
    enddo
    status=0

    return
  end subroutine get_clean_header
  !====================================================================

  
  subroutine fatal_error_msg (msg)
    character(len=*), intent(in) :: msg
       print *,'Fatal error: ', trim(msg)
    call exit_with_status(1)
  end subroutine fatal_error_msg

  subroutine fatal_error_womsg
      print *,'Fatal error'
    call exit_with_status(1)
  end subroutine fatal_error_womsg


  !-----------------------------------------------------
  subroutine assert (testval,msg,errcode)
    logical, intent(in) :: testval
    character(len=*), intent(in), optional :: msg
    integer(i4b), intent(in), optional :: errcode

    if (testval) return

    print *,"Assertion failed: "
    if (present(msg)) print *, trim(msg)
    if (present(errcode)) call exit_with_status (errcode)
    call exit_with_status(1)
  end subroutine assert

  function strupcase(instr) result(outstr)
    ! turns a character string to upper case
    character(len=*), intent(in)  :: instr
    character(len=FILENAMELEN)    :: outstr

    integer(i4b) :: i, ascii, la, ua

    la = iachar('a')
    ua = iachar('A')

    outstr = instr
    do i = 1, min(len_trim(instr),len_trim(outstr))
       ascii = iachar( instr(i:i) )
       if (ascii >= la .and. ascii < la+26) then ! in [a,z]
          outstr(i:i) = achar( ascii - la + ua )
       endif
    enddo

    return
  end function strupcase

    ! ===========================================================
    subroutine exit_with_status (code, msg)
      ! ===========================================================
      integer, intent(in) :: code
      character (len=*), intent(in), optional :: msg
      ! ===========================================================

      if (present(msg)) print *,trim(msg)
      print *,'program exits with exit code ', code

      call exit (code)
    end subroutine exit_with_status



subroutine map_bad_pixels_d(map, fin, fout, nbads, verbose)
  use long_intrinsic_smw, only: long_size
  real(DP), dimension(0:,1:), intent(inout) :: map
  real(DP),                   intent(in)    :: fin, fout
  integer(I8B),                 intent(out)   :: nbads
  logical(LGT),     optional,   intent(in)    :: verbose
  integer(I8B) :: i, npix
  integer(I8B), dimension(1:100) :: imiss
  integer(I4B) :: imap, nmaps
  logical(LGT) :: be_verbose
  real(DP) :: threshold
  !-----------------------------------------

  npix  = long_size(map, 1)
  nmaps = long_size(map, 2)
  threshold = abs(fin * 1.e-5_DP)

  imiss(:) = 0
  do imap = 1, nmaps
     do i=0,npix-1
        if ( abs( map(i,imap) - fin ) < threshold ) then
           map(i,imap) = fout
           imiss(imap) = imiss(imap)+1
        endif
     enddo
  enddo
  nbads = sum(imiss)

  be_verbose = .false.
  if (present(verbose)) be_verbose=verbose
  if (be_verbose) then
     write(*,'(a,1pe11.4)') 'blank value : ' ,fin
     do imap = 1, nmaps
        if (imiss(imap) > 0) then
           write(*,'(i7,a,f7.3,a,1pe11.4)') &
                &           imiss(imap),' missing pixels (', &
                &           (100.0_DP*imiss(imap))/npix,' %),'// &
                &           ' have been set to : ',fout
        endif
     enddo
  endif
 
  return
end subroutine map_bad_pixels_d



subroutine map_bad_pixels_s(map, fin, fout, nbads, verbose)
  use long_intrinsic_smw, only: long_size
  real(SP), dimension(0:,1:), intent(inout) :: map
  real(SP),                   intent(in)    :: fin, fout
  integer(I8B),                 intent(out)   :: nbads
  logical(LGT),     optional,   intent(in)    :: verbose
  integer(I8B) :: i, npix
  integer(I8B), dimension(1:100) :: imiss
  integer(I4B) :: imap, nmaps
  logical(LGT) :: be_verbose
  real(SP) :: threshold
  !-----------------------------------------

  npix  = long_size(map, 1)
  nmaps = long_size(map, 2)
  threshold = abs(fin * 1.e-5_SP)

  imiss(:) = 0
  do imap = 1, nmaps
     do i=0,npix-1
        if ( abs( map(i,imap) - fin ) < threshold ) then
           map(i,imap) = fout
           imiss(imap) = imiss(imap)+1
        endif
     enddo
  enddo
  nbads = sum(imiss)

  be_verbose = .false.
  if (present(verbose)) be_verbose=verbose
  if (be_verbose) then
     write(*,'(a,1pe11.4)') 'blank value : ' ,fin
     do imap = 1, nmaps
        if (imiss(imap) > 0) then
           write(*,'(i7,a,f7.3,a,1pe11.4)') &
                &           imiss(imap),' missing pixels (', &
                &           (100.0_SP*imiss(imap))/npix,' %),'// &
                &           ' have been set to : ',fout
        endif
     enddo
  endif

  return
end subroutine map_bad_pixels_s


end module 
