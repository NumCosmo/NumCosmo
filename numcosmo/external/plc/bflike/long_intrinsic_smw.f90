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

!module long_intrinsic
module long_intrinsic_smw
!lpl: renamed to avoid possible conflicts

  ! If correctly compiled this module redefines 
  ! - the intrisic function SIZE  (renamed long_size)
  ! - the intrisic function COUNT (renamed long_count)
  ! to return 8-byte = 64-bit integer (I8B) results instead of default type integer.
  !
  ! It must be compiled with the option
  ! Intel ifort:    -i8   or   -integer-size 64
  ! GNU gfortran: -fdefault-integer-8
  ! IBM xlf:      -qintsize=8
  ! NAG f95:      -double
  ! Pathscale pathf90:     -i8  or -default64
  ! g95:                   -i8
  !
  ! 2009-02-25: version 1.0, EH/IAP
  !---------------------------------------------------------------------

  use healpix_types, only: i4b, i8b, sp, dp, lgt
  implicit none

  private
  public:: long_size, long_count

  interface long_size
     module procedure size_i_1, size_i_2, &
          &           size_r_1, size_r_2, &
          &           size_d_1, size_d_2
  end interface

  interface long_count
     module procedure count_l_1
  end interface


contains
  !============================================
  !  integer I4B arrays
  !============================================
  function size_i_1(array, dim) result (mysize)
    integer(i4b), dimension(:), intent(in)   :: array
    integer(i4b), optional,     intent(in)   :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_i_1

  function size_i_2(array, dim) result (mysize)
    integer(i4b), dimension(:,:), intent(in) :: array
    integer(i4b), optional ,      intent(in) :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_i_2

  !============================================
  !  integer I8B arrays
  !============================================
  function size_j_1(array, dim) result (mysize)
    integer(i8b), dimension(:), intent(in)   :: array
    integer(i4b), optional,     intent(in)   :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_j_1

  function size_j_2(array, dim) result (mysize)
    integer(i8b), dimension(:,:), intent(in) :: array
    integer(i4b), optional ,      intent(in) :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_j_2
  !============================================
  !  real SP arrays
  !============================================
  function size_r_1(array, dim) result (mysize)
    real(sp),     dimension(:), intent(in)   :: array
    integer(i4b), optional,     intent(in)   :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_r_1

  function size_r_2(array, dim) result (mysize)
    real(sp),     dimension(:,:), intent(in) :: array
    integer(i4b), optional ,      intent(in) :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_r_2

  !============================================
  !  real DP arrays
  !============================================
  function size_d_1(array, dim) result (mysize)
    real(dp),     dimension(:), intent(in)   :: array
    integer(i4b), optional,     intent(in)   :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_d_1

  function size_d_2(array, dim) result (mysize)
    real(dp),     dimension(:,:), intent(in) :: array
    integer(i4b), optional ,      intent(in) :: dim
    integer(i8b)                             :: mysize
    if (present(dim)) then
       mysize = size(array, dim=dim)
    else
       mysize = size(array) ! required for xlf
    endif
    return
  end function size_d_2

  !*************** COUNT ***********************
  !============================================
  !  logical array
  !============================================
  function count_l_1(array) result (mycount)
    logical(lgt), dimension(:), intent(in)   :: array
    integer(i8b)                             :: mycount
    mycount = count(array)
    return
  end function count_l_1


end module
