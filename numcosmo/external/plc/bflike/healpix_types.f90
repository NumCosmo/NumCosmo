!-----------------------------------------------------------------------------
!
!  Copyright(C) 1997-2005 Krzysztof M. Gorski, Eric Hivon,
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

!Module healpix_types
Module healpix_types
!lpl: renamed to avoid possible conflicts

    ! This module sets the types used in the Fortran 90 modules
    ! of the HEALPIX distribution and follows the example of Numerical Recipes
    !
    ! Benjamin D. Wandelt October 1997
    ! Eric Hivon June 1998
    ! Eric Hivon Oct  2001, edited to be compatible with 'F' compiler
    ! Eric Hivon July 2002, addition of i8b, i2b, i1b
    !                       addition of max_i8b, max_i2b and max_i1b
    !            Jan 2005, explicit form of max_i1b because of ifc 8.1.021
    !            June 2005, redefine i8b as 16 digit integer because of Nec f90 compiler

    Integer, Parameter, Public :: i8b = Selected_int_kind(16)
    Integer, Parameter, Public :: i4b = Selected_int_kind(9)
    Integer, Parameter, Public :: i2b = Selected_int_kind(4)
    Integer, Parameter, Public :: i1b = Selected_int_kind(2)
    Integer, Parameter, Public :: sp = Selected_real_kind(5, 30)
    Integer, Parameter, Public :: dp = Selected_real_kind(12, 200)
    Integer, Parameter, Public :: lgt = Kind(.True.)
    Integer, Parameter, Public :: spc = Kind((1.0_sp, 1.0_sp))
    Integer, Parameter, Public :: dpc = Kind((1.0_dp, 1.0_dp))
    !
    Integer(i8b), Parameter, Public :: max_i8b = Huge(1_i8b)
    Integer, Parameter, Public :: max_i4b = Huge(1_i4b)
    Integer, Parameter, Public :: max_i2b = Huge(1_i2b)
    Integer, Parameter, Public :: max_i1b = 127
    Real(Kind=sp), Parameter, Public :: max_sp = Huge(1.0_sp)
    Real(Kind=dp), Parameter, Public :: max_dp = Huge(1.0_dp)

    ! Numerical Constant(Double precision)
    Real(Kind=dp), Parameter, Public :: QUARTPI = 0.785398163397448309615660845819875721049_dp
    Real(Kind=dp), Parameter, Public :: HALFPI = 1.570796326794896619231321691639751442099_dp
    Real(Kind=dp), Parameter, Public :: PI = 3.141592653589793238462643383279502884197_dp
    Real(Kind=dp), Parameter, Public :: TWOPI = 6.283185307179586476925286766559005768394_dp
    Real(Kind=dp), Parameter, Public :: FOURPI = 12.56637061435917295385057353311801153679_dp
    Real(Kind=dp), Parameter, Public :: SQRT2 = 1.41421356237309504880168872420969807856967_dp
    Real(Kind=dp), Parameter, Public :: SQRT2_INV = 0.707106781186547524400844362105_dp
    Real(Kind=dp), Parameter, Public :: EULER = 0.5772156649015328606065120900824024310422_dp
    Real(Kind=dp), Parameter, Public :: SQ4PI_INV = 0.2820947917738781434740397257803862929220_dp
    Real(Kind=dp), Parameter, Public :: TWOTHIRD = 0.6666666666666666666666666666666666666666_dp

    Real(Kind=dp), Parameter, Public :: RAD2DEG = 180.0_dp / PI
    Real(Kind=dp), Parameter, Public :: DEG2RAD = PI / 180.0_dp
    Real(Kind=sp), Parameter, Public :: hpx_sbadval = - 1.6375e30_sp
    Real(Kind=dp), Parameter, Public :: hpx_dbadval = - 1.6375e30_dp

    ! Maximum length of filenames
    Integer, Parameter :: filenamelen = 1024


    ! ---- Normalisation and convention ----
    ! normalisation of spin weighted functions
    Real(Kind=dp), Parameter, Public :: KvS = 1.0_dp ! 1.0 : CMBFAST(Healpix 1.2)
    ! sign of Q
    Real(Kind=dp), Parameter, Public :: sgQ = - 1.0_dp ! -1 : CMBFAST(Healpix 1.2)
    ! sign of spin weighted function !
    Real(Kind=dp), Parameter, Public :: SW1 = - 1.0_dp ! -1 : Healpix 1.2, bug correction

    !  ! normalisation of spin weighted functions
    !  real(kind=dp), parameter, public ::  KvS = 2.0_dp ! 2.0 : KKS  (Healpix 1.1)
    !  ! sign of Q
    !  real(kind=dp), parameter, public :: sgQ = +1.0_dp ! +1 : KKS(Healpix 1.1)
    !  ! sign of spin weighted function !
    !  real(kind=dp), parameter, public :: SW1 = +1.0_dp ! +1 : Healpix 1.1

    Real(Kind=dp), Parameter, Public :: iKvS = 1.0_dp / KvS ! inverse of KvS

  End Module healpix_types
