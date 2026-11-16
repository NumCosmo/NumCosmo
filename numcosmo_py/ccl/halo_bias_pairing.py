#
# halo_bias_pairing.py
#
# Mon Nov 16 2026
# Copyright  2026  Disrael Camargo
# <disrael.camargo@gmail.com>
#
# halo_bias_pairing.py
# Copyright (C) 2026 Disrael Camargo <disrael.camargo@gmail.com>
#
# numcosmo is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# numcosmo is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program.  If not, see <http://www.gnu.org/licenses/>.

r"""Halo bias pairings between CCL and NumCosmo.

This module records, for each halo bias model, which CCL bias object corresponds to
which NumCosmo :obj:`~Nc.HaloBias`, and with which parameters. Neither library
provides this correspondence.

A NumCosmo bias object is built on an :obj:`~Nc.HaloMassFunction`, which needs the
caller's distance and power-spectrum filter. Each constructor therefore returns a
:class:`HaloBiasPair`: the CCL bias object, the NumCosmo multiplicity function to
build that mass function on (from
:mod:`numcosmo_py.ccl.multiplicity_func_pairing`), and a function that builds the
NumCosmo bias object from the mass function.

The pairs are:

* Sheth & Tormen (1999): :class:`pyccl.halos.HaloBiasSheth99` and
  :obj:`~Nc.HaloBiasSTSpher`;
* Sheth, Mo & Tormen (2001): :class:`pyccl.halos.HaloBiasSheth01` and
  :obj:`~Nc.HaloBiasSTEllip`;
* Tinker et al. (2010): :class:`pyccl.halos.HaloBiasTinker10` and
  :obj:`~Nc.HaloBiasTinker`, at critical mass definitions, on the Tinker (2008)
  multiplicity function. Its linear interpolation matches
  :func:`numcosmo_py.ccl.multiplicity_func_pairing.tinker08_pair` and has no effect on
  :math:`b(\sigma)`.

NumCosmo's :math:`\delta_c` is set to CCL's value, 1.68647019984, in all three
pairs, and NumCosmo's :math:`a` to CCL's 0.707 in the Sheth & Tormen (1999) pair.
Every other parameter is passed at the value both codes use by default.

By default (``mass_def_strict=True``), CCL accepts only the friends-of-friends mass
definition for the two Sheth & Tormen models. On NumCosmo's side they are built on
the Sheth & Tormen multiplicity function, whose only mass definition is MEAN.
Neither code's :math:`b(\sigma)` for these models reads the mass definition.

Mean-density mass definitions of Tinker et al. (2010) are not included.
"""

from collections.abc import Callable
from typing import NamedTuple

import pyccl

import numcosmo_py.ccl.multiplicity_func_pairing as mfp
from numcosmo_py import Nc


class HaloBiasPair(NamedTuple):
    """A CCL halo bias and the parts that build the matching NumCosmo halo bias.

    ``nc_bias_new`` takes an :obj:`~Nc.HaloMassFunction` built on ``nc_mulf`` and
    returns the NumCosmo bias object.
    """

    ccl_bias: pyccl.halos.HaloBias
    nc_mulf: Nc.MultiplicityFunc
    nc_bias_new: Callable[[Nc.HaloMassFunction], Nc.HaloBias]


# Linear collapse threshold that CCL's three bias models use: get_delta_c with
# kind='EdS', delta_c = (3/20)(12 pi)^(2/3) (pyccl/halos/halo_model_base.py:28-30;
# used at hbias/sheth99.py:44, hbias/sheth01.py:32, hbias/tinker10.py:36).
DELTA_C_EDS = 1.68647019984

# Sheth & Tormen (1999), eq. (12), with a = 0.707 and p = 0.3 as given below their
# eq. (10) (pyccl/halos/hbias/sheth99.py:37-38).
ST99_A = 0.707
ST99_P = 0.3

# Sheth, Mo & Tormen (2001), eq. (8), with a = 0.707, b = 0.5 and c = 0.6
# (pyccl/halos/hbias/sheth01.py:28, 30, 31).
SMT01_A = 0.707
SMT01_B = 0.5
SMT01_C = 0.6

# Tinker et al. (2010), eq. (6) and Table 2: B = 0.183, b = 1.5, c = 2.4
# (pyccl/halos/hbias/tinker10.py:33-35).
TINKER10_B_UPPER = 0.183
TINKER10_B = 1.5
TINKER10_C = 2.4

# Tinker et al. (2010): the critical mass definitions paired. The model takes
# Delta_m = Delta_c / Omega_m(z) in both codes.
TINKER10_BIAS_MDEFS: list[tuple[str, Nc.MultiplicityFuncMassDef, int]] = [
    ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, 200),
    ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, 500),
]


def _st_spher_new(hmf: Nc.HaloMassFunction) -> Nc.HaloBias:
    """Build NumCosmo's Sheth & Tormen (1999) bias with CCL's parameters."""
    return Nc.HaloBiasSTSpher.new_full(hmf, DELTA_C_EDS, ST99_A, ST99_P)


def _st_ellip_new(hmf: Nc.HaloMassFunction) -> Nc.HaloBias:
    """Build NumCosmo's Sheth, Mo & Tormen (2001) bias with CCL's parameters."""
    return Nc.HaloBiasSTEllip.new_full(hmf, DELTA_C_EDS, SMT01_A, SMT01_B, SMT01_C)


def _tinker_new(hmf: Nc.HaloMassFunction) -> Nc.HaloBias:
    """Build NumCosmo's Tinker et al. (2010) bias with CCL's parameters."""
    return Nc.HaloBiasTinker.new_full(
        hmf, DELTA_C_EDS, TINKER10_B_UPPER, TINKER10_B, TINKER10_C
    )


def sheth99_pair() -> HaloBiasPair:
    """Build the Sheth & Tormen (1999) pair.

    CCL's side uses the friends-of-friends mass definition, and NumCosmo's is built on
    the Sheth & Tormen multiplicity function.
    """
    _, nc_mulf = mfp.sheth99_pair()

    return HaloBiasPair(
        pyccl.halos.HaloBiasSheth99(mass_def=pyccl.halos.MassDefFof),
        nc_mulf,
        _st_spher_new,
    )


def sheth01_pair() -> HaloBiasPair:
    """Build the Sheth, Mo & Tormen (2001) pair.

    CCL's side uses the friends-of-friends mass definition, and NumCosmo's is built on
    the Sheth & Tormen multiplicity function.
    """
    _, nc_mulf = mfp.sheth99_pair()

    return HaloBiasPair(
        pyccl.halos.HaloBiasSheth01(mass_def=pyccl.halos.MassDefFof),
        nc_mulf,
        _st_ellip_new,
    )


def tinker10_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int
) -> HaloBiasPair:
    """Build a Tinker et al. (2010) pair for one critical mass definition.

    NumCosmo's side is built on the Tinker (2008) multiplicity function at the same
    mass definition, with linear interpolation, to match
    :func:`numcosmo_py.ccl.multiplicity_func_pairing.tinker08_pair`; the interpolation
    has no effect on b(sigma).
    Raises :class:`ValueError` if ``rho_type`` is not ``"critical"``: only critical
    mass definitions are paired.
    """
    if rho_type != "critical":
        raise ValueError(
            f"rho_type = {rho_type!r}: only critical mass definitions are paired "
            "for Tinker et al. (2010)."
        )
    _, nc_mulf = mfp.tinker08_pair(rho_type, nc_mdef, Delta)

    return HaloBiasPair(
        pyccl.halos.HaloBiasTinker10(mass_def=pyccl.halos.MassDef(Delta, rho_type)),
        nc_mulf,
        _tinker_new,
    )
