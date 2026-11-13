#
# multiplicity_func_pairing.py
#
# Tue Oct 06 11:19:25 2026
# Copyright  2026  Disrael Camargo
# <disrael.camargo@gmail.com>
#
# multiplicity_func_pairing.py
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

"""Mass-function pairings between CCL and NumCosmo.

This module records, for each multiplicity-function family, which CCL mass function
and mass definition correspond to which NumCosmo multiplicity function and
:obj:`~Nc.MultiplicityFuncMassDef`. Neither library provides this correspondence.

The module defines two kinds of objects:

* lists of overdensities and option tuples enumerating the covered configurations,
  used to parametrize tests and to build benchmark tables; and
* functions returning a ``(CCL mass function, NumCosmo multiplicity function)`` pair
  for one configuration.

The functions return the multiplicity function, not an :obj:`~Nc.HaloMassFunction`,
because building one needs the caller's distance and power-spectrum filter.
"""

import pyccl

from numcosmo_py import Nc

# A CCL mass function paired with the NumCosmo multiplicity function that
# represents the same mass definition.
MassFunctionPair = tuple[pyccl.halos.MassFunc, Nc.MultiplicityFunc]


# Tinker (2008). The fit coefficients are tabulated at the mean overdensities
# Delta_m = 200, 300, 400, 600, 800, 1200, 1600, 2400 and 3200; other values are
# interpolated in log10(Delta_m). A critical Delta_c enters the fit as
# Delta_m = Delta_c / Omega_m(a), and CCL raises ValueError for Delta_m above 3200.
# The critical entries stop at 800, which stays within the table for Omega_m(a) >= 0.25;
# Delta_c = 1200 leaves it for Omega_m(a) < 0.375.
TINKER08_DELTAS_MATTER: list[int] = [
    200,
    300,
    400,
    500,
    600,
    800,
    1200,
    1600,
    2400,
    3200,
]
TINKER08_DELTAS_CRITICAL: list[int] = [200, 300, 400, 500, 600, 800]

# Tinker (2010), mean overdensity only.
TINKER10_DELTAS: list[int] = [200, 300, 400, 600, 800, 1200, 1600, 2400, 3200]

# Bocquet (2016): the fit is calibrated at 200m, 200c and 500c.
BOCQUET_MDEFS: list[tuple[str, Nc.MultiplicityFuncMassDef, int]] = [
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 200),
    ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, 200),
    ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, 500),
]

# CCL selects the hydrodynamical or dark-matter-only calibration with a boolean,
# NumCosmo with an enumeration.
BOCQUET_HYDRO: list[tuple[bool, Nc.MultiplicityFuncBocquetSim]] = [
    (True, Nc.MultiplicityFuncBocquetSim.HYDRO),
    (False, Nc.MultiplicityFuncBocquetSim.DM),
]

# Watson (2013). CCL writes the friends-of-friends mass definition as the string
# "fof", and NumCosmo as the enumeration value FOF; this definition has no Delta.
WATSON_OPTS: list[tuple[str, Nc.MultiplicityFuncMassDef, int | str]] = [
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 178),
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 200),
    ("matter", Nc.MultiplicityFuncMassDef.MEAN, 500),
    ("matter", Nc.MultiplicityFuncMassDef.FOF, "fof"),
]

# Despali (2016), with and without the ellipsoidal-collapse correction. CCL writes the
# virial mass definition as the string "vir", and NumCosmo as the enumeration value
# VIRIAL; the virial fit depends on delta_c(z), which each code computes from the
# cosmology, and takes no Delta.
DESPALI_DELTAS: list[int] = [200, 300, 400, 500, 600]
DESPALI_OPTS: list[tuple[str, Nc.MultiplicityFuncMassDef, bool, int | str]] = (
    [
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, True, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("matter", Nc.MultiplicityFuncMassDef.MEAN, False, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, True, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.CRITICAL, False, Delta)
        for Delta in DESPALI_DELTAS
    ]
    + [
        ("critical", Nc.MultiplicityFuncMassDef.VIRIAL, True, "vir"),
        ("critical", Nc.MultiplicityFuncMassDef.VIRIAL, False, "vir"),
    ]
)


def press74_pair() -> MassFunctionPair:
    """Build the Press & Schechter (1974) pair.

    CCL accepts only the friends-of-friends mass definition for this function, and
    NumCosmo only MEAN, without an overdensity. The function depends only on
    delta_c / sigma, so both sides compute the same values.
    """
    return (
        pyccl.halos.MassFuncPress74(mass_def=pyccl.halos.MassDefFof),
        Nc.MultiplicityFuncPS.new(),
    )


def sheth99_pair() -> MassFunctionPair:
    """Build the Sheth & Tormen (1999) pair.

    CCL accepts only the friends-of-friends mass definition for this function, and
    NumCosmo only MEAN, without an overdensity. The function depends only on
    delta_c / sigma, so both sides compute the same values.
    """
    return (
        pyccl.halos.MassFuncSheth99(mass_def=pyccl.halos.MassDefFof),
        Nc.MultiplicityFuncST.new(),
    )


def jenkins01_pair() -> MassFunctionPair:
    """Build the Jenkins (2001) pair, friends-of-friends.

    Both libraries default to the friends-of-friends mass definition, so neither side
    is given one explicitly.
    """
    return pyccl.halos.MassFuncJenkins01(), Nc.MultiplicityFuncJenkins.new()


def tinker08_pair(
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    Delta: int,
    *,
    linear_interp: bool = True,
) -> MassFunctionPair:
    """Build a Tinker (2008) pair for one mass definition.

    ``linear_interp`` selects how NumCosmo interpolates the fit coefficients in
    log10(Delta_m): ``True`` is linear, as in CCL, and ``False``, NumCosmo's default,
    is the cubic spline recommended in Tinker (2008), Appendix B. The setting changes
    the result only when Delta_m is not tabulated. For a mean mass definition
    Delta_m = Delta; for a critical one Delta_m = Delta / Omega_m(a).
    """
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(nc_mdef)
    mulf.set_Delta(Delta)
    mulf.set_linear_interp(linear_interp)

    return (
        pyccl.halos.MassFuncTinker08(mass_def=pyccl.halos.MassDef(Delta, rho_type)),
        mulf,
    )


def tinker10_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int
) -> MassFunctionPair:
    """Build a Tinker (2010) pair for one mean mass definition."""
    mulf = Nc.MultiplicityFuncTinkerMeanNormalized.new()
    mulf.set_mdef(nc_mdef)
    mulf.set_Delta(Delta)

    return (
        pyccl.halos.MassFuncTinker10(mass_def=pyccl.halos.MassDef(Delta, rho_type)),
        mulf,
    )


def bocquet_pair(
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    ccl_hydro: bool,
    nc_sim: Nc.MultiplicityFuncBocquetSim,
    Delta: int,
) -> MassFunctionPair:
    """Build a Bocquet (2016) pair for one mass definition and simulation set."""
    return (
        pyccl.halos.MassFuncBocquet16(
            hydro=ccl_hydro, mass_def=pyccl.halos.MassDef(Delta, rho_type)
        ),
        Nc.MultiplicityFuncBocquet.new_full(nc_mdef, nc_sim, Delta),
    )


def watson_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int | str
) -> MassFunctionPair:
    """Build a Watson (2013) pair.

    ``Delta`` is a mean overdensity, or the string ``"fof"`` for the
    friends-of-friends calibration, which carries no overdensity.
    """
    mulf = Nc.MultiplicityFuncWatson.new()
    mulf.set_mdef(nc_mdef)
    if nc_mdef != Nc.MultiplicityFuncMassDef.FOF:
        assert isinstance(Delta, int)
        mulf.set_Delta(Delta)

    return (
        pyccl.halos.MassFuncWatson13(mass_def=pyccl.halos.MassDef(Delta, rho_type)),
        mulf,
    )


def despali_pair(
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    ellipsoidal: bool,
    Delta: int | str,
) -> MassFunctionPair:
    """Build a Despali (2016) pair.

    ``Delta`` is an overdensity, or the string ``"vir"`` for the virial mass
    definition, whose fit depends on delta_c(z) and takes no overdensity.
    ``ellipsoidal`` selects the ellipsoidal-collapse correction on both sides.
    """
    mulf = Nc.MultiplicityFuncDespali.new()
    mulf.set_mdef(nc_mdef)
    mulf.set_eo(ellipsoidal)
    if nc_mdef != Nc.MultiplicityFuncMassDef.VIRIAL:
        assert isinstance(Delta, int)
        mulf.set_Delta(Delta)

    return (
        pyccl.halos.MassFuncDespali16(
            ellipsoidal=ellipsoidal, mass_def=pyccl.halos.MassDef(Delta, rho_type)
        ),
        mulf,
    )
