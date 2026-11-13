#
# test_multiplicity_func_pairing.py
#
# Tue Oct 06 11:19:25 2026
# Copyright  2026  Disrael Camargo
# <disrael.camargo@gmail.com>
#
# test_multiplicity_func_pairing.py
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

"""Unit tests for the CCL/NumCosmo mass-function pairings.

Each constructor must return a CCL mass function and a NumCosmo multiplicity function
with the same mass definition (overdensity and reference density) and the same family
options. Press & Schechter (1974) and Sheth & Tormen (1999) are the exception: their
f(sigma) does not depend on the mass definition, and each code accepts only one for
them. Each configuration's NumCosmo mass definition is also checked against the one
CCL's (rho_type, Delta) names, so a configuration changed on one side only fails.
"""

from collections.abc import Callable

import pytest

pytest.importorskip("pyccl")
# flake8: noqa: E402
import pyccl

import numcosmo_py.ccl.multiplicity_func_pairing as mfp
from numcosmo_py import Nc, Ncm

pytestmark = pytest.mark.ccl
Ncm.cfg_init()

MDef = Nc.MultiplicityFuncMassDef


def _expected_mdef(rho_type: str, Delta: int | str) -> Nc.MultiplicityFuncMassDef:
    """Return the NumCosmo mass definition naming CCL's (rho_type, Delta)."""
    if Delta == "fof":
        return MDef.FOF
    if Delta == "vir":
        return MDef.VIRIAL
    return {"matter": MDef.MEAN, "critical": MDef.CRITICAL}[rho_type]


def _check_mass_def(
    pair: mfp.MassFunctionPair,
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    Delta: int | str,
) -> None:
    """Check that both sides of a pair carry the configuration's mass definition.

    The configuration itself is checked first: its NumCosmo mass definition must be
    the one CCL's (rho_type, Delta) names. A string Delta ("fof", "vir") passes no
    overdensity, so NumCosmo's Delta must be left at its default.
    """
    ccl_mf, nc_mulf = pair

    assert nc_mdef == _expected_mdef(rho_type, Delta)
    assert ccl_mf.mass_def.rho_type == rho_type
    assert ccl_mf.mass_def.Delta == Delta
    assert nc_mulf.get_mdef() == nc_mdef

    if isinstance(Delta, str):
        assert nc_mulf.get_Delta() == type(nc_mulf).new().get_Delta()
    else:
        assert nc_mulf.get_Delta() == float(Delta)


def test_tinker08_critical_deltas_stop_at_800() -> None:
    """Check that critical Deltas stop at 800.

    Delta_c = 800 enters the fit as Delta_m = 800 / Omega_m(a), which stays within
    CCL's table (Delta_m <= 3200) for Omega_m(a) >= 0.25.
    """
    assert max(mfp.TINKER08_DELTAS_CRITICAL) <= 800


@pytest.mark.parametrize(
    "rho_type, nc_mdef, Delta",
    [("matter", MDef.MEAN, Delta) for Delta in mfp.TINKER08_DELTAS_MATTER]
    + [("critical", MDef.CRITICAL, Delta) for Delta in mfp.TINKER08_DELTAS_CRITICAL],
)
def test_tinker08_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int
) -> None:
    """Check the Tinker (2008) pair, linearly interpolated by default."""
    pair = mfp.tinker08_pair(rho_type, nc_mdef, Delta)
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncTinker08)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncTinker)
    _check_mass_def(pair, rho_type, nc_mdef, Delta)
    assert nc_mulf.get_property("linear-interp") is True


@pytest.mark.parametrize("linear_interp", [True, False])
def test_tinker08_pair_linear_interp(linear_interp: bool) -> None:
    """Check that the interpolation choice reaches the NumCosmo side."""
    _, nc_mulf = mfp.tinker08_pair(
        "matter", MDef.MEAN, 500, linear_interp=linear_interp
    )

    assert nc_mulf.get_property("linear-interp") is linear_interp


@pytest.mark.parametrize("Delta", mfp.TINKER10_DELTAS)
def test_tinker10_pair(Delta: int) -> None:
    """Check the Tinker (2010) pair, mean overdensity."""
    pair = mfp.tinker10_pair("matter", MDef.MEAN, Delta)
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncTinker10)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncTinkerMeanNormalized)
    _check_mass_def(pair, "matter", MDef.MEAN, Delta)


@pytest.mark.parametrize("ccl_hydro, nc_sim", mfp.BOCQUET_HYDRO)
@pytest.mark.parametrize("rho_type, nc_mdef, Delta", mfp.BOCQUET_MDEFS)
def test_bocquet_pair(
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    Delta: int,
    ccl_hydro: bool,
    nc_sim: Nc.MultiplicityFuncBocquetSim,
) -> None:
    """Check the Bocquet (2016) pair and its simulation set on both sides."""
    pair = mfp.bocquet_pair(rho_type, nc_mdef, ccl_hydro, nc_sim, Delta)
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncBocquet16)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncBocquet)
    _check_mass_def(pair, rho_type, nc_mdef, Delta)
    assert ccl_mf.hydro is ccl_hydro
    assert nc_mulf.get_sim() == nc_sim
    assert (nc_sim == Nc.MultiplicityFuncBocquetSim.HYDRO) is ccl_hydro


@pytest.mark.parametrize("rho_type, nc_mdef, Delta", mfp.WATSON_OPTS)
def test_watson_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int | str
) -> None:
    """Check the Watson (2013) pair, mean overdensity or friends-of-friends."""
    pair = mfp.watson_pair(rho_type, nc_mdef, Delta)
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncWatson13)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncWatson)
    _check_mass_def(pair, rho_type, nc_mdef, Delta)


@pytest.mark.parametrize("rho_type, nc_mdef, ellipsoidal, Delta", mfp.DESPALI_OPTS)
def test_despali_pair(
    rho_type: str,
    nc_mdef: Nc.MultiplicityFuncMassDef,
    ellipsoidal: bool,
    Delta: int | str,
) -> None:
    """Check the Despali (2016) pair and its ellipsoidal correction on both sides."""
    pair = mfp.despali_pair(rho_type, nc_mdef, ellipsoidal, Delta)
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncDespali16)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncDespali)
    _check_mass_def(pair, rho_type, nc_mdef, Delta)
    assert ccl_mf.ellipsoidal is ellipsoidal
    assert nc_mulf.get_eo() is ellipsoidal


def test_despali_opt_covers_deltas() -> None:
    """Check that DESPALI_OPTS is every Delta, both densities, both corrections."""
    expected = {
        (rho_type, Delta, ellipsoidal)
        for rho_type in ("matter", "critical")
        for Delta in mfp.DESPALI_DELTAS
        for ellipsoidal in (True, False)
    } | {("critical", "vir", True), ("critical", "vir", False)}

    assert {(r, D, e) for r, _, e, D in mfp.DESPALI_OPTS} == expected
    assert len(mfp.DESPALI_OPTS) == len(expected)


def test_jenkins01_pair() -> None:
    """Check the Jenkins (2001) pair: friends-of-friends on both sides."""
    pair = mfp.jenkins01_pair()
    ccl_mf, nc_mulf = pair

    assert isinstance(ccl_mf, pyccl.halos.MassFuncJenkins01)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncJenkins)
    _check_mass_def(pair, "matter", MDef.FOF, "fof")


@pytest.mark.parametrize(
    "constructor, ccl_type, nc_type",
    [
        (mfp.press74_pair, pyccl.halos.MassFuncPress74, Nc.MultiplicityFuncPS),
        (mfp.sheth99_pair, pyccl.halos.MassFuncSheth99, Nc.MultiplicityFuncST),
    ],
)
def test_universal_pair(
    constructor: Callable[[], mfp.MassFunctionPair],
    ccl_type: type[pyccl.halos.MassFunc],
    nc_type: type[Nc.MultiplicityFunc],
) -> None:
    """Check the Press & Schechter and Sheth & Tormen pairs.

    CCL accepts only the friends-of-friends mass definition for both, and NumCosmo
    only MEAN, without an overdensity. Their f(sigma) depends only on delta_c / sigma,
    so the mass definition changes no value on either side.
    """
    ccl_mf, nc_mulf = constructor()

    assert isinstance(ccl_mf, ccl_type)
    assert isinstance(nc_mulf, nc_type)
    assert ccl_mf.mass_def.Delta == "fof"
    assert nc_mulf.get_mdef() == MDef.MEAN
