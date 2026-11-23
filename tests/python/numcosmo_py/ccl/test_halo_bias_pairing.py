#
# test_halo_bias_pairing.py
#
# Mon Nov 16 2026
# Copyright  2026  Disrael Camargo
# <disrael.camargo@gmail.com>
#
# test_halo_bias_pairing.py
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

"""Unit tests for the CCL/NumCosmo halo bias pairings.

Each constructor must return a CCL halo bias, the NumCosmo multiplicity function its
mass function is built on, and a function building the NumCosmo halo bias. Both sides
must carry the configuration's mass definition, and every parameter of the NumCosmo
bias must equal the one CCL's object uses, with delta_c also read from CCL's
get_delta_c. The NumCosmo side is built on a mass function that is constructed but
never prepared, so no cosmology is computed.
"""

from collections.abc import Callable
from typing import Any

import pytest

pytest.importorskip("pyccl")
import pyccl
from pyccl.halos.halo_model_base import get_delta_c

import numcosmo_py.ccl.halo_bias_pairing as hbp
from numcosmo_py import Nc, Ncm

pytestmark = pytest.mark.ccl
Ncm.cfg_init()

MDef = Nc.MultiplicityFuncMassDef

# One call of each constructor of the pairing module.
CONSTRUCTOR_CALLS: list[tuple[Callable[..., hbp.HaloBiasPair], tuple[Any, ...]]] = [
    (hbp.sheth99_pair, ()),
    (hbp.sheth01_pair, ()),
    (hbp.tinker10_pair, ("critical", MDef.CRITICAL, 200)),
]


def _build_nc_bias(pair: hbp.HaloBiasPair) -> tuple[Nc.HaloMassFunction, Nc.HaloBias]:
    """Build the NumCosmo side of a pair on an unprepared mass function."""
    dist = Nc.Distance.new(2.0)
    ps = Nc.PowspecMLTransfer.new(Nc.TransferFuncEH.new())
    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    hmf = Nc.HaloMassFunction.new(dist, psf, pair.nc_mulf)

    return hmf, pair.nc_bias_new(hmf)


def _ccl_delta_c() -> float:
    """Return the delta_c CCL's three bias models use."""
    return get_delta_c(None, None, kind="EdS")


def test_sheth99_pair() -> None:
    """Check the Sheth & Tormen (1999) pair and its parameters."""
    pair = hbp.sheth99_pair()
    ccl_bias, nc_mulf, _ = pair
    hmf, nc_bias = _build_nc_bias(pair)

    assert isinstance(ccl_bias, pyccl.halos.HaloBiasSheth99)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncST)
    assert isinstance(nc_bias, Nc.HaloBiasSTSpher)
    assert ccl_bias.mass_def.Delta == "fof"
    assert nc_mulf.get_mdef() == MDef.MEAN
    assert nc_bias.peek_mass_function() is hmf

    assert ccl_bias.use_delta_c_fit is False
    assert nc_bias.get_delta_c() == _ccl_delta_c()
    assert nc_bias.get_delta_c() == hbp.DELTA_C_EDS
    assert nc_bias.get_a() == ccl_bias.a
    assert nc_bias.get_p() == ccl_bias.p


def test_sheth01_pair() -> None:
    """Check the Sheth, Mo & Tormen (2001) pair and its parameters."""
    pair = hbp.sheth01_pair()
    ccl_bias, nc_mulf, _ = pair
    hmf, nc_bias = _build_nc_bias(pair)

    assert isinstance(ccl_bias, pyccl.halos.HaloBiasSheth01)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncST)
    assert isinstance(nc_bias, Nc.HaloBiasSTEllip)
    assert ccl_bias.mass_def.Delta == "fof"
    assert nc_mulf.get_mdef() == MDef.MEAN
    assert nc_bias.peek_mass_function() is hmf

    assert nc_bias.get_delta_c() == _ccl_delta_c()
    assert nc_bias.get_delta_c() == ccl_bias.dc
    assert nc_bias.get_a() == ccl_bias.a
    assert nc_bias.get_b() == ccl_bias.b
    assert nc_bias.get_c() == ccl_bias.c


@pytest.mark.parametrize("rho_type, nc_mdef, Delta", hbp.TINKER10_BIAS_MDEFS)
def test_tinker10_pair(
    rho_type: str, nc_mdef: Nc.MultiplicityFuncMassDef, Delta: int
) -> None:
    """Check the Tinker et al. (2010) pair, its mass definition and parameters."""
    pair = hbp.tinker10_pair(rho_type, nc_mdef, Delta)
    ccl_bias, nc_mulf, _ = pair
    hmf, nc_bias = _build_nc_bias(pair)

    assert isinstance(ccl_bias, pyccl.halos.HaloBiasTinker10)
    assert isinstance(nc_mulf, Nc.MultiplicityFuncTinker)
    assert isinstance(nc_bias, Nc.HaloBiasTinker)
    assert ccl_bias.mass_def.rho_type == rho_type
    assert ccl_bias.mass_def.Delta == Delta
    assert nc_mulf.get_mdef() == nc_mdef
    assert nc_mulf.get_Delta() == float(Delta)
    assert nc_mulf.get_property("linear-interp") is True
    assert nc_bias.peek_mass_function() is hmf

    assert nc_bias.get_delta_c() == _ccl_delta_c()
    assert nc_bias.get_delta_c() == ccl_bias.dc
    assert nc_bias.get_B() == ccl_bias.B
    assert nc_bias.get_b() == ccl_bias.b
    assert nc_bias.get_c() == ccl_bias.c


def test_unaligned_parameters_are_numcosmo_defaults() -> None:
    """Check that the parameters other than delta_c and a are NumCosmo's defaults.

    The pairs set NumCosmo's delta_c, and a for Sheth & Tormen (1999), to CCL's
    values; every other parameter equals that of the object NumCosmo's plain
    constructor builds.
    """
    hmf, st_spher = _build_nc_bias(hbp.sheth99_pair())
    default_spher = Nc.HaloBiasSTSpher.new(hmf)
    assert isinstance(st_spher, Nc.HaloBiasSTSpher)
    assert st_spher.get_p() == default_spher.get_p()

    hmf, st_ellip = _build_nc_bias(hbp.sheth01_pair())
    default_ellip = Nc.HaloBiasSTEllip.new(hmf)
    assert isinstance(st_ellip, Nc.HaloBiasSTEllip)
    assert st_ellip.get_a() == default_ellip.get_a()
    assert st_ellip.get_b() == default_ellip.get_b()
    assert st_ellip.get_c() == default_ellip.get_c()

    hmf, tinker = _build_nc_bias(hbp.tinker10_pair("critical", MDef.CRITICAL, 200))
    default_tinker = Nc.HaloBiasTinker.new(hmf)
    assert isinstance(tinker, Nc.HaloBiasTinker)
    assert tinker.get_B() == default_tinker.get_B()
    assert tinker.get_b() == default_tinker.get_b()
    assert tinker.get_c() == default_tinker.get_c()


def test_tinker10_bias_mdefs_are_critical() -> None:
    """Check that every Tinker et al. (2010) definition paired is critical."""
    assert hbp.TINKER10_BIAS_MDEFS
    for rho_type, nc_mdef, _ in hbp.TINKER10_BIAS_MDEFS:
        assert rho_type == "critical"
        assert nc_mdef == MDef.CRITICAL


def test_tinker10_pair_rejects_matter() -> None:
    """Check that tinker10_pair raises ValueError for a mean mass definition."""
    with pytest.raises(ValueError, match="only critical mass definitions"):
        hbp.tinker10_pair("matter", MDef.MEAN, 200)


def test_pairs_cover_every_constructor() -> None:
    """Check that CONSTRUCTOR_CALLS calls every *_pair constructor of the module."""
    module_constructors = {
        name
        for name, value in vars(hbp).items()
        if name.endswith("_pair") and callable(value)
    }
    called = {constructor.__name__ for constructor, _ in CONSTRUCTOR_CALLS}
    missing = sorted(module_constructors - called)

    assert not missing, f"CONSTRUCTOR_CALLS does not call: {', '.join(missing)}"


def test_halo_bias_pair() -> None:
    """Check HaloBiasPair's fields and that every constructor returns one.

    compare_halo_bias and compare_halo_bias_sigma take the three fields positionally,
    in the order of HaloBiasPair._fields.
    """
    ccl_bias, nc_mulf, nc_bias_new = hbp.sheth99_pair()
    pair = hbp.HaloBiasPair(ccl_bias, nc_mulf, nc_bias_new)

    assert hbp.HaloBiasPair._fields == ("ccl_bias", "nc_mulf", "nc_bias_new")
    assert pair.ccl_bias is ccl_bias
    assert pair.nc_mulf is nc_mulf
    assert pair.nc_bias_new is nc_bias_new

    for constructor, args in CONSTRUCTOR_CALLS:
        assert isinstance(constructor(*args), hbp.HaloBiasPair)
