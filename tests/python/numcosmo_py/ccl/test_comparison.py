#
# test_comparison.py
#
# Mon Oct 19 2026
# Copyright  2026  Disrael Camargo
# <disrael.camargo@gmail.com>
#
# test_comparison.py
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

"""Unit tests for the halo mass function comparison helpers.

The functions compare_mass_function and compare_multiplicity_function are checked
against CCL on the cosmologies of ccl_cosmo_eh_linear. compare_multiplicity_function
evaluates CCL's f(sigma) through the private method MassFunc._get_fsigma, which every
CCL mass function overrides; its signature is pinned on the base class and on each
family the pairing module builds, so that a pyccl release changing it fails here.
"""

import inspect
from collections.abc import Callable
from typing import Any

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal

pytest.importorskip("pyccl")
# flake8: noqa: E402
import pyccl

import numcosmo_py.ccl.comparison as nc_cmp
import numcosmo_py.ccl.multiplicity_func_pairing as mfp
from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import Cosmology

pytestmark = pytest.mark.ccl
pytest_plugins = ["python.fixtures_ccl"]
Ncm.cfg_init()

MDef = Nc.MultiplicityFuncMassDef
PairConstructor = Callable[..., mfp.MassFunctionPair]

# The parameters compare_multiplicity_function passes, in this order.
GET_FSIGMA_PARAMETERS = ["self", "cosmo", "sigM", "a", "lnM"]

# One call of each constructor of the pairing module; the CCL classes are read from
# the pairs built, so they are the ones the comparison actually calls.
CONSTRUCTOR_CALLS: list[tuple[PairConstructor, tuple[Any, ...]]] = [
    (mfp.press74_pair, ()),
    (mfp.sheth99_pair, ()),
    (mfp.jenkins01_pair, ()),
    (mfp.tinker08_pair, ("matter", MDef.MEAN, 200)),
    (mfp.tinker10_pair, ("matter", MDef.MEAN, 200)),
    (
        mfp.bocquet_pair,
        ("matter", MDef.MEAN, True, Nc.MultiplicityFuncBocquetSim.HYDRO, 200),
    ),
    (mfp.watson_pair, ("matter", MDef.MEAN, 200)),
    (mfp.despali_pair, ("matter", MDef.MEAN, True, 200)),
]
PAIRS = [constructor(*args) for constructor, args in CONSTRUCTOR_CALLS]
CCL_FAMILIES = [type(ccl_mf) for ccl_mf, _ in PAIRS]

# Bound on |f_CCL - f_NumCosmo| / f_NumCosmo for each family, at the same sigma, lnM
# and z: the smallest power of ten at least twice the largest difference measured over
# the three cosmologies of ccl_cosmo_eh_linear at both CCL precision settings,
# z = 0, 0.5, 1, 1.5 and every configuration of the family; 1.0e-13 where both codes
# evaluate the same expression, as in test_massfunc.py.
FSIGMA_RELTOL = {
    "press74": 1.0e-10,
    "sheth99": 1.0e-10,
    "jenkins01": 1.0e-13,
    "tinker08": 1.0e-9,
    "tinker10": 1.0e-13,
    "bocquet": 1.0e-9,
    "watson": 1.0e-8,
    "despali": 1.0e-9,
}
FSIGMA_Z = (0.0, 0.5, 1.0, 1.5)

# Every configuration the pairing module enumerates, as (family, constructor, args),
# with test_massfunc.py's ids.
MULTIPLICITY_CONFIGS = (
    [
        pytest.param("press74", mfp.press74_pair, (), id="press74"),
        pytest.param("sheth99", mfp.sheth99_pair, (), id="sheth99"),
        pytest.param("jenkins01", mfp.jenkins01_pair, (), id="jenkins01"),
    ]
    + [
        pytest.param(
            "tinker08",
            mfp.tinker08_pair,
            ("matter", MDef.MEAN, Delta),
            id=f"tinker08-matter_{Delta}",
        )
        for Delta in mfp.TINKER08_DELTAS_MATTER
    ]
    + [
        pytest.param(
            "tinker08",
            mfp.tinker08_pair,
            ("critical", MDef.CRITICAL, Delta),
            id=f"tinker08-critical_{Delta}",
        )
        for Delta in mfp.TINKER08_DELTAS_CRITICAL
    ]
    + [
        pytest.param(
            "tinker10",
            mfp.tinker10_pair,
            ("matter", MDef.MEAN, Delta),
            id=f"tinker10-mean_{Delta}",
        )
        for Delta in mfp.TINKER10_DELTAS
    ]
    + [
        pytest.param(
            "bocquet",
            mfp.bocquet_pair,
            (rho_type, nc_mdef, ccl_hydro, nc_sim, Delta),
            id=f"bocquet-{rho_type}_hydro_{ccl_hydro}_{Delta}",
        )
        for rho_type, nc_mdef, Delta in mfp.BOCQUET_MDEFS
        for ccl_hydro, nc_sim in mfp.BOCQUET_HYDRO
    ]
    + [
        pytest.param(
            "watson",
            mfp.watson_pair,
            (rho_type, nc_mdef, Delta),
            id=f"watson-{rho_type}_{Delta}",
        )
        for rho_type, nc_mdef, Delta in mfp.WATSON_OPTS
    ]
    + [
        pytest.param(
            "despali",
            mfp.despali_pair,
            (rho_type, nc_mdef, ellipsoidal, Delta),
            id=f"despali-{rho_type}_Ellip_{ellipsoidal}_{Delta}",
        )
        for rho_type, nc_mdef, ellipsoidal, Delta in mfp.DESPALI_OPTS
    ]
)


def _tinker08_200m() -> mfp.MassFunctionPair:
    """Build the Tinker (2008) 200m CCL/NumCosmo pair the tests compare."""
    return mfp.tinker08_pair("matter", MDef.MEAN, 200)


def test_ccl_get_fsigma_signature() -> None:
    """Check that MassFunc._get_fsigma takes the parameters the comparison passes."""
    assert hasattr(pyccl.halos.MassFunc, "_get_fsigma")
    assert callable(pyccl.halos.MassFunc._get_fsigma)
    parameters = inspect.signature(pyccl.halos.MassFunc._get_fsigma).parameters

    assert list(parameters) == GET_FSIGMA_PARAMETERS


@pytest.mark.parametrize("ccl_family", CCL_FAMILIES, ids=lambda cls: cls.__name__)
def test_ccl_family_get_fsigma_signature(
    ccl_family: type[pyccl.halos.MassFunc],
) -> None:
    """Check each family's _get_fsigma, the override the comparison calls."""
    assert issubclass(ccl_family, pyccl.halos.MassFunc)
    parameters = inspect.signature(ccl_family._get_fsigma).parameters

    assert list(parameters) == GET_FSIGMA_PARAMETERS


def test_pairs_cover_every_constructor() -> None:
    """Check that PAIRS builds every *_pair constructor of the pairing module.

    A constructor added to the module must be added to CONSTRUCTOR_CALLS, so that its
    CCL family's _get_fsigma is pinned too.
    """
    module_constructors = {
        name
        for name, value in vars(mfp).items()
        if name.endswith("_pair") and callable(value)
    }
    built = {constructor.__name__ for constructor, _ in CONSTRUCTOR_CALLS}
    missing = sorted(module_constructors - built)

    assert not missing, f"PAIRS does not build: {', '.join(missing)}"


def test_compare_mass_function(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
) -> None:
    """Compare the halo mass function of CCL and NumCosmo with compare_mass_function.

    Checks that M dn/dlog10(M) agrees for Tinker (2008) 200m at z = 0, 0.05, 0.5
    and 1. At CCL defaults the largest difference is at z = 0.05, where CCL's
    sigma(M, a) table (A_SPLINE_NA_SM = 13, A_SPLINE_NLOG_SM = 6) is interpolated
    between its last two scale-factor nodes, a = 0.925 and a = 1.
    """
    m = np.geomspace(1.0e11, 1.0e15, 64)
    reltol_target = 1.0e-6 if ccl_cosmo_eh_linear.high_precision else 5.0e-3

    for z in (0.0, 0.05, 0.5, 1.0):
        ccl_hmf, nc_mulf = _tinker08_200m()
        cmp = nc_cmp.compare_mass_function(
            ccl_cosmo_eh_linear,
            nc_cosmo_eh_linear,
            ccl_hmf,
            nc_mulf,
            m,
            z,
            model="Tinker08 200m",
        )

        assert_allclose(cmp.y1, cmp.y2, rtol=reltol_target)


def test_compare_mass_function_reporting(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
) -> None:
    """Check that compare_mass_function reports through CompareFunc1d.

    The benchmark pages read the returned CompareFunc1d, so it must carry the
    masses as its abscissa, positive mass functions from both codes, and a summary
    row matching the shared table header.
    """
    m = np.geomspace(1.0e11, 1.0e15, 32)
    ccl_hmf, nc_mulf = _tinker08_200m()

    cmp = nc_cmp.compare_mass_function(
        ccl_cosmo_eh_linear,
        nc_cosmo_eh_linear,
        ccl_hmf,
        nc_mulf,
        m,
        0.0,
        model="Tinker08 200m",
    )

    assert isinstance(cmp, nc_cmp.CompareFunc1d)
    assert_allclose(cmp.x, m, rtol=0.0, atol=0.0)
    assert cmp.y1.shape == m.shape
    assert cmp.y2.shape == m.shape
    assert np.all(cmp.y1 > 0.0)
    assert np.all(cmp.y2 > 0.0)
    assert cmp.model == "Tinker08 200m"
    assert cmp.rel_diff_max >= cmp.rel_diff_min
    assert len(cmp.summary_row()) == len(nc_cmp.CompareFunc1d.table_header())


@pytest.mark.parametrize(
    "compare_func",
    [nc_cmp.compare_mass_function, nc_cmp.compare_multiplicity_function],
    ids=["mass_function", "multiplicity_function"],
)
def test_compare_rejects_redshift_outside_filter(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    compare_func: Callable[..., nc_cmp.CompareFunc1d],
) -> None:
    """Check that both comparisons raise ValueError for z above the filter's zf."""
    z = nc_cosmo_eh_linear.psf_tophat.props.zf + 1.0
    ccl_hmf, nc_mulf = _tinker08_200m()

    with pytest.raises(ValueError, match="outside the filter's redshift range"):
        compare_func(
            ccl_cosmo_eh_linear,
            nc_cosmo_eh_linear,
            ccl_hmf,
            nc_mulf,
            np.geomspace(1.0e11, 1.0e15, 4),
            z,
        )


@pytest.mark.parametrize("family, constructor, args", MULTIPLICITY_CONFIGS)
def test_compare_multiplicity_function(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    family: str,
    constructor: PairConstructor,
    args: tuple[Any, ...],
) -> None:
    """Compare the multiplicity function f(sigma) of CCL and NumCosmo.

    Both codes evaluate f at the same sigma, lnM and z, so differences in
    sigma(M, z) do not enter. Where the family uses them, each code computes its own
    delta_c, Delta and Omega_m(a). Each family's bound is FSIGMA_RELTOL's.
    """
    m = np.geomspace(1.0e11, 1.0e15, 64)

    for z in FSIGMA_Z:
        ccl_hmf, nc_mulf = constructor(*args)
        cmp = nc_cmp.compare_multiplicity_function(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ccl_hmf, nc_mulf, m, z
        )

        assert_allclose(cmp.y1, cmp.y2, rtol=FSIGMA_RELTOL[family], atol=0.0)


@pytest.mark.parametrize(
    "constructor, args",
    CONSTRUCTOR_CALLS,
    ids=[constructor.__name__ for constructor, _ in CONSTRUCTOR_CALLS],
)
def test_compare_multiplicity_function_ccl_inputs(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
    constructor: PairConstructor,
    args: tuple[Any, ...],
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """Check that CCL receives NumCosmo's sigma, the scale factor and lnM.

    CCL's _get_fsigma is wrapped to record its arguments and result.
    """
    m = np.geomspace(1.0e11, 1.0e15, 64)

    # CCL objects are immutable, so the method is wrapped on the family's class.
    calls: list[tuple[Any, ...]] = []
    ccl_family = type(constructor(*args)[0])
    get_fsigma = ccl_family._get_fsigma

    def recording_get_fsigma(
        self: pyccl.halos.MassFunc,
        cosmo_ccl: pyccl.Cosmology,
        sigM: np.ndarray,
        a: float,
        lnM: np.ndarray,
    ) -> np.ndarray:
        fsigma = get_fsigma(self, cosmo_ccl, sigM, a, lnM)
        calls.append((cosmo_ccl, np.copy(sigM), a, np.copy(lnM), np.copy(fsigma)))
        return fsigma

    monkeypatch.setattr(ccl_family, "_get_fsigma", recording_get_fsigma)

    for z in FSIGMA_Z:
        ccl_hmf, nc_mulf = constructor(*args)
        calls.clear()
        cmp = nc_cmp.compare_multiplicity_function(
            ccl_cosmo_eh_linear, nc_cosmo_eh_linear, ccl_hmf, nc_mulf, m, z
        )

        assert len(calls) == 1
        cosmo_ccl, sigM, a, lnM, fsigma = calls[0]
        assert cosmo_ccl is ccl_cosmo_eh_linear
        assert_array_equal(sigM, cmp.x)
        assert a == 1.0 / (1.0 + z)
        assert_array_equal(lnM, np.log(m))
        assert_array_equal(cmp.y1, fsigma)


def test_compare_multiplicity_function_reporting(
    ccl_cosmo_eh_linear: pyccl.Cosmology,
    nc_cosmo_eh_linear: Cosmology,
) -> None:
    """Check that compare_multiplicity_function reports through CompareFunc1d.

    The returned object must carry sigma as its abscissa, positive multiplicity
    functions from both codes, and a summary row matching the shared table header.
    """
    m = np.geomspace(1.0e11, 1.0e15, 32)
    ccl_hmf, nc_mulf = _tinker08_200m()

    cmp = nc_cmp.compare_multiplicity_function(
        ccl_cosmo_eh_linear,
        nc_cosmo_eh_linear,
        ccl_hmf,
        nc_mulf,
        m,
        0.0,
        model="Tinker08 200m",
    )

    assert isinstance(cmp, nc_cmp.CompareFunc1d)
    assert cmp.x.shape == m.shape
    assert cmp.y1.shape == m.shape
    assert cmp.y2.shape == m.shape
    assert np.all(cmp.x > 0.0)
    assert np.all(np.diff(cmp.x) < 0.0)
    assert np.all(cmp.y1 > 0.0)
    assert np.all(cmp.y2 > 0.0)
    assert cmp.model == "Tinker08 200m"
    assert cmp.rel_diff_max >= cmp.rel_diff_min
    assert len(cmp.summary_row()) == len(nc_cmp.CompareFunc1d.table_header())
