#!/usr/bin/env python
#
# test_galaxy_shape_factor_var_add.py
#
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Tests for the shape factor calculator with the variance-add method.

``NcGalaxyShapeFactor`` owns the whole HSM measurement engine (generation,
geometry caches, frame bookkeeping); subclasses supply only the
intrinsic-ellipticity marginalization. ``NcGalaxyShapeFactorVarAdd`` is the
legacy variance-add Gaussian approximation, validated for golden parity
against the pristine ``NcGalaxySDShapeHSMGaussGlobal`` oracle (same math, the
intrinsic width coming from the ``NcGalaxyShapePopGauss`` model in the mset
instead of a parameter of the shape model itself).

FROZEN REFERENCE VALUES: the parity documented above was proven by running
both engines live and is captured, not re-derived, in
``test_integ_parity_legacy``/``test_gen_parity_legacy``/
``test_direct_estimate_parity_legacy`` below. Values were captured from an
actual passing run of this file's original legacy-comparison code, at git
rev ``77313f22`` (2026-07-16), then legacy
(``NcGalaxySDShapeHSMGaussGlobal``/``NcGalaxySDObsRedshiftSpec``/
``NcGalaxySDPositionFlat``/``NcGalaxySDTrueRedshiftLSSTSRD``) construction
was removed so these tests no longer depend on legacy at runtime -- legacy
is slated for deletion in a follow-up PR. Each frozen assertion keeps the
tolerance (``rtol``/``atol``/bit-exact ``==``) that the original live
comparison used. The ``_INTEG_PARITY_FROZEN`` and ``_GEN_PARITY_FROZEN``
sequences are stored as ``Ncm.Matrix`` binfiles (``data/truth_tables/wl/``)
rather than inline literals; see ``_load_integ_parity_golden`` and
``_load_gen_parity_golden``.
"""

import math

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

SIGMA_INT = 0.3

# (ra, dec, z, epsilon_obs_1, epsilon_obs_2, std_noise, c1, c2, m)
_GALAXIES = [
    (0.03, 0.02, 0.90, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, 0.45, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10),
    (0.05, -0.08, 0.15, 0.02, 0.03, 0.04, 0.0, 0.0, 0.0),  # z below the cluster
]

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]
_CONV_NAMES = {
    Nc.GalaxyWLObsEllipConv.TRACE: "TRACE",
    Nc.GalaxyWLObsEllipConv.TRACE_DET: "TRACE_DET",
}

# Frozen integ() sequences, one row per (ellip_conv, galaxy, use_lnp) case,
# 100 values per row (matching np.linspace(0.05, 1.5, 100)). Stored as a
# flat (len(_CONVS) * len(_GALAXIES) * 2, 100) matrix, blocked by ellip_conv
# (matching _CONVS order), then by galaxy (matching _GALAXIES order), then
# by use_lnp (False, True).
_INTEG_PARITY_GOLDEN_FILE = (
    "truth_tables/wl/nc_galaxy_shape_factor_var_add_integ_parity.bin"
)


def _load_integ_parity_golden() -> np.ndarray:
    """Load the frozen integ() sequences as a (len(_CONVS), len(_GALAXIES),
    2, 100) array."""
    path = Ncm.cfg_get_data_filename(_INTEG_PARITY_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CONVS), len(_GALAXIES), 2, 100)


# Frozen gen() draw sequences, one row per (ellip_conv, galaxy) case, 50
# (eps_int_1, eps_int_2, eps_obs_1, eps_obs_2) tuples per row. Stored as a
# flat (len(_CONVS) * len(_GALAXIES) * 50, 4) matrix, blocked by ellip_conv
# (matching _CONVS order), then by galaxy (matching _GALAXIES order).
_GEN_PARITY_GOLDEN_FILE = (
    "truth_tables/wl/nc_galaxy_shape_factor_var_add_gen_parity.bin"
)


def _load_gen_parity_golden() -> np.ndarray:
    """Load the frozen gen() draw sequences as a (len(_CONVS),
    len(_GALAXIES), 50, 4) array."""
    path = Ncm.cfg_get_data_filename(_GEN_PARITY_GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CONVS), len(_GALAXIES), 50, 4)


def _build_mset():
    """One mset serving the new calculator: lens models plus the population
    model it reads."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()

    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    pop.param_set_by_name("sigma", SIGMA_INT)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)

    # Redshift models needed only to build the composed redshift fragment.
    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    return mset


def _build_new(mset, ellip_conv):
    """VarAdd factor + per-galaxy data wired from position/redshift fragments."""
    posf = Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2)
    pos_data = Nc.GalaxyPositionFactorData.new(posf, mset)
    zf = Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0)
    z_data = Nc.GalaxyRedshiftFactorData.new(zf, mset)

    gsf = Nc.GalaxyShapeFactorVarAdd.new(ellip_conv)
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    return gsf, data, pos_data, z_data


def _set_galaxy(new, galaxy):
    """Feeds one galaxy's raw values into the new engine's data holders."""
    ra, dec, z, e1, e2, std_noise, c1, c2, m = galaxy
    gsf, s_data, pos_data, z_data = new

    pos_data.ra = ra
    pos_data.dec = dec
    z_data.z = z

    gsf.data_set(s_data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)


_INTEG_PARITY_FROZEN = _load_integ_parity_golden()


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_parity_legacy(ellip_conv, galaxy, use_lnp):
    """The variance-add integrand is checked against the frozen legacy
    oracle (rtol=1e-12: the integrand calls exp()/complex division, not
    bit-exact across platforms -- see module docstring)."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new

    _set_galaxy(new, galaxy)

    gsf.prepare_data_array(mset, [s_data], True, True)

    new_integ = gsf.integ(mset, use_lnp)

    frozen = _INTEG_PARITY_FROZEN[
        _CONVS.index(ellip_conv), _GALAXIES.index(galaxy), int(use_lnp)
    ]
    for z, expected in zip(np.linspace(0.05, 1.5, 100), frozen):
        assert_allclose(new_integ.eval(z, s_data), expected, rtol=1.0e-12, atol=1.0e-12)


_GEN_PARITY_FROZEN = _load_gen_parity_golden()


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_gen_parity_legacy(ellip_conv, galaxy):
    """Same seed => identical intrinsic draws and observed ellipticities,
    checked against a frozen legacy draw sequence (rtol=1e-12: gen() draws
    intrinsic noise via ncm_rng_gaussian_gen(), not bit-exact across
    platforms -- see module docstring)."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new

    _set_galaxy(new, galaxy)

    rng_new = Ncm.RNG.seeded_new(None, 42)

    frozen = _GEN_PARITY_FROZEN[_CONVS.index(ellip_conv), _GALAXIES.index(galaxy)]
    for eps_int_1, eps_int_2, eps_obs_1, eps_obs_2 in frozen:
        gsf.gen(mset, s_data, rng_new)
        assert_allclose(
            [
                s_data.epsilon_int_1,
                s_data.epsilon_int_2,
                s_data.epsilon_obs_1,
                s_data.epsilon_obs_2,
            ],
            [eps_int_1, eps_int_2, eps_obs_1, eps_obs_2],
            rtol=1.0e-12,
            atol=1.0e-12,
        )


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_eval_at_nodes_matches_integrand(ellip_conv, galaxy):
    """The fixed-node likelihood path agrees with the integrand path.

    Below the cluster redshift the integrand applies a smooth step
    suppression while the node path uses gt = 0, hence the tolerance.
    """
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new

    _set_galaxy(new, galaxy)

    z_nodes = Ncm.Vector.new_array(np.linspace(0.05, 1.5, 25).tolist())
    out_new = Ncm.Vector.new(z_nodes.len())

    gsf.prepare_data_array_at_nodes(mset, [s_data], [z_nodes], True, True, True)
    gsf.eval_at_nodes(mset, s_data, z_nodes, out_new)

    gsf.prepare_data_array(mset, [s_data], True, True)
    integ = gsf.integ(mset, False)
    expected = [integ.eval(z, s_data) for z in z_nodes.dup_array()]

    assert_allclose(out_new.dup_array(), expected, rtol=1.0e-9)


_DIRECT_ESTIMATE_FROZEN = {
    "TRACE": (
        0.0026892385874255896,
        -0.037018753682268174,
        0.18446442899466103,
        0.1677352221394491,
        0.20442287348760513,
    ),
    "TRACE_DET": (
        -3.317093720609368e-05,
        -0.06977125438926611,
        0.3354157805889639,
        0.30585159452347266,
        0.20456109871796177,
    ),
}


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_direct_estimate_parity_legacy(ellip_conv):
    """Direct shear estimate is checked against the frozen legacy oracle
    (the intrinsic rms is computed from the population model, algebraically
    identical to the legacy std_shape_from_sigma but not bit-identical; see
    module docstring)."""
    mset = _build_mset()
    gsf, _, _, _ = _build_new(mset, ellip_conv)

    rng_new = Ncm.RNG.seeded_new(None, 7)
    new_array = []

    for galaxy in _GALAXIES * 10:
        new_i = _build_new(mset, ellip_conv)
        _set_galaxy(new_i, galaxy)

        gsf.gen(mset, new_i[1], rng_new)
        new_array.append(new_i[1])

    new_est = gsf.direct_estimate(mset, new_array)
    frozen = _DIRECT_ESTIMATE_FROZEN[_CONV_NAMES[ellip_conv]]
    # gt/gx are near-zero cancellation sums over many intrinsic-noise draws,
    # so a sub-ULP cross-platform summation-order difference can blow up a
    # pure rtol check even though it's numerically negligible (see
    # test_gen_parity_legacy above for the same rationale).
    assert_allclose(
        (new_est.gt, new_est.gx, new_est.sigma_t, new_est.sigma_x, new_est.rho),
        frozen,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_ln_marginal_consistency(ellip_conv, galaxy):
    """The ln integrand is the log of the linear integrand."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new
    _set_galaxy(new, galaxy)

    gsf.prepare_data_array(mset, [s_data], True, True)
    lin = gsf.integ(mset, False)
    lnp = gsf.integ(mset, True)

    for z in np.linspace(0.05, 1.5, 20):
        p = lin.eval(z, s_data)
        assert p > 0.0
        assert_allclose(lnp.eval(z, s_data), math.log(p), rtol=1.0e-12)


def test_required_columns():
    """The factor requires its own columns plus the upstream fragments'."""
    mset = _build_mset()
    gsf, s_data, _, _ = _build_new(mset, Nc.GalaxyWLObsEllipConv.TRACE_DET)

    cols = Nc.GalaxyShapeFactorData.required_columns(s_data)
    own = [
        "epsilon_int_1",
        "epsilon_int_2",
        "epsilon_obs_1",
        "epsilon_obs_2",
        "std_noise",
        "c1",
        "c2",
        "m",
    ]
    assert cols[: len(own)] == own
    # Upstream position and redshift fragments append theirs.
    for col in ("ra", "dec", "zp"):
        assert col in cols


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_eval_marginal_dispatch(ellip_conv):
    """The marginal hooks are callable directly and mutually consistent."""
    mset = _build_mset()
    gsf, s_data, _, _ = _build_new(mset, ellip_conv)
    pop = mset.peek(Nc.GalaxyShapePop.id())

    gsf.data_set(
        s_data, 0.04, -0.03, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
    )

    p = gsf.eval_marginal(pop, s_data, 0.02, 0.01, 0.04, -0.03)
    lnp = gsf.eval_ln_marginal(pop, s_data, 0.02, 0.01, 0.04, -0.03)
    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
