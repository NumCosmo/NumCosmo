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
legacy variance-add Gaussian approximation, validated here for golden parity
against the pristine ``NcGalaxySDShapeHSMGaussGlobal`` oracle (same math, the
intrinsic width coming from the ``NcGalaxyShapePopGauss`` model in the mset
instead of a parameter of the shape model itself).
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


def _build_mset():
    """One mset serving both engines: lens models shared by construction, the
    population model read only by the new calculator."""
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


def _build_legacy(ellip_conv):
    """Legacy HSMGaussGlobal oracle + its nested data chain."""
    gsdtr = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    gsdor = Nc.GalaxySDObsRedshiftSpec.new(gsdtr, 0.0, 2.0)
    sdz_data = Nc.GalaxySDObsRedshiftData.new(gsdor)
    lpos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)
    sdpos_data = Nc.GalaxySDPositionData.new(lpos, sdz_data)

    lshape = Nc.GalaxySDShapeHSMGaussGlobal.new(ellip_conv)
    lshape.param_set_by_name("sigma", SIGMA_INT)
    data = Nc.GalaxySDShapeData.new(lshape, sdpos_data)
    return lshape, data, sdpos_data, sdz_data


def _set_galaxy(new, legacy, galaxy):
    """Feeds the same galaxy into both engines."""
    ra, dec, z, e1, e2, std_noise, c1, c2, m = galaxy
    gsf, s_data, pos_data, z_data = new
    lshape, ls_data, sdpos_data, sdz_data = legacy

    pos_data.ra = ra
    pos_data.dec = dec
    z_data.z = z
    sdpos_data.ra = ra
    sdpos_data.dec = dec
    sdz_data.z = z

    gsf.data_set(s_data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)
    lshape.data_set(ls_data, e1, e2, std_noise, c1, c2, m)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_parity_legacy(ellip_conv, galaxy, use_lnp):
    """The variance-add integrand reproduces the legacy oracle bit for bit."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    legacy = _build_legacy(ellip_conv)
    gsf, s_data, _, _ = new
    lshape, ls_data, _, _ = legacy

    _set_galaxy(new, legacy, galaxy)

    gsf.prepare_data_array(mset, [s_data], True, True)
    lshape.prepare_data_array(mset, [ls_data], True, True)

    new_integ = gsf.integ(mset, use_lnp)
    old_integ = lshape.integ(use_lnp)
    old_integ.prepare(mset)

    for z in np.linspace(0.05, 1.5, 100):
        assert new_integ.eval(z, s_data) == old_integ.eval(z, ls_data)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_gen_parity_legacy(ellip_conv, galaxy):
    """Same seed => identical intrinsic draws and observed ellipticities."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    legacy = _build_legacy(ellip_conv)
    gsf, s_data, _, _ = new
    lshape, ls_data, _, _ = legacy

    _set_galaxy(new, legacy, galaxy)
    _, _, _, _, std_noise, c1, c2, m = galaxy[1:]

    rng_new = Ncm.RNG.seeded_new(None, 42)
    rng_old = Ncm.RNG.seeded_new(None, 42)

    for _ in range(50):
        gsf.gen(mset, s_data, rng_new)
        lshape.gen(
            mset, ls_data, std_noise, c1, c2, m,
            Nc.WLEllipticityFrame.CELESTIAL, rng_old,
        )
        e1, e2, *_ = lshape.data_get(ls_data)
        assert s_data.epsilon_int_1 == ls_data.epsilon_int_1
        assert s_data.epsilon_int_2 == ls_data.epsilon_int_2
        assert s_data.epsilon_obs_1 == e1
        assert s_data.epsilon_obs_2 == e2


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_eval_at_nodes_matches_integrand(ellip_conv, galaxy):
    """The fixed-node likelihood path agrees with the integrand path.

    The legacy eval_at_nodes is not callable through introspection (its @out
    annotation crashes), so the cross-engine check lives in the integrand
    parity test; here the new crit-cache path is validated against the new
    optzs path. Below the cluster redshift the integrand applies a smooth
    step suppression while the node path uses gt = 0, hence the tolerance.
    """
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new

    _set_galaxy(new, _build_legacy(ellip_conv), galaxy)

    z_nodes = Ncm.Vector.new_array(np.linspace(0.05, 1.5, 25).tolist())
    out_new = Ncm.Vector.new(z_nodes.len())

    gsf.prepare_data_array_at_nodes(mset, [s_data], [z_nodes], True, True, True)
    gsf.eval_at_nodes(mset, s_data, z_nodes, out_new)

    gsf.prepare_data_array(mset, [s_data], True, True)
    integ = gsf.integ(mset, False)
    expected = [integ.eval(z, s_data) for z in z_nodes.dup_array()]

    assert_allclose(out_new.dup_array(), expected, rtol=1.0e-9)


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_direct_estimate_parity_legacy(ellip_conv):
    """Direct shear estimate matches the legacy oracle (the intrinsic rms is
    computed from the population model, algebraically identical to the legacy
    std_shape_from_sigma but not bit-identical)."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    legacy = _build_legacy(ellip_conv)
    gsf, _, _, _ = new
    lshape, _, _, _ = legacy

    rng_new = Ncm.RNG.seeded_new(None, 7)
    rng_old = Ncm.RNG.seeded_new(None, 7)
    new_array = []
    old_array = []

    for galaxy in _GALAXIES * 10:
        new_i = _build_new(mset, ellip_conv)
        legacy_i = (lshape, *_build_legacy(ellip_conv)[1:])
        _set_galaxy(new_i, legacy_i, galaxy)
        _, _, _, _, std_noise, c1, c2, m = galaxy[1:]

        gsf.gen(mset, new_i[1], rng_new)
        lshape.gen(
            mset, legacy_i[1], std_noise, c1, c2, m,
            Nc.WLEllipticityFrame.CELESTIAL, rng_old,
        )
        new_array.append(new_i[1])
        old_array.append(legacy_i[1])

    new_est = gsf.direct_estimate(mset, new_array)
    old_est = lshape.direct_estimate(mset, old_array)
    assert_allclose(new_est, old_est, rtol=1.0e-12)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_ln_marginal_consistency(ellip_conv, galaxy):
    """The ln integrand is the log of the linear integrand."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, _, _ = new
    _set_galaxy(new, _build_legacy(ellip_conv), galaxy)

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
        "epsilon_int_1", "epsilon_int_2", "epsilon_obs_1", "epsilon_obs_2",
        "std_noise", "c1", "c2", "m",
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

    gsf.data_set(s_data, 0.04, -0.03, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL)

    p = gsf.eval_marginal(pop, s_data, 0.02, 0.01, 0.04, -0.03)
    lnp = gsf.eval_ln_marginal(pop, s_data, 0.02, 0.01, 0.04, -0.03)
    assert p > 0.0
    assert_allclose(lnp, math.log(p), rtol=1.0e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
