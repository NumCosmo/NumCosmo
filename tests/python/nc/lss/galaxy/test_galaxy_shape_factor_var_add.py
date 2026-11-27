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

``NcGalaxyShapeFactor`` evaluates and draws the observed ellipticity of a galaxy at a
given reduced shear; subclasses supply only the intrinsic-ellipticity
marginalization. The shear of a galaxy comes from ``NcGalaxyHaloWL`` and
``NcWLClusterLens`` (``_galaxy_shear`` below). ``NcGalaxyShapeFactorVarAdd`` is the
legacy variance-add Gaussian approximation, validated for parity
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
rather than inline literals; see ``_load_integ_parity_truth_table`` and
``_load_gen_parity_truth_table``.
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
_INTEG_PARITY_TRUTH_TABLE_FILE = (
    "truth_tables/wl/nc_galaxy_shape_factor_var_add_integ_parity.bin"
)


def _load_integ_parity_truth_table() -> np.ndarray:
    """Load the frozen integ() sequences as a (len(_CONVS), len(_GALAXIES),
    2, 100) array."""
    path = Ncm.cfg_get_data_filename(_INTEG_PARITY_TRUTH_TABLE_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(len(_CONVS), len(_GALAXIES), 2, 100)


# Frozen gen() draw sequences, one row per (ellip_conv, galaxy) case, 50
# (eps_int_1, eps_int_2, eps_obs_1, eps_obs_2) tuples per row. Stored as a
# flat (len(_CONVS) * len(_GALAXIES) * 50, 4) matrix, blocked by ellip_conv
# (matching _CONVS order), then by galaxy (matching _GALAXIES order).
# The eps_obs columns were regenerated from the current engine (2026-09-25)
# after ncm_trivec_get_spherical_coord() moved from acos to atan2 and
# ncm_quaternion_set_to_rotate_to_z() to atan2 half-angles: they moved by at
# most 1.8e-12 absolute; the eps_int columns are unchanged.
_GEN_PARITY_TRUTH_TABLE_FILE = (
    "truth_tables/wl/nc_galaxy_shape_factor_var_add_gen_parity.bin"
)


def _load_gen_parity_truth_table() -> np.ndarray:
    """Load the frozen gen() draw sequences as a (len(_CONVS),
    len(_GALAXIES), 50, 4) array."""
    path = Ncm.cfg_get_data_filename(_GEN_PARITY_TRUTH_TABLE_FILE, True)
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
    data = Nc.GalaxyShapeFactorData.new(gsf, mset)
    return gsf, data, pos_data, z_data


def _galaxy_shear(gsf, mset, pos_data, data):
    """The reduced shear of a galaxy, at position pos_data and with shape data data,
    as a function of its source redshift,
    in the frame of data, as NcDataClusterWLFactor builds it from NcGalaxyHaloWL and
    NcWLClusterLens: g = (1 + m) g_t(z) exp(2 i phi) + c. In front of the lens g_t is
    that just behind it times exp((z - z_cl) / 0.001). Also prepares gsf."""
    smd = mset.peek(Nc.WLSurfaceMassDensity.id())
    cosmo = mset.peek(Nc.HICosmo.id())
    hwl = Nc.GalaxyHaloWL.new(Nc.WLClusterLens.new(smd.props.distance))
    h_data = Nc.GalaxyHaloWLData.new()
    lens = hwl.peek_lens()

    gsf.prepare_data_array(mset, [data])
    hwl.prepare(mset)
    hwl.update_data_radius(h_data, pos_data)
    hwl.update_data_profile(mset.peek(Nc.HaloDensityProfile.id()), h_data)
    sigma, mean_sigma = h_data.get_sigma()
    rot_re, rot_im = h_data.get_rotation()
    if data.coord == Nc.WLEllipticityFrame.CARTESIAN:
        rot_im = -rot_im
    z_cl = lens.get_z_cl()
    _, _, _, c1, c2, m = gsf.data_get(data)

    def shear(z):
        if z > z_cl:
            gt = lens.reduced_shear_at(cosmo, sigma, mean_sigma, z)
        else:
            gt0 = lens.reduced_shear_at(
                cosmo, sigma, mean_sigma, z_cl * (1.0 + np.finfo(float).eps)
            )
            gt = math.exp((z - z_cl) / 0.001) * gt0
        return (1.0 + m) * gt * rot_re + c1, (1.0 + m) * gt * rot_im + c2

    return shear, z_cl


def _shape_at_z(gsf, mset, pos_data, data, use_lnp):
    """The shape factor of the galaxy of data as a function of the source redshift:
    eval_shear at the shear of _galaxy_shear."""
    shear, _ = _galaxy_shear(gsf, mset, pos_data, data)
    pop = mset.peek(Nc.GalaxyShapePop.id())

    def f(z):
        g_1, g_2 = shear(z)
        return gsf.eval_shear(pop, data, g_1, g_2, use_lnp)

    return f


def _gen(gsf, mset, pos_data, z_data, data, rng):
    """Draws the observed ellipticity of the galaxy at pos_data and z_data."""
    shear, z_cl = _galaxy_shear(gsf, mset, pos_data, data)
    z = z_data.z
    lensed = z > z_cl
    if lensed:
        g_1, g_2 = shear(z)
    else:
        _, _, _, g_1, g_2, _ = gsf.data_get(data)
    gsf.gen_shear(mset.peek(Nc.GalaxyShapePop.id()), data, lensed, g_1, g_2, rng)


def _set_galaxy(new, galaxy):
    """Feeds one galaxy's raw values into the new engine's data holders."""
    ra, dec, z, e1, e2, std_noise, c1, c2, m = galaxy
    gsf, s_data, pos_data, z_data = new

    pos_data.ra = ra
    pos_data.dec = dec
    z_data.z = z

    gsf.data_set(s_data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)


_INTEG_PARITY_FROZEN = _load_integ_parity_truth_table()


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_parity_legacy(ellip_conv, galaxy, use_lnp):
    """The variance-add integrand is checked against the frozen legacy
    oracle (rtol=1e-12: the integrand calls exp()/complex division, not
    bit-exact across platforms -- see module docstring)."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, pos_data, _ = new

    _set_galaxy(new, galaxy)

    new_integ = _shape_at_z(gsf, mset, pos_data, s_data, use_lnp)

    frozen = _INTEG_PARITY_FROZEN[
        _CONVS.index(ellip_conv), _GALAXIES.index(galaxy), int(use_lnp)
    ]
    for z, expected in zip(np.linspace(0.05, 1.5, 100), frozen):
        assert_allclose(new_integ(z), expected, rtol=1.0e-12, atol=1.0e-12)


_GEN_PARITY_FROZEN = _load_gen_parity_truth_table()


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_gen_parity_legacy(ellip_conv, galaxy):
    """Same seed => identical intrinsic draws and observed ellipticities,
    checked against a frozen legacy draw sequence (rtol=1e-12: gen() draws
    intrinsic noise via ncm_rng_gaussian_gen(), not bit-exact across
    platforms -- see module docstring)."""
    mset = _build_mset()
    new = _build_new(mset, ellip_conv)
    gsf, s_data, pos_data, z_data = new

    _set_galaxy(new, galaxy)

    rng_new = Ncm.RNG.seeded_new(None, 42)

    frozen = _GEN_PARITY_FROZEN[_CONVS.index(ellip_conv), _GALAXIES.index(galaxy)]
    for eps_int_1, eps_int_2, eps_obs_1, eps_obs_2 in frozen:
        _gen(gsf, mset, pos_data, z_data, s_data, rng_new)
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


# gt and gx were recomputed (2026-11-27) when direct_estimate began rotating c1, c2 into the
# tangential/cross frame; sigma_t, sigma_x and rho are the legacy values.
_DIRECT_ESTIMATE_FROZEN = {
    "TRACE": (
        0.0024560963847448133,
        -0.03980888102683978,
        0.18446442899466103,
        0.1677352221394491,
        0.20442287348760513,
    ),
    "TRACE_DET": (
        -0.0002663131398755384,
        -0.07256138173382227,
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
    rotations = []
    hwl = Nc.GalaxyHaloWL.new(
        Nc.WLClusterLens.new(mset.peek(Nc.WLSurfaceMassDensity.id()).props.distance)
    )
    hwl.prepare(mset)

    for galaxy in _GALAXIES * 10:
        new_i = _build_new(mset, ellip_conv)
        _set_galaxy(new_i, galaxy)

        _gen(gsf, mset, new_i[2], new_i[3], new_i[1], rng_new)
        new_array.append(new_i[1])

        h_data = Nc.GalaxyHaloWLData.new()
        hwl.update_data_radius(h_data, new_i[2])
        rotations.extend(h_data.get_rotation())

    new_est = gsf.direct_estimate(mset, new_array, rotations)
    frozen = _DIRECT_ESTIMATE_FROZEN[_CONV_NAMES[ellip_conv]]
    # gt/gx are near-zero cancellation sums over many intrinsic-noise draws,
    # so a sub-ULP cross-platform summation-order difference can amplify a
    # pure rtol check even though it's numerically negligible (see
    # test_gen_parity_legacy above for the same rationale).
    assert_allclose(
        (new_est.gt, new_est.gx, new_est.sigma_t, new_est.sigma_x, new_est.rho),
        frozen,
        rtol=1.0e-12,
        atol=1.0e-12,
    )


def test_required_columns():
    """The factor requires its own columns, then those of the population fragment; the
    position and redshift columns belong to their own fragments."""
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
    for col in ("ra", "dec", "zp"):
        assert col not in cols


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
