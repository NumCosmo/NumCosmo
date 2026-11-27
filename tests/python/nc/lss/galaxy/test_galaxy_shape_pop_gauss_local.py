#!/usr/bin/env python
#
# test_galaxy_shape_pop_gauss_local.py
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

"""Global vs per-object intrinsic variance for the truncated-Gaussian Pop.

``NcGalaxyShapePopGauss`` ("Global") draws its width sigma from a single
shared model parameter; ``NcGalaxyShapePopGaussLocal`` resolves it per galaxy
from an input RMS ellipticity carried on the Data fragment (as legacy
``NcGalaxySDShapeHSMGauss`` reads a per-galaxy ``std_shape`` catalog column,
contrasted with ``NcGalaxySDShapeHSMGaussGlobal``'s single ``sigma``
parameter). These tests check: (1) Global and Local agree when tuned to the
same population width, and (2) ``NcGalaxyShapeFactorVarAdd`` built on Local
reproduces the legacy per-galaxy oracle.

FROZEN REFERENCE VALUES: the parity documented in (2) above was proven by
running both engines live and is captured, not re-derived, in
``test_integ_parity_legacy_per_galaxy`` below. Values were captured from an
actual passing run of this file's original legacy-comparison code, at git
rev ``77313f22`` (2026-07-16), then legacy (``NcGalaxySDShapeHSMGauss``/
``NcGalaxySDObsRedshiftSpec``/``NcGalaxySDPositionFlat``/
``NcGalaxySDTrueRedshiftLSSTSRD``) construction was removed so this test no
longer depends on legacy at runtime -- legacy is slated for deletion in a
follow-up PR. The comparison keeps the tolerance (``rtol=1e-8``) that the
original live comparison used -- not bit-identical, since the sigma-from-
e_rms inversion uses an independent bisection here rather than legacy's
Newton solve on a differently rearranged (but algebraically equivalent)
formula. The captured sequence is stored as an ``Ncm.Matrix`` binfile
(``data/truth_tables/wl/``) rather than inline literals; see ``_load_truth_table``.
"""

import math
import subprocess
import sys

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]

# (ra, dec, z, epsilon_obs_1, epsilon_obs_2, std_shape, std_noise, c1, c2, m)
_GALAXIES = [
    (0.03, 0.02, 0.90, 0.05, -0.02, 0.28, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, 0.45, -0.04, 0.01, 0.35, 0.05, -0.002, 0.004, -0.10),
    (0.05, -0.08, 0.15, 0.02, 0.03, 0.22, 0.04, 0.0, 0.0, 0.0),  # z below cluster
]


def _sigma_for_e_rms(pop, data, target, lo=1.0e-3, hi=5.0, tol=1.0e-13):
    """Bisects the model's own sigma parameter until e_rms(data) == target."""
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        pop.param_set_by_name("sigma", mid)
        pop.prepare(data)
        if pop.e_rms(data) < target:
            lo = mid
        else:
            hi = mid
        if (hi - lo) < tol * hi:
            break
    return 0.5 * (lo + hi)


def _build_mset(pop):
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    dp = Nc.HaloDensityProfileNFW.new(hms)
    hp = Nc.HaloPosition.new(dist)
    smd = Nc.WLSurfaceMassDensity.new(dist)

    hms.param_set_by_name("log10MDelta", 14.0)
    hp.param_set_by_name("z", 0.2)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (cosmo, dp, hp, smd, pop):
        mset.set(model)

    mset.set(Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source())
    mset.set(Nc.GalaxyRedshiftObsGauss.new())

    return mset


def _build_new(mset, ellip_conv):
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


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_local_matches_global_at_equal_e_rms(ellip_conv):
    """VarAdd + Local reproduces VarAdd + Global when tuned to the same e_rms."""
    target_e_rms = 0.28

    pop_global = Nc.GalaxyShapePopGauss.new()
    probe_data = Nc.GalaxyShapePopData.new(pop_global)
    _sigma_for_e_rms(pop_global, probe_data, target_e_rms)

    mset_global = _build_mset(pop_global)
    mset_local = _build_mset(Nc.GalaxyShapePopGaussLocal.new())

    gsf_g, data_g, _, _ = _build_new(mset_global, ellip_conv)
    gsf_l, data_l, _, _ = _build_new(mset_local, ellip_conv)

    for gsf, data in ((gsf_g, data_g), (gsf_l, data_l)):
        gsf.data_set(
            data, 0.05, -0.02, 0.03, 0.0, 0.0, 0.0, Nc.WLEllipticityFrame.CELESTIAL
        )

    pop_local = mset_local.peek(Nc.GalaxyShapePop.id())
    Nc.GalaxyShapePopGaussLocal.data_set(pop_local, data_l.pop_data, target_e_rms)
    gsf_g.prepare_data_array(mset_global, [data_g])
    gsf_l.prepare_data_array(mset_local, [data_l])
    pop_global_m = mset_global.peek(Nc.GalaxyShapePop.id())

    for g_1, g_2 in ((0.0, 0.0), (0.05, -0.02), (0.2, 0.1), (-0.3, 0.25), (0.6, -0.4)):
        assert_allclose(
            gsf_l.eval_shear(pop_local, data_l, g_1, g_2, False),
            gsf_g.eval_shear(pop_global_m, data_g, g_1, g_2, False),
            rtol=1.0e-8,
        )


_TRUTH_TABLE_FILE = "truth_tables/wl/nc_galaxy_shape_pop_gauss_local_integ_parity.bin"
_TRUTH_TABLE_N = 50


def _load_truth_table() -> np.ndarray:
    """Load the frozen per-galaxy integ() sequences as a
    (len(_CONVS), len(_GALAXIES), 2, _TRUTH_TABLE_N) array, blocked by ellip_conv
    (matching _CONVS order), then by galaxy (matching _GALAXIES order), then
    by use_lnp (False, True)."""
    path = Ncm.cfg_get_data_filename(_TRUTH_TABLE_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    matrix = ser.from_binfile(path)
    assert isinstance(matrix, Ncm.Matrix)
    return np.array(matrix.dup_array()).reshape(
        len(_CONVS), len(_GALAXIES), 2, _TRUTH_TABLE_N
    )


_INTEG_PARITY_PER_GALAXY_FROZEN = _load_truth_table()


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_parity_legacy_per_galaxy(ellip_conv, galaxy, use_lnp):
    """VarAdd + Local is checked against the frozen legacy per-galaxy
    HSMGauss oracle (see module docstring for the tolerance rationale)."""
    ra, dec, z, e1, e2, std_shape, std_noise, c1, c2, m = galaxy

    mset = _build_mset(Nc.GalaxyShapePopGaussLocal.new())
    gsf, data, pos_data, z_data = _build_new(mset, ellip_conv)

    pos_data.ra, pos_data.dec = ra, dec
    z_data.z = z

    gsf.data_set(data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)

    pop = mset.peek(Nc.GalaxyShapePop.id())
    Nc.GalaxyShapePopGaussLocal.data_set(pop, data.pop_data, std_shape)

    new_integ = _shape_at_z(gsf, mset, pos_data, data, use_lnp)

    frozen = _INTEG_PARITY_PER_GALAXY_FROZEN[
        _CONVS.index(ellip_conv), _GALAXIES.index(galaxy), int(use_lnp)
    ]
    zs = np.linspace(0.05, 1.5, 50)
    new_vals = [new_integ(zz) for zz in zs]
    assert_allclose(new_vals, frozen, rtol=1.0e-8)


def test_ldata_read_row_native_e_rms():
    """A catalog with the current 'e_rms' column reads it directly -- no
    legacy fallback, no warning."""
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE,
        Nc.WLEllipticityFrame.CELESTIAL,
        1,
        ["e_rms"],
        [Ncm.CatalogColType.DOUBLE],
    )
    obs.set("e_rms", 0, 0.31)

    pop = Nc.GalaxyShapePopGaussLocal.new()
    data = Nc.GalaxyShapePopData.new(pop)
    data.read_row(obs, 0)

    assert data.e_rms == 0.31


def test_ldata_read_row_legacy_std_shape_fallback():
    """A catalog with only the pre-rename 'std_shape' column (no 'e_rms')
    falls back to it -- same value, plus a one-time warning naming both
    columns and the cache path to delete for a redownload. Checked in a
    subprocess: the fallback's own warning-once guard is a process-lifetime
    static, so a fresh process is the only way to reliably observe it fire
    (see nc_galaxy_shape_pop_gauss_local.c's
    _nc_galaxy_shape_pop_gauss_local_ldata_read_row)."""
    script = (
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "obs = Nc.GalaxyWLObs.new(\n"
        "    Nc.GalaxyWLObsEllipConv.TRACE, Nc.WLEllipticityFrame.CELESTIAL, 1,\n"
        "    ['std_shape'], [Ncm.CatalogColType.DOUBLE],\n"
        ")\n"
        "obs.set('std_shape', 0, 0.28)\n"
        "pop = Nc.GalaxyShapePopGaussLocal.new()\n"
        "data = Nc.GalaxyShapePopData.new(pop)\n"
        "data.read_row(obs, 0)\n"
        "assert data.e_rms == 0.28, data.e_rms\n"
        "print('OK')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert "OK" in result.stdout
    assert "no 'e_rms' column" in result.stderr
    assert "'std_shape'" in result.stderr


def test_ldata_read_row_no_e_rms_or_legacy_column_errors():
    """Neither 'e_rms' nor the legacy 'std_shape' present: falls through to
    the original (pre-fallback) behaviour -- a fatal, real process abort
    from ncm_catalog_get()'s column lookup, not a silently wrong value.
    Checked in a subprocess as above."""
    script = (
        "from numcosmo_py import Nc, Ncm\n"
        "Ncm.cfg_init()\n"
        "obs = Nc.GalaxyWLObs.new(\n"
        "    Nc.GalaxyWLObsEllipConv.TRACE, Nc.WLEllipticityFrame.CELESTIAL, 1,\n"
        "    ['ra'], [Ncm.CatalogColType.DOUBLE],\n"
        ")\n"
        "obs.set('ra', 0, 0.0)\n"
        "pop = Nc.GalaxyShapePopGaussLocal.new()\n"
        "data = Nc.GalaxyShapePopData.new(pop)\n"
        "data.read_row(obs, 0)\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script],
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode != 0
    assert "Column 'e_rms' not found" in result.stderr


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
