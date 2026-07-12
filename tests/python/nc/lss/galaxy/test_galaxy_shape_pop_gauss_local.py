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
"""

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
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    return gsf, data


def _build_legacy_local(ellip_conv):
    """Legacy per-galaxy HSMGauss (non-Global) oracle."""
    gsdtr = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    gsdor = Nc.GalaxySDObsRedshiftSpec.new(gsdtr, 0.0, 2.0)
    sdz_data = Nc.GalaxySDObsRedshiftData.new(gsdor)
    lpos = Nc.GalaxySDPositionFlat.new(-0.2, 0.2, -0.2, 0.2)
    sdpos_data = Nc.GalaxySDPositionData.new(lpos, sdz_data)

    lshape = Nc.GalaxySDShapeHSMGauss.new(ellip_conv)
    data = Nc.GalaxySDShapeData.new(lshape, sdpos_data)
    return lshape, data, sdpos_data, sdz_data


@pytest.mark.parametrize("ellip_conv", _CONVS)
def test_local_matches_global_at_equal_e_rms(ellip_conv):
    """VarAdd + Local reproduces VarAdd + Global when tuned to the same e_rms."""
    target_e_rms = 0.28

    pop_global = Nc.GalaxyShapePopGauss.new()
    probe_data = Nc.GalaxyShapePopData.new(pop_global)
    _sigma_for_e_rms(pop_global, probe_data, target_e_rms)

    mset_global = _build_mset(pop_global)
    mset_local = _build_mset(Nc.GalaxyShapePopGaussLocal.new())

    gsf_g, data_g = _build_new(mset_global, ellip_conv)
    gsf_l, data_l = _build_new(mset_local, ellip_conv)

    for gsf, data, mset in ((gsf_g, data_g, mset_global), (gsf_l, data_l, mset_local)):
        pos_data = data.pos_data
        z_data = data.z_data
        pos_data.ra, pos_data.dec = 0.03, 0.02
        z_data.z = 0.9
        gsf.data_set(data, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05, Nc.WLEllipticityFrame.CELESTIAL)

    pop_local = mset_local.peek(Nc.GalaxyShapePop.id())
    pop_local_data = data_l.pop_data
    Nc.GalaxyShapePopGaussLocal.data_set(pop_local, pop_local_data, target_e_rms)

    gsf_g.prepare_data_array(mset_global, [data_g], True, True)
    gsf_l.prepare_data_array(mset_local, [data_l], True, True)

    integ_g = gsf_g.integ(mset_global, False)
    integ_l = gsf_l.integ(mset_local, False)

    for z in np.linspace(0.05, 1.5, 50):
        assert_allclose(integ_l.eval(z, data_l), integ_g.eval(z, data_g), rtol=1.0e-8)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("galaxy", _GALAXIES)
@pytest.mark.parametrize("use_lnp", [False, True])
def test_integ_parity_legacy_per_galaxy(ellip_conv, galaxy, use_lnp):
    """VarAdd + Local matches the legacy per-galaxy HSMGauss oracle.

    Not bit-identical: the sigma-from-e_rms inversion uses an independent
    bisection here rather than legacy's Newton solve on a differently
    rearranged (but algebraically equivalent) formula.
    """
    ra, dec, z, e1, e2, std_shape, std_noise, c1, c2, m = galaxy

    mset = _build_mset(Nc.GalaxyShapePopGaussLocal.new())
    gsf, data = _build_new(mset, ellip_conv)
    lshape, ls_data, sdpos_data, sdz_data = _build_legacy_local(ellip_conv)

    data.pos_data.ra, data.pos_data.dec = ra, dec
    data.z_data.z = z
    sdpos_data.ra, sdpos_data.dec = ra, dec
    sdz_data.z = z

    gsf.data_set(data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)
    lshape.data_set(ls_data, e1, e2, std_shape, std_noise, c1, c2, m)

    pop = mset.peek(Nc.GalaxyShapePop.id())
    Nc.GalaxyShapePopGaussLocal.data_set(pop, data.pop_data, std_shape)

    gsf.prepare_data_array(mset, [data], True, True)
    lshape.prepare_data_array(mset, [ls_data], True, True)

    new_integ = gsf.integ(mset, use_lnp)
    old_integ = lshape.integ(use_lnp)
    old_integ.prepare(mset)

    for zz in np.linspace(0.05, 1.5, 50):
        assert_allclose(new_integ.eval(zz, data), old_integ.eval(zz, ls_data), rtol=1.0e-8)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
