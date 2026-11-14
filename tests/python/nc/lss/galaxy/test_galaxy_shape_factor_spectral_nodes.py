#!/usr/bin/env python
#
# test_galaxy_shape_factor_spectral_nodes.py
#
# Copyright (C) 2026 Caio Lima de Oliveira <caiolimadeoliveira@pm.me>
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

"""Tests of the node path of NcGalaxyShapeFactorSpectral.

nc_galaxy_shape_factor_eval_at_nodes() and nc_galaxy_shape_factor_eval_ln_at_nodes()
evaluate the nodes of a galaxy through the eval_nodes hook, which Spectral overrides to
check its interpolant once per galaxy. The values must be those of the per-point hooks:
the integrand path evaluates the same marginal at each redshift, with the shear computed
by another route.
"""

import numpy as np
import pytest
from numpy.testing import assert_allclose

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

Z_CL = 0.2

# (ra, dec, epsilon_obs_1, epsilon_obs_2, std_noise, c1, c2, m)
_GALAXIES = [
    (0.03, 0.02, 0.05, -0.02, 0.03, 0.005, -0.003, 0.05),
    (-0.10, 0.15, -0.04, 0.01, 0.05, -0.002, 0.004, -0.10),
    # Close to a massive halo: the nodes far behind it have |g| > interp-g-max
    (0.002, 0.0, 0.30, -0.20, 0.10, 0.0, 0.0, 0.0),
]

_CONVS = [Nc.GalaxyWLObsEllipConv.TRACE, Nc.GalaxyWLObsEllipConv.TRACE_DET]


def _build(ellip_conv, use_interp):
    """A Spectral factor, its model set and the data of one galaxy."""
    cosmo = Nc.HICosmoDEXcdm.new()
    dist = Nc.Distance.new(100.0)
    hms = Nc.HaloCMParam.new(Nc.HaloMassSummaryMassDef.MEAN, 200.0)
    hp = Nc.HaloPosition.new(dist)
    pop = Nc.GalaxyShapePopGauss.new()

    hms.param_set_by_name("log10MDelta", 15.3)
    hp.param_set_by_name("z", Z_CL)
    pop.param_set_by_name("sigma", 0.3)
    hp.prepare(cosmo)

    mset = Ncm.MSet.empty_new()
    for model in (
        cosmo,
        Nc.HaloDensityProfileNFW.new(hms),
        hp,
        Nc.WLSurfaceMassDensity.new(dist),
        pop,
        Nc.GalaxyRedshiftPopLSSTSRD.new_y1_source(),
        Nc.GalaxyRedshiftObsGauss.new(),
    ):
        mset.set(model)

    gsf = Nc.GalaxyShapeFactorSpectral.new(ellip_conv)
    gsf.set_use_interp(use_interp)
    pos_data = Nc.GalaxyPositionFactorData.new(
        Nc.GalaxyPositionFactorFlat.new(-0.2, 0.2, -0.2, 0.2), mset
    )
    z_data = Nc.GalaxyRedshiftFactorData.new(
        Nc.GalaxyRedshiftFactorComposed.new(0.0, 20.0), mset
    )
    data = Nc.GalaxyShapeFactorData.new(gsf, mset, pos_data, z_data)
    return gsf, mset, data, pos_data


def _set_galaxy(gsf, data, pos_data, galaxy):
    ra, dec, e1, e2, std_noise, c1, c2, m = galaxy
    pos_data.ra = ra
    pos_data.dec = dec
    gsf.data_set(data, e1, e2, std_noise, c1, c2, m, Nc.WLEllipticityFrame.CELESTIAL)


def _eval_nodes(gsf, mset, data, z, use_lnp):
    z_nodes = Ncm.Vector.new_array(list(z))
    out = Ncm.Vector.new(len(z))
    gsf.prepare_data_array_at_nodes(mset, [data], [z_nodes], True, True, True)
    if use_lnp:
        gsf.eval_ln_at_nodes(mset, data, z_nodes, out)
    else:
        gsf.eval_at_nodes(mset, data, z_nodes, out)
    return np.array(out.dup_array())


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("use_interp", [False, True])
@pytest.mark.parametrize("galaxy", _GALAXIES)
def test_eval_at_nodes_matches_integrand(ellip_conv, use_interp, galaxy):
    """Behind the lens every node equals the integrand at its redshift."""
    gsf, mset, data, pos_data = _build(ellip_conv, use_interp)
    _set_galaxy(gsf, data, pos_data, galaxy)

    z = np.linspace(0.25, 3.0, 70)
    nodes = _eval_nodes(gsf, mset, data, z, False)

    gsf.prepare_data_array(mset, [data], True, True)
    integ = gsf.integ(mset, False)
    expected = np.array([integ.eval(z_j, data) for z_j in z])

    assert np.all(np.isfinite(nodes))
    assert_allclose(nodes, expected, rtol=1.0e-9)


@pytest.mark.parametrize("ellip_conv", _CONVS)
@pytest.mark.parametrize("use_interp", [False, True])
def test_eval_ln_at_nodes_is_log(ellip_conv, use_interp):
    """The logarithmic node path is the logarithm of the linear one, also across the
    lens redshift and across interp-g-max."""
    gsf, mset, data, pos_data = _build(ellip_conv, use_interp)
    _set_galaxy(gsf, data, pos_data, _GALAXIES[2])

    z = np.linspace(0.05, 3.0, 150)
    lin = _eval_nodes(gsf, mset, data, z, False)
    ln = _eval_nodes(gsf, mset, data, z, True)

    assert_allclose(ln, np.log(lin), rtol=1.0e-13, atol=1.0e-13)


@pytest.mark.parametrize("use_interp", [False, True])
def test_eval_at_nodes_front_of_lens(use_interp):
    """Every node at or in front of the lens has the value of the shear-free point."""
    gsf, mset, data, pos_data = _build(Nc.GalaxyWLObsEllipConv.TRACE, use_interp)
    _set_galaxy(gsf, data, pos_data, _GALAXIES[0])

    z = np.array([0.05, 0.1, 0.15, Z_CL, 0.4])
    nodes = _eval_nodes(gsf, mset, data, z, False)
    _, _, e1, e2, _, c1, c2, _ = _GALAXIES[0]
    expected = gsf.eval_marginal(
        mset.peek_by_name("NcGalaxyShapePop"), data, c1, c2, e1, e2
    )

    assert_allclose(nodes[:4], expected, rtol=0.0, atol=0.0)
    assert nodes[4] != nodes[0]


def test_eval_at_nodes_g_max_fallback():
    """With the interpolant, the nodes with |g| > interp-g-max are integrated: they equal
    the values without the interpolant exactly, and the other nodes agree with them
    within the tolerance of the interpolant. The galaxy has nodes of both kinds."""
    z = np.linspace(0.25, 3.0, 70)
    vals = []
    for use_interp in (False, True):
        gsf, mset, data, pos_data = _build(Nc.GalaxyWLObsEllipConv.TRACE, use_interp)
        _set_galaxy(gsf, data, pos_data, _GALAXIES[2])
        vals.append(_eval_nodes(gsf, mset, data, z, True))

    direct, interp = vals
    fallback = direct == interp

    assert 0 < np.count_nonzero(fallback) < len(z)
    assert_allclose(interp, direct, rtol=0.0, atol=1.0e-5)
