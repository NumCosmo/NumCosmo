#!/usr/bin/env python
#
# test_galaxy_sd_redshift_observable.py
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

"""Standalone unit tests for the photometric-redshift observable model models.

These exercise ``NcGalaxyRedshiftObs`` (Gaussian variant) in isolation from
any calculator: the conditional ``P(z_phot|z)`` is validated for normalization
over ``z_phot``, the analytic peak, consistency of ``gen()`` with the kernel,
and the per-galaxy fragment read/write round-trip.
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="model")
def fixture_model():
    return Nc.GalaxyRedshiftObsGauss.new()


@pytest.mark.parametrize("z_true,sigma0", [(0.3, 0.05), (0.8, 0.03), (1.5, 0.1)])
def test_normalization(model, z_true, sigma0):
    """Integral of P(z_phot|z) over z_phot is unity at fixed true z."""
    data = Nc.GalaxyRedshiftObsData.new(model)
    sz = sigma0 * (1.0 + z_true)
    zps = np.linspace(z_true - 10.0 * sz, z_true + 10.0 * sz, 8001)
    p = np.empty_like(zps)
    for k, zp in enumerate(zps):
        model.data_set(data, zp, sigma0)
        p[k] = model.eval(data, z_true)
    assert_allclose(np.trapezoid(p, zps), 1.0, rtol=1.0e-6)


@pytest.mark.parametrize("z_true,sigma0", [(0.3, 0.05), (0.8, 0.03)])
def test_peak(model, z_true, sigma0):
    """P(z_phot=z|z) equals the Gaussian peak 1/(sqrt(2 pi) sigma_z)."""
    data = Nc.GalaxyRedshiftObsData.new(model)
    sz = sigma0 * (1.0 + z_true)
    model.data_set(data, z_true, sigma0)
    assert_allclose(model.eval(data, z_true), 1.0 / (np.sqrt(2.0 * np.pi) * sz), rtol=1.0e-12)


@pytest.mark.parametrize("z_true,sigma0", [(0.3, 0.05), (0.8, 0.03), (1.5, 0.1)])
def test_gen_consistency(model, z_true, sigma0):
    """gen() samples z_phot ~ N(z, sigma0 (1+z))."""
    data = Nc.GalaxyRedshiftObsData.new(model)
    sz = sigma0 * (1.0 + z_true)
    model.data_set(data, 0.0, sigma0)
    rng = Ncm.RNG.seeded_new(None, 7)
    n_samples = 300000
    s = np.array([model.gen(data, z_true, rng) for _ in range(n_samples)])
    # ~5 sigma bands on mean and std of the mean.
    assert_allclose(s.mean(), z_true, atol=5.0 * sz / np.sqrt(n_samples))
    assert_allclose(s.std(), sz, rtol=0.02)


@pytest.mark.parametrize("z_true,sigma0", [(0.3, 0.05), (0.8, 0.03), (1.5, 0.1)])
def test_window_mass(model, z_true, sigma0):
    """window_mass equals the trapezoid integral of eval over [lo, hi]."""
    data = Nc.GalaxyRedshiftObsData.new(model)
    sz = sigma0 * (1.0 + z_true)
    lo, hi = z_true - 1.3 * sz, z_true + 0.7 * sz
    # Reference: integrate the kernel itself over the window.
    zps = np.linspace(lo, hi, 6001)
    ref = np.empty_like(zps)
    for k, zp in enumerate(zps):
        model.data_set(data, zp, sigma0)
        ref[k] = model.eval(data, z_true)
    got = model.window_mass(data, z_true, lo, hi)
    assert_allclose(got, np.trapezoid(ref, zps), rtol=1.0e-6)
    # Full real line -> unit mass.
    full = model.window_mass(data, z_true, z_true - 40.0 * sz, z_true + 40.0 * sz)
    assert_allclose(full, 1.0, rtol=1.0e-9)


def test_fragment_roundtrip(model):
    """set -> write_row -> read_row -> get is exact for the observation fragment."""
    data = Nc.GalaxyRedshiftObsData.new(model)
    cols = Nc.GalaxyRedshiftObsData.required_columns(data)
    assert set(cols) == {"zp", "sigma0"}

    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE_DET, Nc.WLEllipticityFrame.CELESTIAL, 1, cols
    )
    model.data_set(data, 0.723, 0.041)
    data.write_row(obs, 0)

    data2 = Nc.GalaxyRedshiftObsData.new(model)
    data2.read_row(obs, 0)
    zp, sigma0 = model.data_get(data2)
    assert_allclose([zp, sigma0], [0.723, 0.041], rtol=0.0, atol=0.0)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
