#!/usr/bin/env python
#
# test_galaxy_sd_redshift_observable_population.py
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

"""Tests for the population-level photo-z observable distribution.

``NcGalaxyRedshiftObsSel`` exposes the pointwise density
``eval(z, obs)`` and its window integral ``window_mass(z, lo, hi)``; the two must
be consistent (window_mass = int eval dobs). The Gaussian subclass matches
N(z, sigma0*(1+z)).
"""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Nc

_CASES = [(0.05, 0.5), (0.03, 1.2), (0.1, 2.0)]


def _build(sigma0):
    obs_pop = Nc.GalaxyRedshiftObsSelGauss.new()
    obs_pop.param_set_by_name("sigma0", sigma0)
    return obs_pop


@pytest.mark.parametrize("sigma0,z", _CASES)
def test_eval_matches_gaussian(sigma0, z):
    """eval(z, obs) is N(obs | z, sigma0*(1+z))."""
    obs_pop = _build(sigma0)
    sigmaz = sigma0 * (1.0 + z)
    obs = np.linspace(z - 4.0 * sigmaz, z + 4.0 * sigmaz, 200)
    got = np.array([obs_pop.eval(z, o) for o in obs])
    expected = np.exp(-0.5 * ((obs - z) / sigmaz) ** 2) / (
        np.sqrt(2.0 * np.pi) * sigmaz
    )
    assert_allclose(got, expected, rtol=1.0e-12)


@pytest.mark.parametrize("sigma0,z", _CASES)
def test_eval_integrates_to_window_mass(sigma0, z):
    """int_lo^hi eval(z, obs) dobs equals window_mass(z, lo, hi)."""
    obs_pop = _build(sigma0)
    sigmaz = sigma0 * (1.0 + z)
    lo, hi = z - 1.5 * sigmaz, z + 2.5 * sigmaz
    obs = np.linspace(lo, hi, 40000)
    dens = np.array([obs_pop.eval(z, o) for o in obs])
    quad = np.trapezoid(dens, obs)
    assert_allclose(quad, obs_pop.window_mass(z, lo, hi), rtol=1.0e-5)


def test_full_window_mass_is_physical_norm():
    """window_mass(z, 0, inf) is the zp >= 0 normalization 1/2 erfc(-z/sqrt2 sigmaz)."""
    from scipy.special import erfc

    obs_pop = _build(0.05)
    for z in (0.2, 1.0, 2.5):
        sigmaz = 0.05 * (1.0 + z)
        expected = 0.5 * erfc(-z / (np.sqrt(2.0) * sigmaz))
        assert_allclose(obs_pop.window_mass(z, 0.0, np.inf), expected, rtol=1.0e-12)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
