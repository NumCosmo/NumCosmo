#
# test_wz_reconstruction.py
#
# Sat Jun 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# test_wz_reconstruction.py
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
# with this program.  If not, see <http://www.gnu.org/licenses/>.

"""Tests for the truncated-basis function-space sampler."""

import pytest
import numpy as np
from numpy.testing import assert_allclose

from numcosmo_py import Ncm
from numcosmo_py.experiments.wz_reconstruction import (
    TruncatedBasisSampler,
    BasisType,
)

Ncm.cfg_init()

GEO = Ncm.SplineCurvatureType.GEOMETRIC
D2 = Ncm.SplineCurvatureType.D2


def make_sampler(**kwargs) -> TruncatedBasisSampler:
    """A default Chebyshev sampler over [0, ln(1+2.33)]."""
    defaults = dict(
        x_min=0.0,
        x_max=float(np.log1p(2.33)),
        n_modes=8,
        amplitude=0.1,
        decay=2.0,
    )
    defaults.update(kwargs)
    return TruncatedBasisSampler(**defaults)


@pytest.mark.parametrize("basis", [BasisType.CHEBYSHEV, BasisType.FOURIER])
def test_design_matrix_shape(basis: BasisType) -> None:
    """Design matrix has one column per mode."""
    sampler = make_sampler(basis=basis)
    x = np.linspace(sampler.x_min, sampler.x_max, 31)
    phi = sampler.design_matrix(x)
    assert phi.shape == (31, sampler.n_modes)


def test_mode_sigma_decay() -> None:
    """Mode sigmas follow the power-law spectral decay."""
    sampler = make_sampler(amplitude=0.3, decay=1.5, n_modes=5)
    k = np.arange(1, 6)
    assert_allclose(sampler.mode_sigma(), 0.3 * k**-1.5)


def test_sampling_reproducible() -> None:
    """Same seed yields identical coefficients; different seed differs."""
    sampler = make_sampler()
    a1 = sampler.sample_coeffs(np.random.default_rng(7))
    a2 = sampler.sample_coeffs(np.random.default_rng(7))
    a3 = sampler.sample_coeffs(np.random.default_rng(8))
    assert_allclose(a1, a2)
    assert not np.allclose(a1, a3)


def test_evaluate_matches_manual_chebyshev() -> None:
    """evaluate() reproduces an explicit Chebyshev sum with offset."""
    sampler = make_sampler(n_modes=4, offset=-1.0)
    coeffs = np.array([0.1, -0.2, 0.05, 0.03])
    x = np.linspace(sampler.x_min, sampler.x_max, 17)
    u = 2.0 * (x - sampler.x_min) / (sampler.x_max - sampler.x_min) - 1.0
    manual = -1.0 + sum(
        c * np.polynomial.chebyshev.Chebyshev.basis(k + 1)(u)
        for k, c in enumerate(coeffs)
    )
    assert_allclose(sampler.evaluate(coeffs, x), manual)


def test_offset_is_flat_zero_curvature() -> None:
    """A pure offset (zero coefficients) has zero curvature."""
    sampler = make_sampler(offset=-1.0)
    coeffs = np.zeros(sampler.n_modes)
    assert_allclose(sampler.curvature_lp(coeffs, D2, 2.0), 0.0, atol=1e-9)
    assert_allclose(sampler.curvature_lp(coeffs, GEO, 2.0), 0.0, atol=1e-9)


def test_curvature_grows_with_amplitude() -> None:
    """Scaling all coefficients by c scales the D2 curvature by |c|."""
    sampler = make_sampler()
    coeffs = sampler.sample_coeffs(np.random.default_rng(3))
    base = sampler.curvature_lp(coeffs, D2, 2.0)
    scaled = sampler.curvature_lp(2.5 * coeffs, D2, 2.0)
    assert_allclose(scaled, 2.5 * base, rtol=1e-6)


def test_curvature_grows_as_decay_decreases() -> None:
    """Less spectral decay (more high-k power) yields more expected curvature."""
    rng = np.random.default_rng(11)
    n_draws = 64

    def mean_curv(decay: float) -> float:
        sampler = make_sampler(decay=decay, amplitude=0.1)
        return float(
            np.mean(
                [
                    sampler.curvature_lp(sampler.sample_coeffs(rng), D2, 2.0)
                    for _ in range(n_draws)
                ]
            )
        )

    assert mean_curv(3.0) < mean_curv(1.0)


def test_invalid_configuration() -> None:
    """Constructor validates its arguments."""
    with pytest.raises(ValueError):
        make_sampler(x_min=1.0, x_max=0.0)
    with pytest.raises(ValueError):
        make_sampler(n_modes=0)
    with pytest.raises(ValueError):
        make_sampler(amplitude=-1.0)
