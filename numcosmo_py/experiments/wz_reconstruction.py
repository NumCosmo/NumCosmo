#
# wz_reconstruction.py
#
# Sat Jun 21 2026
# Copyright  2026  Sandro Dias Pinto Vitenti <vitenti@uel.br>
#
# wz_reconstruction.py
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

"""Function-space tools for the q(z)/w(z) curvature reconstruction study.

This module defines a *systematic* function space for injected truths, replacing
the hand-picked fiducial models of the previous q(z) reconstruction. Deviations
from a fiducial are expanded on a truncated orthogonal basis with a
spectral-decay prior on the coefficients,

    f(x) = offset + sum_{k=1}^{K} a_k phi_k(x),   a_k ~ N(0, sigma_k^2),

with sigma_k = amplitude * k^{-decay}. The amplitude sets the overall feature
size and the decay sets the smoothness / feature scale (small decay -> more
high-frequency power -> sharper, higher-curvature features). The induced
curvature is measured with the same NcmSpline curvature backend used by the
reconstruction priors, so injected and recovered curvature live on a common
footing.
"""

from enum import StrEnum, auto
import dataclasses

import numpy as np
import numpy.typing as npt

from numcosmo_py import Ncm


class BasisType(StrEnum):
    """Orthogonal basis for the truncated function-space expansion."""

    CHEBYSHEV = auto()
    FOURIER = auto()


@dataclasses.dataclass(kw_only=True)
class TruncatedBasisSampler:
    """Sampler of functions from a truncated-basis, spectral-decay prior.

    The domain is ``[x_min, x_max]`` in the reconstruction's natural variable
    (``alpha = ln(1+z)`` for w(z), ``z`` for q(z)). Coefficients are drawn as
    ``a_k ~ N(0, (amplitude * k^{-decay})^2)`` for ``k = 1 .. n_modes``.
    """

    x_min: float
    x_max: float
    n_modes: int
    amplitude: float
    decay: float
    basis: BasisType = BasisType.CHEBYSHEV
    offset: float = 0.0
    curvature_grid: int = 1025

    def __post_init__(self) -> None:
        if self.x_max <= self.x_min:
            raise ValueError("x_max must be greater than x_min.")
        if self.n_modes < 1:
            raise ValueError("n_modes must be >= 1.")
        if self.amplitude < 0.0:
            raise ValueError("amplitude must be non-negative.")

    def mode_sigma(self) -> npt.NDArray[np.float64]:
        """Per-mode prior standard deviations ``sigma_k``."""
        k = np.arange(1, self.n_modes + 1, dtype=float)
        return self.amplitude * k ** (-self.decay)

    def design_matrix(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Basis matrix ``phi_k(x)`` with shape ``(len(x), n_modes)``."""
        x = np.asarray(x, dtype=float)
        if self.basis is BasisType.CHEBYSHEV:
            u = 2.0 * (x - self.x_min) / (self.x_max - self.x_min) - 1.0
            # chebvander gives T_0..T_{n_modes}; drop the constant T_0.
            return np.polynomial.chebyshev.chebvander(u, self.n_modes)[:, 1:]
        # Fourier: sine modes vanishing at the endpoints.
        t = (x - self.x_min) / (self.x_max - self.x_min)
        k = np.arange(1, self.n_modes + 1, dtype=float)
        return np.sin(np.pi * np.outer(t, k))

    def sample_coeffs(self, rng: np.random.Generator) -> npt.NDArray[np.float64]:
        """Draw a coefficient vector from the spectral-decay prior."""
        return rng.normal(0.0, self.mode_sigma())

    def evaluate(
        self,
        coeffs: npt.NDArray[np.float64],
        x: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate ``offset + sum_k a_k phi_k(x)`` at ``x``."""
        coeffs = np.asarray(coeffs, dtype=float)
        if coeffs.shape != (self.n_modes,):
            raise ValueError(f"coeffs must have shape ({self.n_modes},).")
        return self.offset + self.design_matrix(x) @ coeffs

    def as_spline(self, coeffs: npt.NDArray[np.float64]) -> Ncm.Spline:
        """Build a fine not-a-knot spline of the drawn function for curvature."""
        x = np.linspace(self.x_min, self.x_max, self.curvature_grid)
        y = self.evaluate(coeffs, x)
        spline = Ncm.SplineCubicNotaknot.new_full(
            Ncm.Vector.new_array(x.tolist()),
            Ncm.Vector.new_array(y.tolist()),
            True,
        )
        return spline

    def curvature_lp(
        self,
        coeffs: npt.NDArray[np.float64],
        ctype: Ncm.SplineCurvatureType,
        p: float,
    ) -> float:
        """Domain-normalized ``L_p`` curvature norm of the drawn function.

        Uses the same NcmSpline backend as the reconstruction priors, so this is
        the injected-curvature counterpart of the recovered statistic.
        """
        return self.as_spline(coeffs).curvature_lp_norm(
            ctype, p, self.x_min, self.x_max
        )
