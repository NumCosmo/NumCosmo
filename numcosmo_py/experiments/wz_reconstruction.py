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
from abc import ABC, abstractmethod
import dataclasses
import math

import numpy as np
import numpy.typing as npt

from numcosmo_py import Ncm, Nc


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

    def deviation(
        self,
        coeffs: npt.NDArray[np.float64],
        x: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate the deviation ``sum_k a_k phi_k(x)`` (no offset/fiducial)."""
        coeffs = np.asarray(coeffs, dtype=float)
        if coeffs.shape != (self.n_modes,):
            raise ValueError(f"coeffs must have shape ({self.n_modes},).")
        return self.design_matrix(x) @ coeffs

    def evaluate(
        self,
        coeffs: npt.NDArray[np.float64],
        x: npt.NDArray[np.float64],
    ) -> npt.NDArray[np.float64]:
        """Evaluate ``offset + sum_k a_k phi_k(x)`` at ``x``."""
        return self.offset + self.deviation(coeffs, x)

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


@dataclasses.dataclass(frozen=True)
class ProjectionResult:
    """Outcome of projecting a sampled truth onto a reconstruction's knots."""

    injected_curvature: float
    projected_curvature: float
    projection_error: float


class ReconstructionTarget(ABC):
    """A spline reconstruction model fed by the function-space sampler.

    Concrete targets (q(z), w(z)) know their natural variable, fiducial, knot
    placement, and how to read/write knot values, so the sampler and the study
    drivers stay model-agnostic. The sampler supplies the *deviation* from the
    fiducial; the target adds the fiducial and owns the curvature measurement.
    """

    error_grid: int = 1025

    @property
    @abstractmethod
    def domain(self) -> tuple[float, float]:
        """Reconstruction range ``(x_min, x_max)`` in the natural variable."""

    @abstractmethod
    def new_cosmo(self) -> Nc.HICosmo:
        """Build a fresh model instance."""

    @abstractmethod
    def knot_x(self, cosmo: Nc.HICosmo) -> npt.NDArray[np.float64]:
        """Knot positions in the natural variable."""

    @abstractmethod
    def fiducial(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        """Fiducial function values at ``x``."""

    @abstractmethod
    def set_knots(self, cosmo: Nc.HICosmo, values: npt.NDArray[np.float64]) -> None:
        """Set the model's knot values."""

    @abstractmethod
    def evaluate(
        self, cosmo: Nc.HICosmo, x: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        """Evaluate the model's reconstructed function at ``x``."""

    @abstractmethod
    def curvature_lp(
        self, cosmo: Nc.HICosmo, ctype: Ncm.SplineCurvatureType, p: float
    ) -> float:
        """Domain-normalized ``L_p`` curvature of the model's function."""

    def set_truth(
        self,
        cosmo: Nc.HICosmo,
        sampler: TruncatedBasisSampler,
        coeffs: npt.NDArray[np.float64],
    ) -> None:
        """Set knots to ``fiducial + deviation`` for the drawn coefficients."""
        x = self.knot_x(cosmo)
        self.set_knots(cosmo, self.fiducial(x) + sampler.deviation(coeffs, x))

    def _truth_spline(
        self,
        sampler: TruncatedBasisSampler,
        coeffs: npt.NDArray[np.float64],
    ) -> Ncm.Spline:
        x_min, x_max = self.domain
        x = np.linspace(x_min, x_max, self.error_grid)
        y = self.fiducial(x) + sampler.deviation(coeffs, x)
        return Ncm.SplineCubicNotaknot.new_full(
            Ncm.Vector.new_array(x.tolist()),
            Ncm.Vector.new_array(y.tolist()),
            True,
        )

    def project(
        self,
        sampler: TruncatedBasisSampler,
        coeffs: npt.NDArray[np.float64],
        ctype: Ncm.SplineCurvatureType,
        p: float,
    ) -> ProjectionResult:
        """Measure representational (knot) bias for a drawn truth.

        Returns the injected curvature (of the fine ``fiducial + deviation``
        truth), the projected curvature (of the model after its knots are set to
        the truth), and the max absolute difference between the truth and the
        knot-limited model over a fine grid. The projection error is the curvature
        the N-knot model *cannot* represent, independent of any data noise.
        """
        x_min, x_max = self.domain
        x = np.linspace(x_min, x_max, self.error_grid)

        truth = self._truth_spline(sampler, coeffs)
        injected = truth.curvature_lp_norm(ctype, p, x_min, x_max)

        cosmo = self.new_cosmo()
        self.set_truth(cosmo, sampler, coeffs)
        projected = self.curvature_lp(cosmo, ctype, p)

        truth_y = np.array([truth.eval(xi) for xi in x])
        model_y = self.evaluate(cosmo, x)
        error = float(np.max(np.abs(truth_y - model_y)))

        return ProjectionResult(
            injected_curvature=injected,
            projected_curvature=projected,
            projection_error=error,
        )


@dataclasses.dataclass
class WSplineTarget(ReconstructionTarget):
    """w(z) reconstruction target with a constant ``w_fid`` fiducial."""

    n_knots: int
    z_max: float
    w_fid: float = -1.0

    @property
    def domain(self) -> tuple[float, float]:
        return (0.0, math.log1p(self.z_max))

    def new_cosmo(self) -> Nc.HICosmoDEWSpline:
        return Nc.HICosmoDEWSpline.new(self.n_knots, self.z_max)

    def knot_x(self, cosmo: Nc.HICosmo) -> npt.NDArray[np.float64]:
        return np.array(cosmo.get_alpha().dup_array())

    def fiducial(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return np.full_like(np.asarray(x, dtype=float), self.w_fid)

    def set_knots(self, cosmo: Nc.HICosmo, values: npt.NDArray[np.float64]) -> None:
        cosmo.props.w = Ncm.Vector.new_array(np.asarray(values).tolist())

    def evaluate(
        self, cosmo: Nc.HICosmo, x: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        # x is alpha = ln(1+z); the model evaluates w as a function of z.
        return np.array([cosmo.w_de(math.expm1(xi)) for xi in np.asarray(x)])

    def curvature_lp(
        self, cosmo: Nc.HICosmo, ctype: Ncm.SplineCurvatureType, p: float
    ) -> float:
        return cosmo.lp_norm(ctype, p)


@dataclasses.dataclass
class QSplineTarget(ReconstructionTarget):
    """q(z) reconstruction target with a fiducial-cosmology q(z)."""

    n_knots: int
    z_max: float
    fiducial_cosmo: Nc.HICosmo

    @property
    def domain(self) -> tuple[float, float]:
        return (0.0, self.z_max)

    def new_cosmo(self) -> Nc.HICosmoQSpline:
        return Nc.HICosmoQSpline.new(
            Ncm.SplineCubicNotaknot.new(), self.n_knots, self.z_max
        )

    def knot_x(self, cosmo: Nc.HICosmo) -> npt.NDArray[np.float64]:
        return np.array(cosmo.props.spline.peek_xv().dup_array())

    def fiducial(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return np.array([self.fiducial_cosmo.q(float(xi)) for xi in np.asarray(x)])

    def set_knots(self, cosmo: Nc.HICosmo, values: npt.NDArray[np.float64]) -> None:
        cosmo.props.qparam = Ncm.Vector.new_array(np.asarray(values).tolist())

    def evaluate(
        self, cosmo: Nc.HICosmo, x: npt.NDArray[np.float64]
    ) -> npt.NDArray[np.float64]:
        return np.array([cosmo.q(float(xi)) for xi in np.asarray(x)])

    def curvature_lp(
        self, cosmo: Nc.HICosmo, ctype: Ncm.SplineCurvatureType, p: float
    ) -> float:
        return cosmo.lp_norm(ctype, p)


def lcdm_fiducial(omega_c0: float = 0.25, omega_b0: float = 0.05) -> Nc.HICosmo:
    """A flat LCDM fiducial cosmology for q(z) reconstruction targets."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo["Omegak"] = 0.0
    cosmo["Omegac"] = omega_c0
    cosmo["Omegab"] = omega_b0
    cosmo["w"] = -1.0
    return cosmo
