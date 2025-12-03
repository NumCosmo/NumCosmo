#
# _truncated_stats.py
#
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Truncated normal distribution statistics.

This module provides functions for computing statistics of truncated
normal distributions and inverting them to recover underlying parameters.
"""

import sys
from typing import TypedDict

import numpy as np
from scipy.stats import norm
from scipy.optimize import least_squares, OptimizeResult


class InversionResult(TypedDict, total=False):
    """Result from truncated statistics inversion."""

    success: bool
    message: str
    mu: float | None
    sigma: float | None
    cost: float
    nfev: int
    optimality: float


def mean_lnR_truncated(
    mu: float | np.ndarray, sigma: float | np.ndarray, lnR_cut: float
) -> float | np.ndarray:
    """Calculate the mean of a truncated normal distribution.

    For a normal distribution N(mu, sigma^2) truncated from below at lnR_cut,
    compute the expected value of the truncated distribution.

    :param mu: Mean of the underlying normal distribution
    :param sigma: Standard deviation of the underlying normal distribution
    :param lnR_cut: Lower truncation point
    :return: Mean of the truncated distribution
    """
    A = (lnR_cut - mu) / sigma
    phi = norm.pdf(A)
    S = norm.sf(A)
    # avoid division by zero
    T = phi / np.maximum(S, 1e-300)
    return mu + sigma * T


def std_lnR_truncated(
    mu: float | np.ndarray, sigma: float | np.ndarray, lnR_cut: float
) -> float | np.ndarray:
    """Calculate the standard deviation of a truncated normal distribution.

    For a normal distribution N(mu, sigma^2) truncated from below at lnR_cut,
    compute the standard deviation of the truncated distribution.

    :param mu: Mean of the underlying normal distribution
    :param sigma: Standard deviation of the underlying normal distribution
    :param lnR_cut: Lower truncation point
    :return: Standard deviation of the truncated distribution
    """
    A = (lnR_cut - mu) / sigma
    phi = norm.pdf(A)
    S = norm.sf(A)
    T = phi / np.maximum(S, 1e-300)
    corr = np.sqrt(np.maximum(1.0 + A * T - T * T, 0.0))
    return sigma * corr


def invert_truncated_stats(
    m_obs: float,
    s_obs: float,
    lnR_cut: float,
    mu_bounds: tuple[float, float] = (-2.0, 6.0),
    sigma_bounds: tuple[float, float] = (1e-4, 5.0),
    tries: list[tuple[float, float]] | None = None,
) -> InversionResult:
    """Invert truncated normal statistics to recover underlying parameters.

    Solve for (mu, sigma) such that truncated-mean(mu, sigma) = m_obs and
    truncated-std(mu, sigma) = s_obs, where the distribution is truncated
    from below at lnR_cut.

    :param m_obs: Observed truncated mean
    :param s_obs: Observed truncated standard deviation
    :param lnR_cut: Lower truncation point
    :param mu_bounds: Bounds for mu search (default: (-2.0, 6.0))
    :param sigma_bounds: Bounds for sigma search (default: (1e-4, 5.0))
    :param tries: Initial guesses as list of (mu0, log_sigma0) tuples (default: None)
    :return: Dictionary with keys: mu, sigma, success, message, and optimization info
    """

    # objective on parameters x = [mu, t] with sigma = exp(t) to enforce positivity
    def residuals(x):
        mu = x[0]
        sigma = np.exp(x[1])
        r1 = mean_lnR_truncated(mu, sigma, lnR_cut) - m_obs
        r2 = std_lnR_truncated(mu, sigma, lnR_cut) - s_obs
        return np.array([r1, r2])

    # default initial guesses (several tries for robustness)
    if tries is None:
        tries = [
            (m_obs, s_obs),
        ]

    # bounds in transformed space
    mu_lo, mu_hi = mu_bounds
    sig_lo, sig_hi = sigma_bounds
    # convert sigma bounds to t bounds
    t_lo, t_hi = np.log(sig_lo), np.log(sig_hi)
    lower = np.array([mu_lo, t_lo])
    upper = np.array([mu_hi, t_hi])

    best: OptimizeResult | None = None
    for i, (mu0, logsigma0) in enumerate(tries):
        x0 = np.array([mu0, logsigma0])
        r0 = np.linalg.norm(residuals(x0))
        for my_try, logsigma_try in zip(
            np.random.uniform(mu_lo, mu_hi, 20), np.random.uniform(t_lo, t_hi, 20)
        ):
            r = residuals([my_try, logsigma_try])
            if np.linalg.norm(r) < r0:
                r0 = np.linalg.norm(r)
                x0 = np.array([my_try, logsigma_try])

        try:
            sol = least_squares(
                residuals, x0, xtol=2.23e-16, bounds=(lower, upper), max_nfev=2000
            )
            # pick best by cost
            if best is None:
                best = sol
            elif sol.cost < best.cost:
                best = sol
        except (ValueError, RuntimeError) as e:
            print(
                (
                    f"  -- optimizer failed on try {i} -- exception: {e} "
                    f"-- initial guess: x0={x0} bounds=({lower},{upper})"
                ),
                file=sys.stderr,
                flush=True,
            )

    if best is None:
        return {
            "success": False,
            "message": "optimizer failed",
            "mu": None,
            "sigma": None,
        }

    mu_hat = best.x[0]
    sigma_hat = float(np.exp(best.x[1]))
    return {
        "success": bool(best.success),
        "message": best.message,
        "mu": float(mu_hat),
        "sigma": float(sigma_hat),
        "cost": float(best.cost),
        "nfev": int(best.nfev),
        "optimality": float(best.optimality),
    }


def invert_truncated_stats_mu_from_sample(
    lnR_sample: np.ndarray, lnR_cut: float
) -> float | None:
    """Recover underlying mu from a sample of truncated data.

    :param lnR_sample: Array of log-richness values from truncated distribution
    :param lnR_cut: Lower truncation point
    :return: Estimated underlying mu, or None if inversion failed
    """
    assert len(lnR_sample) > 0, "Sample must have at least one element"
    result = invert_truncated_stats(np.mean(lnR_sample), np.std(lnR_sample), lnR_cut)
    return result["mu"]


def invert_truncated_stats_sigma_from_sample(
    lnR_sample: np.ndarray, lnR_cut: float
) -> float | None:
    """Recover underlying sigma from a sample of truncated data.

    :param lnR_sample: Array of log-richness values from truncated distribution
    :param lnR_cut: Lower truncation point
    :return: Estimated underlying sigma, or None if inversion failed
    """
    assert len(lnR_sample) > 0, "Sample must have at least one element"
    result = invert_truncated_stats(np.mean(lnR_sample), np.std(lnR_sample), lnR_cut)
    return result["sigma"]
