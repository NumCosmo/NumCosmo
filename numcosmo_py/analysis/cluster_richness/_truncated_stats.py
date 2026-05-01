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

This module provides functions for working with truncated normal distributions
in the context of mass-richness relations.

Model: lnR ~ N(mu(lnM, z), sigma(lnM, z)^2) truncated at lnR >= lnR_cut

Independent variables: lnM (log mass), z (redshift)
Model parameters: mu(lnM, z), sigma(lnM, z) - Gaussian distribution parameters
Dependent variable: lnR (log richness)

Functions provided:
- Forward conversion: (mu, sigma) → (truncated mean, truncated std)
- Inverse conversion: (truncated mean, truncated std) → (mu, sigma)
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
    """Calculate the truncated mean from model parameters.

    Forward conversion: Given model parameters mu(lnM, z) and sigma(lnM, z),
    compute the expected value (mean) of the truncated distribution for lnR.

    For lnR ~ N(mu, sigma^2) truncated at lnR >= lnR_cut:
    E[lnR | lnR >= lnR_cut] = mu + sigma * lambda(alpha)

    where alpha = (lnR_cut - mu) / sigma and lambda is the inverse Mills ratio.

    :param mu: Model parameter mu(lnM, z) - location of the underlying Gaussian
    :param sigma: Model parameter sigma(lnM, z) - scale of the underlying Gaussian
    :param lnR_cut: Lower truncation point in log richness
    :return: Truncated mean of lnR
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
    """Calculate the truncated standard deviation from model parameters.

    Forward conversion: Given model parameters mu(lnM, z) and sigma(lnM, z),
    compute the standard deviation of the truncated distribution for lnR.

    For lnR ~ N(mu, sigma^2) truncated at lnR >= lnR_cut:
    Var[lnR | lnR >= lnR_cut] = sigma^2 * [1 + alpha*lambda(alpha) - lambda(alpha)^2]

    where alpha = (lnR_cut - mu) / sigma and lambda is the inverse Mills ratio.

    :param mu: Model parameter mu(lnM, z) - location of the underlying Gaussian
    :param sigma: Model parameter sigma(lnM, z) - scale of the underlying Gaussian
    :param lnR_cut: Lower truncation point in log richness
    :return: Truncated standard deviation of lnR
    """
    A = (lnR_cut - mu) / sigma
    phi = norm.pdf(A)
    S = norm.sf(A)
    T = phi / np.maximum(S, 1e-300)
    corr = np.sqrt(np.maximum(1.0 + A * T - T * T, 0.0))
    return sigma * corr


def invert_truncated_stats(
    mean_obs: float,
    std_obs: float,
    lnR_cut: float,
    mu_bounds: tuple[float, float] = (-2.0, 6.0),
    sigma_bounds: tuple[float, float] = (1e-4, 5.0),
    tries: list[tuple[float, float]] | None = None,
) -> InversionResult:
    """Recover model parameters from observed truncated statistics.

    Inverse conversion: Given observed sample statistics (mean, std) of lnR,
    solve for the underlying model parameters (mu, sigma) that would produce
    these truncated statistics.

    Solves the system:
      mean_lnR_truncated(mu, sigma, lnR_cut) = mean_obs
      std_lnR_truncated(mu, sigma, lnR_cut) = std_obs

    :param mean_obs: Observed sample mean of lnR (truncated)
    :param std_obs: Observed sample standard deviation of lnR (truncated)
    :param lnR_cut: Lower truncation point in log richness
    :param mu_bounds: Search bounds for mu parameter (default: (-2.0, 6.0))
    :param sigma_bounds: Search bounds for sigma parameter (default: (1e-4, 5.0))
    :param tries: Initial guesses as list of (mu0, log_sigma0) tuples (default: None)
    :return: Dictionary with recovered mu, sigma, success flag, and optimization info
    """

    # objective on parameters x = [mu, t] with sigma = exp(t) to enforce positivity
    def residuals(x):
        mu = x[0]
        sigma = np.exp(x[1])
        r1 = mean_lnR_truncated(mu, sigma, lnR_cut) - mean_obs
        r2 = std_lnR_truncated(mu, sigma, lnR_cut) - std_obs
        return np.array([r1, r2])

    # default initial guesses (several tries for robustness)
    if tries is None:
        tries = [(1.0, 0.5)]

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
            np.random.uniform(0.0, mu_hi, 200), np.random.uniform(t_lo, t_hi, 200)
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
    """Recover model parameter mu from observed lnR sample.

    Convenience function that computes sample statistics and inverts them
    to recover the underlying mu parameter.

    :param lnR_sample: Array of observed log-richness values (truncated sample)
    :param lnR_cut: Lower truncation point in log richness
    :return: Estimated model parameter mu, or None if inversion failed
    """
    assert len(lnR_sample) > 0, "Sample must have at least one element"
    result = invert_truncated_stats(np.mean(lnR_sample), np.std(lnR_sample), lnR_cut)
    return result["mu"]


def invert_truncated_stats_sigma_from_sample(
    lnR_sample: np.ndarray, lnR_cut: float
) -> float | None:
    """Recover model parameter sigma from observed lnR sample.

    Convenience function that computes sample statistics and inverts them
    to recover the underlying sigma parameter.

    :param lnR_sample: Array of observed log-richness values (truncated sample)
    :param lnR_cut: Lower truncation point in log richness
    :return: Estimated model parameter sigma, or None if inversion failed
    """
    assert len(lnR_sample) > 0, "Sample must have at least one element"
    result = invert_truncated_stats(np.mean(lnR_sample), np.std(lnR_sample), lnR_cut)
    return result["sigma"]


def truncated_to_normal(
    x: float | np.ndarray,
    mu: float | np.ndarray,
    sigma: float | np.ndarray,
    a: float,
) -> float | np.ndarray:
    """Map samples from one-sided truncated to standard normal.

    It maps truncated N(mu, sigma^2) on [a, ∞) to standard normal N(0,1).

    Numerically stable transformation using the CDF method for one-sided truncation.
    Uses different strategies based on whether the truncation point is in the left
    or right tail to avoid catastrophic cancellation.

    For x ~ TruncatedNormal(mu, sigma, [a, ∞)), returns z ~ Normal(0, 1) such that
    the quantile positions are preserved.

    :param x: Sample(s) from truncated normal distribution (scalar or array)
    :param mu: Location parameter of underlying Gaussian (scalar or array matching x)
    :param sigma: Scale parameter of underlying Gaussian (scalar or array matching x)
    :param a: Lower truncation point (scalar)
    :return: Corresponding standard normal sample(s)
    """
    # Standardize the limits and data
    alpha = (a - mu) / sigma
    t = (x - mu) / sigma

    # Handle both scalar and array cases for alpha
    # Use vectorized operations that work for both scalars and arrays
    SF_alpha = norm.sf(alpha)
    SF_t = norm.sf(t)
    Phi_alpha = norm.cdf(alpha)
    Phi_t = norm.cdf(t)

    # For alpha > 0 (right tail): u = 1 - SF(t)/SF(alpha)
    # For alpha <= 0 (left tail): u = [Φ(t) - Φ(alpha)] / [1 - Φ(alpha)]
    u_right = 1.0 - SF_t / SF_alpha
    u_left = (Phi_t - Phi_alpha) / (1.0 - Phi_alpha)

    # Use np.where to select the appropriate formula based on alpha
    u = np.where(alpha > 0, u_right, u_left)

    # Map uniform [0,1] to standard normal
    # Clamp u to avoid numerical issues at boundaries
    # u = np.clip(u, 1e-15, 1.0 - 1e-15)
    z = norm.ppf(u)

    return z
