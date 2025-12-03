#
# _diagnostics.py
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

"""Diagnostic plotting functions for cluster richness analysis.

This module provides visualization tools for analyzing the truncated
mass-richness relation and validating model fits.
"""

from typing import cast

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.axes import Axes
from scipy.stats import binned_statistic

from ._truncated_stats import (
    std_lnR_truncated,
    invert_truncated_stats_mu_from_sample,
    invert_truncated_stats_sigma_from_sample,
)


def compute_binned_statistics(
    mu: np.ndarray,
    lnR: np.ndarray,
    mu_theo: np.ndarray,
    sigma_theo: np.ndarray,
    lnR_cut: float,
    nbins: int = 40,
) -> dict[str, np.ndarray]:
    """Compute binned statistics for diagnostic analysis.

    :param mu: Predicted truncated mean values
    :param lnR: Observed log-richness values
    :param mu_theo: Theoretical (untruncated) mean values
    :param sigma_theo: Theoretical (untruncated) sigma values
    :param lnR_cut: Richness cut value (log units)
    :param nbins: Number of bins (default: 40)
    :return: Dictionary of computed statistics
    """
    # Predicted truncated sigma
    sigma_pred = std_lnR_truncated(mu_theo, sigma_theo, lnR_cut)

    # Compute bins using quantiles for equal-count bins
    mu_bins = np.quantile(mu, np.linspace(0, 1, nbins + 1))
    mu_theo_bins = np.quantile(mu_theo, np.linspace(0, 1, nbins + 1))

    mu_bin_center = 0.5 * (mu_bins[:-1] + mu_bins[1:])
    mu_theo_bin_center = 0.5 * (mu_theo_bins[:-1] + mu_theo_bins[1:])

    # Compute empirical variance in bins
    var_lnR_emp, _, _ = binned_statistic(
        mu, lnR, statistic=lambda x: np.var(x, ddof=1), bins=mu_bins
    )

    # Predicted sigma in same bins - cast to numpy array to satisfy type checker
    sigma_pred_arr = np.asarray(sigma_pred)
    sigma_bin_pred, _, _ = binned_statistic(
        mu, sigma_pred_arr, statistic="mean", bins=mu_bins
    )
    mu_count, _, _ = binned_statistic(mu, mu, statistic="count", bins=mu_bins)
    lnR_mean, _, _ = binned_statistic(mu, lnR, statistic="mean", bins=mu_bins)

    # Recovered parameters via inversion - use wrapper to handle None
    def _invert_sigma(x: np.ndarray) -> float:
        result = invert_truncated_stats_sigma_from_sample(x, lnR_cut)
        return result if result is not None else np.nan

    def _invert_mu(x: np.ndarray) -> float:
        result = invert_truncated_stats_mu_from_sample(x, lnR_cut)
        return result if result is not None else np.nan

    sigma_rec, _, _ = binned_statistic(
        mu_theo,
        lnR,
        statistic=_invert_sigma,
        bins=mu_theo_bins,
    )
    mu_rec, _, _ = binned_statistic(
        mu_theo,
        lnR,
        statistic=_invert_mu,
        bins=mu_theo_bins,
    )
    mu_rec_count, _, _ = binned_statistic(
        mu_theo, mu_theo, statistic="count", bins=mu_theo_bins
    )
    sigma_theo_mean, _, _ = binned_statistic(
        mu_theo, sigma_theo, statistic="mean", bins=mu_theo_bins
    )

    return {
        "mu_bins": mu_bins,
        "mu_bin_center": mu_bin_center,
        "mu_theo_bins": mu_theo_bins,
        "mu_theo_bin_center": mu_theo_bin_center,
        "var_lnR_emp": var_lnR_emp,
        "sigma_bin_pred": sigma_bin_pred,
        "mu_count": mu_count,
        "lnR_mean": lnR_mean,
        "sigma_rec": sigma_rec,
        "mu_rec": mu_rec,
        "mu_rec_count": mu_rec_count,
        "sigma_theo_mean": sigma_theo_mean,
    }


def plot_mu_recovery(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot recovered mu vs theoretical mu.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    ax.plot(stats["mu_theo_bin_center"], stats["mu_rec"], "o-", label="Recovered μ")
    ax.plot(
        stats["mu_theo_bin_center"], stats["mu_theo_bin_center"], "k--", label="True μ"
    )
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\mu$")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_sigma_recovery(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot recovered sigma vs theoretical sigma.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    ax.plot(stats["mu_theo_bin_center"], stats["sigma_rec"], "o-", label="Recovered σ")
    ax.plot(
        stats["mu_theo_bin_center"], stats["sigma_theo_mean"], "s-", label="Predicted σ"
    )
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\sigma$")
    ax.legend()
    ax.grid(True)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_bin_counts(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot number of elements per bin.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    ax.plot(stats["mu_bin_center"], stats["mu_count"], "o-")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel("Number of elements per bin")
    ax.set_yscale("log")
    ax.grid(True)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_mean_lnR(
    stats: dict[str, np.ndarray],
    with_errors: bool = True,
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot mean ln(R) per bin.

    :param stats: Dictionary from compute_binned_statistics
    :param with_errors: Whether to include error bars (default: True)
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    if with_errors:
        yerr = np.sqrt(stats["var_lnR_emp"] / stats["mu_count"])
    else:
        yerr = np.sqrt(stats["var_lnR_emp"])

    ax.errorbar(stats["mu_bin_center"], stats["lnR_mean"], yerr=yerr, fmt="o-")
    ax.plot(stats["mu_bin_center"], stats["mu_bin_center"], "k--")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"$\ln R$")
    ax.grid(True)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_empirical_vs_model_sigma(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot empirical sigma(lnR) vs model prediction.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    ax.plot(
        stats["mu_bin_center"],
        np.sqrt(stats["var_lnR_emp"]),
        "o",
        label=r"Empirical $\sigma(\ln R)$",
    )
    ax.plot(stats["mu_bin_center"], stats["sigma_bin_pred"], "-", label=r"Model σ(μ)")
    ax.set_xlabel(r"$\mu$ (mean ln richness)")
    ax.set_ylabel(r"$\sigma(\ln R)$")
    ax.legend()
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_sigma_residuals(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot residuals between empirical and model sigma.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    residuals = np.sqrt(stats["var_lnR_emp"]) - stats["sigma_bin_pred"]

    ax.axhline(0, color="gray", lw=1)
    ax.plot(stats["mu_bin_center"], residuals, "o-")
    ax.set_xlabel(r"$\mu$")
    ax.set_ylabel(r"Residual: empirical $\sigma$ − model $\sigma$")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_diagnostic_summary(
    mu: np.ndarray,
    lnR: np.ndarray,
    mu_theo: np.ndarray,
    sigma_theo: np.ndarray,
    lnR_cut: float,
    nbins: int = 40,
    show: bool = True,
) -> tuple[Figure, np.ndarray]:
    """Create a summary figure with all diagnostic plots.

    :param mu: Predicted truncated mean values
    :param lnR: Observed log-richness values
    :param mu_theo: Theoretical (untruncated) mean values
    :param sigma_theo: Theoretical (untruncated) sigma values
    :param lnR_cut: Richness cut value (log units)
    :param nbins: Number of bins (default: 40)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and array of Axes objects
    """
    stats = compute_binned_statistics(mu, lnR, mu_theo, sigma_theo, lnR_cut, nbins)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    plot_mu_recovery(stats, ax=axes[0, 0], show=False)
    plot_sigma_recovery(stats, ax=axes[0, 1], show=False)
    plot_bin_counts(stats, ax=axes[0, 2], show=False)
    plot_mean_lnR(stats, ax=axes[1, 0], show=False)
    plot_empirical_vs_model_sigma(stats, ax=axes[1, 1], show=False)
    plot_sigma_residuals(stats, ax=axes[1, 2], show=False)

    fig.suptitle(f"Diagnostic Summary (cut = {np.exp(lnR_cut):.1f})")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, axes
