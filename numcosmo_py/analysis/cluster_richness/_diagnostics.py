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
    invert_truncated_stats_mu_from_sample,
    invert_truncated_stats_sigma_from_sample,
)


def compute_binned_statistics(
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    lnR_cut: float,
    nbins: int = 40,
) -> dict[str, np.ndarray]:
    """Compute binned statistics for diagnostic analysis.

    :param mu: Predicted truncated mean values
    :param std_pred: Predicted truncated sigma values
    :param lnR: Observed log-richness values
    :param mu_theo: Theoretical (untruncated) mean values
    :param sigma_theo: Theoretical (untruncated) sigma values
    :param lnR_cut: Richness cut value (log units)
    :param nbins: Number of bins (default: 40)
    :return: Dictionary of computed statistics
    """
    # Compute bins using quantiles for equal-count bins
    mean_bins = np.quantile(mean_pred, np.linspace(0, 1, nbins + 1))
    mean_bin_centers = 0.5 * (mean_bins[:-1] + mean_bins[1:])

    # Compute bins using quantiles for equal-count bins
    mu_bins = np.quantile(mu, np.linspace(0, 1, nbins + 1))
    mu_bin_centers = 0.5 * (mu_bins[:-1] + mu_bins[1:])

    # Compute empirical mean and variance in bins
    sample_count, _, _ = binned_statistic(
        mean_pred, lnR, statistic="count", bins=mean_bins
    )
    sample_mean, _, _ = binned_statistic(
        mean_pred, lnR, statistic="mean", bins=mean_bins
    )
    sample_std, _, _ = binned_statistic(
        mean_pred, lnR, statistic=lambda x: np.std(x, ddof=1), bins=mean_bins
    )
    bin_mean_pred, _, _ = binned_statistic(
        mean_pred, mean_pred, statistic="mean", bins=mean_bins
    )
    bin_std_pred, _, _ = binned_statistic(
        mean_pred, std_pred, statistic="mean", bins=mean_bins
    )

    def _invert_sigma(x: np.ndarray) -> float:
        if len(x) == 0:
            return np.nan
        result = invert_truncated_stats_sigma_from_sample(x, lnR_cut)
        return result if result is not None else np.nan

    def _invert_mu(x: np.ndarray) -> float:
        if len(x) == 0:
            return np.nan
        result = invert_truncated_stats_mu_from_sample(x, lnR_cut)
        return result if result is not None else np.nan

    # Compute predicted mean and variance in bins
    mu_count, _, _ = binned_statistic(mu, mu, statistic="count", bins=mu_bins)
    bin_mu, _, _ = binned_statistic(mu, mu, statistic="mean", bins=mu_bins)
    bin_sigma, _, _ = binned_statistic(mu, sigma, statistic="mean", bins=mu_bins)

    sample_mu, _, _ = binned_statistic(mu, lnR, statistic=_invert_mu, bins=mu_bins)
    sample_sigma, _, _ = binned_statistic(
        mu, lnR, statistic=_invert_sigma, bins=mu_bins
    )

    return {
        # Bins based on truncated mean predictions
        "mean_bins": mean_bins,
        "mean_bin_centers": mean_bin_centers,
        # Empirical statistics in mean_pred bins
        "sample_count": sample_count,
        "sample_mean": sample_mean,
        "sample_std": sample_std,
        # Predicted values in mean_pred bins
        "bin_mean_pred": bin_mean_pred,
        "bin_std_pred": bin_std_pred,
        # Bins based on theoretical (untruncated) mu
        "mu_bins": mu_bins,
        "mu_bin_centers": mu_bin_centers,
        # Theoretical statistics in mu bins
        "mu_count": mu_count,
        "bin_mu": bin_mu,
        "bin_sigma": bin_sigma,
        # Recovered parameters from inversion in mu bins
        "sample_mu": sample_mu,
        "sample_sigma": sample_sigma,
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

    ax.plot(
        stats["mu_bin_centers"],
        stats["sample_mu"],
        "o-",
        markersize=4,
        label=r"Recovered $\mu$ (from inversion)",
    )
    ax.plot(
        stats["mu_bin_centers"],
        stats["bin_mu"],
        "k--",
        lw=1.5,
        label=r"True $\mu$ (untruncated)",
    )
    ax.set_xlabel(r"Theoretical mean $\mu$")
    ax.set_ylabel(r"Recovered mean $\mu$")
    ax.set_title(r"$\mu$ Recovery from Truncated Distribution")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3)
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

    ax.plot(
        stats["mu_bin_centers"],
        stats["sample_sigma"],
        "o-",
        markersize=4,
        label=r"Recovered $\sigma$ (from inversion)",
    )
    ax.plot(
        stats["mu_bin_centers"],
        stats["bin_sigma"],
        "s--",
        markersize=4,
        label=r"True $\sigma$ (untruncated)",
    )
    ax.set_xlabel(r"Theoretical mean $\mu$")
    ax.set_ylabel(r"Recovered dispersion $\sigma$")
    ax.set_title(r"$\sigma$ Recovery from Truncated Distribution")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3)
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

    ax.plot(stats["mean_bin_centers"], stats["sample_count"], "o-", markersize=4)
    ax.set_xlabel(r"Predicted truncated mean $\langle \ln R \rangle$")
    ax.set_ylabel("Clusters per bin")
    ax.set_title("Bin Population")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3)
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
        yerr = stats["sample_std"] / np.sqrt(stats["sample_count"])
    else:
        yerr = stats["sample_std"]

    ax.errorbar(
        stats["mean_bin_centers"],
        stats["sample_mean"],
        yerr=yerr,
        fmt="o",
        markersize=4,
        capsize=2,
        label=r"Sample mean $\langle \ln R \rangle$",
    )
    ax.plot(
        stats["mean_bin_centers"],
        stats["bin_mean_pred"],
        "k--",
        lw=1.5,
        label="Predicted truncated mean",
    )
    ax.set_xlabel(r"Predicted truncated mean $\langle \ln R \rangle$")
    ax.set_ylabel(r"Observed $\ln R$")
    ax.set_title("Mean Richness: Data vs Model")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3)
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
        stats["mean_bin_centers"],
        stats["sample_std"],
        "o",
        markersize=4,
        label=r"Sample std $\sigma(\ln R)$",
    )
    ax.plot(
        stats["mean_bin_centers"],
        stats["bin_std_pred"],
        "-",
        lw=1.5,
        label="Predicted truncated std",
    )
    ax.set_xlabel(r"Predicted truncated mean $\langle \ln R \rangle$")
    ax.set_ylabel(r"Dispersion $\sigma(\ln R)$")
    ax.set_title("Dispersion: Data vs Model")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3)
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

    residuals = stats["sample_std"] - stats["bin_std_pred"]

    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.plot(stats["mean_bin_centers"], residuals, "o-", markersize=4)
    ax.set_xlabel(r"Predicted truncated mean $\langle \ln R \rangle$")
    ax.set_ylabel(r"$\sigma_{\rm sample} - \sigma_{\rm model}$")
    ax.set_title("Dispersion Residuals")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_diagnostic_summary(
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
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
    stats = compute_binned_statistics(
        mean_pred, std_pred, lnR, mu, sigma, lnR_cut, nbins
    )

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))

    plot_mu_recovery(stats, ax=axes[0, 0], show=False)
    plot_sigma_recovery(stats, ax=axes[0, 1], show=False)
    plot_bin_counts(stats, ax=axes[0, 2], show=False)
    plot_mean_lnR(stats, ax=axes[1, 0], show=False)
    plot_empirical_vs_model_sigma(stats, ax=axes[1, 1], show=False)
    plot_sigma_residuals(stats, ax=axes[1, 2], show=False)

    fig.suptitle(
        f"Diagnostic Summary (richness cut $\\lambda \\geq$ {np.exp(lnR_cut):.1f})",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout()

    if show:
        plt.show()

    return fig, axes
