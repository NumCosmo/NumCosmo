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

This module provides visualization tools for validating mass-richness models
by comparing model predictions with observed data.

Two validation approaches:
1. Forward: Model (mu, sigma) → Truncated (mean, std) → Compare with sample stats
2. Inverse: Sample (mean, std) → Invert to (mu, sigma) → Compare with model

The diagnostics help verify that the model parameters mu(lnM, z) and sigma(lnM, z)
correctly predict the distribution of observed lnR values.
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

    Bins data in two ways:
    1. By predicted truncated mean - for comparing predicted vs observed lnR statistics
    2. By model parameter mu - for validating parameter recovery via inversion

    :param mean_pred: Predicted truncated mean from model (forward: mu → mean)
    :param std_pred: Predicted truncated std from model (forward: sigma → std)
    :param lnR: Observed log-richness values from sample
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param nbins: Number of bins (default: 40)
    :return: Dictionary of computed statistics for both binning schemes
    """
    # Compute bins using quantiles for equal-count bins
    if False:
        mean_bins = np.quantile(mean_pred, np.linspace(0, 1, nbins + 1))
        mu_bins = np.quantile(mu, np.linspace(0, 1, nbins + 1))
    else:
        mean_bins = np.linspace(np.min(mean_pred), np.max(mean_pred), nbins + 1)
        mu_bins = np.linspace(np.min(mu), np.max(mu), nbins + 1)

    mean_bin_centers = 0.5 * (mean_bins[:-1] + mean_bins[1:])
    mu_bin_centers = 0.5 * (mu_bins[:-1] + mu_bins[1:])

    # Binning by predicted truncated mean: Compare model predictions vs observations
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

    # Binning by model parameter mu: Validate parameter recovery via inversion
    mu_count, _, _ = binned_statistic(mu, mu, statistic="count", bins=mu_bins)
    bin_mu, _, _ = binned_statistic(mu, mu, statistic="mean", bins=mu_bins)
    bin_sigma, _, _ = binned_statistic(mu, sigma, statistic="mean", bins=mu_bins)

    sample_mu, _, _ = binned_statistic(mu, lnR, statistic=_invert_mu, bins=mu_bins)
    sample_sigma, _, _ = binned_statistic(
        mu, lnR, statistic=_invert_sigma, bins=mu_bins
    )

    return {
        # Binning scheme 1: By predicted truncated mean (forward validation)
        "mean_bins": mean_bins,
        "mean_bin_centers": mean_bin_centers,
        # Observed sample statistics in mean_pred bins
        "sample_count": sample_count,
        "sample_mean": sample_mean,
        "sample_std": sample_std,
        # Model predictions (forward: mu,sigma → mean,std) in mean_pred bins
        "bin_mean_pred": bin_mean_pred,
        "bin_std_pred": bin_std_pred,
        # Binning scheme 2: By model parameter mu (inverse validation)
        "mu_bins": mu_bins,
        "mu_bin_centers": mu_bin_centers,
        # Model parameters in mu bins
        "mu_count": mu_count,
        "bin_mu": bin_mu,
        "bin_sigma": bin_sigma,
        # Recovered parameters (inverse: mean,std → mu,sigma) in mu bins
        "sample_mu": sample_mu,
        "sample_sigma": sample_sigma,
    }


def plot_mu_recovery(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot recovered mu vs model mu (inverse validation).

    Validates inverse conversion: sample statistics → model parameters.
    Compares mu recovered from sample (mean, std) against model mu(lnM, z).

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]

    ax.bar(
        stats["mu_bin_centers"],
        stats["sample_mu"],
        width=widths,
        edgecolor="C0",
        facecolor="C0",
        alpha=0.7,
        label=r"Recovered $\mu$ (inverse: sample → model)",
    )
    ax.stairs(
        stats["bin_mu"],
        mu_bins,
        color="black",
        lw=2,
        baseline=None,
        label=r"Model $\mu(\ln M, z)$",
    )
    ax.set_xlabel(r"Model parameter $\mu(\ln M, z)$")
    ax.set_ylabel(r"Recovered $\mu$")
    ax.set_title(r"Model Parameter $\mu$ Recovery (Inverse Validation)")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_sigma_recovery(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot recovered sigma vs model sigma (inverse validation).

    Validates inverse conversion: sample statistics → model parameters.
    Compares sigma recovered from sample (mean, std) against model sigma(lnM, z).

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]

    ax.bar(
        stats["mu_bin_centers"],
        stats["sample_sigma"],
        width=widths,
        edgecolor="C0",
        facecolor="C0",
        alpha=0.7,
        label=r"Recovered $\sigma$ (inverse: sample → model)",
    )
    ax.stairs(
        stats["bin_sigma"],
        mu_bins,
        color="black",
        lw=2,
        baseline=None,
        label=r"Model $\sigma(\ln M, z)$",
    )
    ax.set_xlabel(r"Model parameter $\mu(\ln M, z)$")
    ax.set_ylabel(r"Recovered $\sigma$")
    ax.set_title(r"Model Parameter $\sigma$ Recovery (Inverse Validation)")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_bin_counts(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot cluster count per bin.

    Shows the distribution of clusters across bins of predicted truncated mean.
    Helps assess statistical power in different regions.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    mean_bins = stats["mean_bins"]
    widths = mean_bins[1:] - mean_bins[:-1]

    ax.bar(
        stats["mean_bin_centers"],
        stats["sample_count"],
        width=widths,
        edgecolor="black",
        alpha=0.7,
        label="Clusters per bin",
    )
    ax.set_xlabel(r"Predicted truncated mean $\mathrm{E}[\ln R]$ from model")
    ax.set_ylabel("Clusters per bin")
    ax.set_title("Cluster Count Distribution")
    ax.set_yscale("log")
    ax.grid(True, alpha=0.3, axis="y")
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
    """Plot observed vs predicted truncated mean (forward validation).

    Validates forward conversion: model parameters → truncated statistics.
    Compares predicted truncated mean from model against observed sample mean.

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

    mean_bins = stats["mean_bins"]
    widths = mean_bins[1:] - mean_bins[:-1]

    if with_errors:
        yerr = stats["sample_std"] / np.sqrt(stats["sample_count"])
    else:
        yerr = stats["sample_std"]

    # Bar plot for observed sample mean with error bars
    ax.bar(
        stats["mean_bin_centers"],
        stats["sample_mean"],
        width=widths,
        yerr=yerr,
        edgecolor="C0",
        facecolor="C0",
        alpha=0.7,
        capsize=2,
        label=r"Observed sample mean of $\ln R$",
    )

    # Step plot for model prediction (forward: mu,sigma → mean)
    ax.stairs(
        stats["bin_mean_pred"],
        mean_bins,
        color="black",
        lw=2,
        baseline=None,
        label=r"Predicted $\mathrm{E}[\ln R]$ (forward: model → truncated)",
    )
    ax.set_xlabel(r"Predicted truncated mean $\mathrm{E}[\ln R]$ from model")
    ax.set_ylabel(r"Mean $\ln R$")
    ax.set_title("Truncated Mean Comparison (Forward Validation)")
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
    """Plot observed vs predicted truncated std (forward validation).

    Validates forward conversion: model parameters → truncated statistics.
    Compares predicted truncated std from model against observed sample std.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    mean_bins = stats["mean_bins"]
    widths = mean_bins[1:] - mean_bins[:-1]

    # Bar plot for observed sample std
    ax.bar(
        stats["mean_bin_centers"],
        stats["sample_std"],
        width=widths,
        edgecolor="C0",
        facecolor="C0",
        alpha=0.7,
        label=r"Observed sample std of $\ln R$",
    )

    # Step plot for model prediction (forward: mu,sigma → std)
    ax.stairs(
        stats["bin_std_pred"],
        mean_bins,
        color="C1",
        lw=2,
        baseline=None,
        label=r"Predicted $\mathrm{Std}[\ln R]$ (forward: model → truncated)",
    )
    ax.set_xlabel(r"Predicted truncated mean $\mathrm{E}[\ln R]$ from model")
    ax.set_ylabel(r"Std $\ln R$")
    ax.set_title("Truncated Std Comparison (Forward Validation)")
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
    """Plot residuals between observed and predicted truncated std.

    Shows deviations between model predictions (forward) and observed data.
    Helps identify systematic biases in the model.

    :param stats: Dictionary from compute_binned_statistics
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 4))
    else:
        fig = cast(Figure, ax.get_figure())

    mean_bins = stats["mean_bins"]
    widths = mean_bins[1:] - mean_bins[:-1]
    residuals = stats["sample_std"] - stats["bin_std_pred"]

    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.bar(
        stats["mean_bin_centers"],
        residuals,
        width=widths,
        edgecolor="black",
        alpha=0.7,
        color="C2",
    )
    ax.set_xlabel(r"Predicted truncated mean $\mathrm{E}[\ln R]$ from model")
    ax.set_ylabel(r"$\mathrm{Std}_{\rm observed} - \mathrm{Std}_{\rm predicted}$")
    ax.set_title("Truncated Std Residuals")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_diagnostic_summary(
    stats: dict[str, np.ndarray],
    lnR_cut: float,
    show: bool = True,
) -> tuple[Figure, np.ndarray]:
    """Create a summary figure with all diagnostic plots.

    :param stats: Dictionary from compute_binned_statistics
    :param lnR_cut: Richness cut value (log units) for title
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and array of Axes objects
    """
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
