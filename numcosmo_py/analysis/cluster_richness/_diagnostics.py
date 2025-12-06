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
from scipy.stats import norm

from ._truncated_stats import (
    invert_truncated_stats_mu_from_sample,
    invert_truncated_stats_sigma_from_sample,
    truncated_to_normal,
)


def compute_binned_statistics(
    mean_pred: np.ndarray,
    std_pred: np.ndarray,
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    sigma_lnR: np.ndarray,
    lnR_cut: float,
    nbins: int = 40,
    use_quantiles: bool = False,
) -> dict[str, np.ndarray]:
    """Compute binned statistics for diagnostic analysis.

    Bins data in two ways:
    1. By predicted truncated mean - for comparing predicted vs observed lnR statistics
    2. By model parameter mu - for validating parameter recovery via inversion

    The predictions account for total scatter: sigma_total = sqrt(sigma^2 + sigma_lnR^2),
    which combines the intrinsic model scatter with catalog measurement uncertainties.

    :param mean_pred: Predicted truncated mean from model (forward: mu → mean)
    :param std_pred: Predicted truncated std from model (forward: sigma → std)
    :param lnR: Observed log-richness values from sample
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param sigma_lnR: Catalog uncertainty in ln(richness) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param nbins: Number of bins (default: 40)
    :return: Dictionary of computed statistics for both binning schemes
    """
    # Compute total scatter: intrinsic model + catalog uncertainty
    sigma_total = np.sqrt(sigma**2 + sigma_lnR**2)
    # Compute bins using quantiles for equal-count bins
    if use_quantiles:
        mu_bins = np.quantile(mu, np.linspace(0, 1, nbins + 1))
    else:
        mu_bins = np.linspace(np.min(mu), np.max(mu), nbins + 1)

    mu_bin_centers = 0.5 * (mu_bins[:-1] + mu_bins[1:])

    # Binning by predicted truncated mean: Compare model predictions vs observations
    sample_mean, _, _ = binned_statistic(mu, lnR, statistic="mean", bins=mu_bins)
    sample_std, _, _ = binned_statistic(mu, lnR, statistic="std", bins=mu_bins)
    bin_mean_pred, _, _ = binned_statistic(
        mu, mean_pred, statistic="mean", bins=mu_bins
    )
    bin_std_pred, _, _ = binned_statistic(mu, std_pred, statistic="mean", bins=mu_bins)

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
    bin_sigma_lnR, _, _ = binned_statistic(
        mu, sigma_lnR, statistic="mean", bins=mu_bins
    )
    bin_sigma_total, _, _ = binned_statistic(
        mu, sigma_total, statistic="mean", bins=mu_bins
    )

    sample_mu, _, _ = binned_statistic(mu, lnR, statistic=_invert_mu, bins=mu_bins)
    sample_sigma, _, _ = binned_statistic(
        mu, lnR, statistic=_invert_sigma, bins=mu_bins
    )

    return {
        # Binning scheme 1: By predicted truncated mean (forward validation)
        # Observed sample statistics in mean_pred bins
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
        "bin_sigma_lnR": bin_sigma_lnR,
        "bin_sigma_total": bin_sigma_total,
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

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]

    ax.bar(
        stats["mu_bin_centers"],
        stats["mu_count"],
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

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]

    if with_errors:
        yerr = stats["sample_std"] / np.sqrt(stats["mu_count"])
    else:
        yerr = stats["sample_std"]

    # Bar plot for observed sample mean with error bars
    ax.bar(
        stats["mu_bin_centers"],
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
        mu_bins,
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

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]

    # Bar plot for observed sample std
    ax.bar(
        stats["mu_bin_centers"],
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
        mu_bins,
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

    mu_bins = stats["mu_bins"]
    widths = mu_bins[1:] - mu_bins[:-1]
    residuals = stats["sample_std"] - stats["bin_std_pred"]

    ax.axhline(0, color="gray", lw=1, ls="--")
    ax.bar(
        stats["mu_bin_centers"],
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


def plot_scatter_components(
    stats: dict[str, np.ndarray],
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot scatter components: intrinsic model sigma, catalog uncertainty, and total.

    Shows how catalog uncertainties contribute to the total observed scatter.
    The total scatter is what the likelihood uses: sigma_total = sqrt(sigma^2 + sigma_lnR^2).

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

    ax.stairs(
        stats["bin_sigma"],
        mu_bins,
        color="C0",
        lw=2,
        baseline=None,
        label=r"Intrinsic model $\sigma(\ln M, z)$",
    )
    ax.stairs(
        stats["bin_sigma_lnR"],
        mu_bins,
        color="C2",
        lw=2,
        baseline=None,
        label=r"Catalog uncertainty $\sigma_{\ln R}$",
    )
    ax.stairs(
        stats["bin_sigma_total"],
        mu_bins,
        color="black",
        lw=2,
        baseline=None,
        label=r"Total scatter $\sqrt{\sigma^2 + \sigma_{\ln R}^2}$",
    )
    ax.stairs(
        stats["sample_sigma"],
        mu_bins,
        color="C1",
        lw=2,
        linestyle="--",
        baseline=None,
        label=r"Recovered $\sigma$ (inverse: sample → model)",
    )
    ax.set_xlabel(r"Model parameter $\mu(\ln M, z)$")
    ax.set_ylabel(r"Scatter $\sigma$")
    ax.set_title("Scatter Components")
    ax.legend(loc="best", fontsize="small")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_standardized_residuals_qq(
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    lnR_cut: float,
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot Q-Q plot of standardized residuals.

    Transforms truncated normal samples to standard normal using the model
    parameters (mu, sigma) and compares against theoretical N(0,1) quantiles.
    If the model is correct, points should fall on the diagonal.

    :param lnR: Observed log-richness values
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(6, 6))
    else:
        fig = cast(Figure, ax.get_figure())

    # Transform truncated samples to standard normal
    z_obs = truncated_to_normal(lnR, mu, sigma, lnR_cut)

    # Sort observed z-scores
    z_sorted = np.sort(z_obs)

    # Theoretical quantiles from N(0,1)
    n = len(z_sorted)
    p = (np.arange(1, n + 1) - 0.5) / n
    z_theoretical = np.percentile(np.random.standard_normal(100000), p * 100)

    # Q-Q plot
    ax.scatter(z_theoretical, z_sorted, alpha=0.5, s=10, edgecolors="none")

    # Add diagonal reference line
    lim = max(
        abs(z_theoretical.min()),
        abs(z_theoretical.max()),
        abs(z_sorted.min()),
        abs(z_sorted.max()),
    )
    ax.plot([-lim, lim], [-lim, lim], "r--", lw=2, label="Perfect fit")

    ax.set_xlabel("Theoretical N(0,1) Quantiles")
    ax.set_ylabel("Observed Standardized Residuals")
    ax.set_title("Q-Q Plot: Standardized Residuals vs N(0,1)")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    ax.set_aspect("equal")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_standardized_residuals_histogram(
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    lnR_cut: float,
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot histogram of standardized residuals compared to N(0,1).

    Transforms truncated normal samples to standard normal using the model
    parameters and overlays the expected N(0,1) density. If the model is
    correct, the histogram should match the theoretical distribution.

    :param lnR: Observed log-richness values
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = cast(Figure, ax.get_figure())

    # Transform truncated samples to standard normal
    z_obs = truncated_to_normal(lnR, mu, sigma, lnR_cut)

    # Histogram of standardized residuals
    ax.hist(
        z_obs,
        bins=50,
        density=True,
        alpha=0.7,
        edgecolor="black",
        label="Observed residuals",
    )

    assert isinstance(z_obs, np.ndarray)
    z_range = np.linspace(z_obs.min() - 0.5, z_obs.max() + 0.5, 200)
    ax.plot(z_range, norm.pdf(z_range), "r-", lw=2, label="Theoretical N(0,1)")

    ax.set_xlabel("Standardized Residuals")
    ax.set_ylabel("Density")
    ax.set_title("Distribution of Standardized Residuals")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3, axis="y")
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_standardized_residuals_vs_predicted(
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    lnR_cut: float,
    ax: Axes | None = None,
    show: bool = True,
) -> tuple[Figure, Axes]:
    """Plot standardized residuals vs predicted values.

    Helps identify systematic patterns or heteroscedasticity in the residuals.
    Points should be randomly scattered around zero with no trends.

    :param lnR: Observed log-richness values
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param ax: Matplotlib axes to plot on (default: None, creates new figure)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and Axes objects
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 5))
    else:
        fig = cast(Figure, ax.get_figure())

    # Transform truncated samples to standard normal
    z_obs = truncated_to_normal(lnR, mu, sigma, lnR_cut)

    # Plot residuals vs predicted mu
    ax.scatter(mu, z_obs, alpha=0.5, s=10, edgecolors="none")
    ax.axhline(0, color="r", linestyle="--", lw=2, label="Zero residual")

    # Add reference bands at ±2 sigma
    ax.axhline(2, color="gray", linestyle=":", lw=1, alpha=0.5)
    ax.axhline(-2, color="gray", linestyle=":", lw=1, alpha=0.5)

    ax.set_xlabel(r"Predicted $\mu(\ln M, z)$")
    ax.set_ylabel("Standardized Residuals")
    ax.set_title("Residuals vs Predicted Values")
    ax.legend(loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if show:
        plt.show()

    return fig, ax


def plot_residuals_summary(
    lnR: np.ndarray,
    mu: np.ndarray,
    sigma: np.ndarray,
    lnR_cut: float,
    show: bool = True,
) -> tuple[Figure, np.ndarray]:
    """Create a summary figure with standardized residual diagnostics.

    :param lnR: Observed log-richness values
    :param mu: Model parameter mu(lnM, z) for each cluster
    :param sigma: Model parameter sigma(lnM, z) for each cluster
    :param lnR_cut: Richness cut value (log units)
    :param show: Whether to call plt.show() (default: True)
    :return: Figure and array of Axes objects
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    plot_standardized_residuals_qq(lnR, mu, sigma, lnR_cut, ax=axes[0], show=False)
    plot_standardized_residuals_histogram(
        lnR, mu, sigma, lnR_cut, ax=axes[1], show=False
    )
    plot_standardized_residuals_vs_predicted(
        lnR, mu, sigma, lnR_cut, ax=axes[2], show=False
    )

    fig.suptitle(
        f"Standardized Residuals Diagnostics (richness cut $\\lambda \\geq$ {np.exp(lnR_cut):.1f})",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout()

    if show:
        plt.show()

    return fig, axes


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
    fig, axes = plt.subplots(3, 3, figsize=(15, 12))

    plot_mu_recovery(stats, ax=axes[0, 0], show=False)
    plot_sigma_recovery(stats, ax=axes[0, 1], show=False)
    plot_scatter_components(stats, ax=axes[0, 2], show=False)
    plot_bin_counts(stats, ax=axes[1, 0], show=False)
    plot_mean_lnR(stats, ax=axes[1, 1], show=False)
    plot_empirical_vs_model_sigma(stats, ax=axes[1, 2], show=False)
    plot_sigma_residuals(stats, ax=axes[2, 0], show=False)
    # Leave axes[2, 1] and axes[2, 2] empty for future additions
    axes[2, 1].axis("off")
    axes[2, 2].axis("off")

    fig.suptitle(
        f"Diagnostic Summary (richness cut $\\lambda \\geq$ {np.exp(lnR_cut):.1f})",
        fontsize=14,
        fontweight="bold",
    )
    fig.tight_layout()

    if show:
        plt.show()

    return fig, axes
