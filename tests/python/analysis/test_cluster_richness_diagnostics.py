#!/usr/bin/env python
#
# test_cluster_richness_diagnostics.py
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

"""Tests for cluster richness diagnostic functions."""

import numpy as np
from numpy.random import RandomState  # pylint: disable=no-name-in-module
import pytest
import matplotlib
import matplotlib.pyplot as plt
from typing import Any

# Use non-interactive backend for testing
# flake8: noqa: E402
# pylint: disable=wrong-import-position
matplotlib.use("Agg")

from numcosmo_py import Ncm
from numcosmo_py.analysis.cluster_richness import (
    compute_binned_statistics,
    plot_mu_recovery,
    plot_sigma_recovery,
    plot_bin_counts,
    plot_mean_lnR,
    plot_empirical_vs_model_sigma,
    mean_lnR_truncated,
    std_lnR_truncated,
)
from numcosmo_py.analysis.cluster_richness._diagnostics import (
    plot_sigma_residuals,
    plot_scatter_components,
    plot_standardized_residuals_qq,
    plot_standardized_residuals_histogram,
    plot_standardized_residuals_vs_predicted,
    plot_residuals_summary,
    plot_diagnostic_summary,
)

Ncm.cfg_init()


@pytest.fixture(name="mock_diagnostic_data")
def fixture_mock_diagnostic_data() -> dict[str, Any]:
    """Create mock data for diagnostic testing.

    Returns a dictionary with:
    - Model parameters (mu, sigma)
    - Observed data (lnR)
    - Predicted statistics (mean_pred, std_pred)
    - Catalog uncertainties (sigma_lnR)
    - Cut value (lnR_cut)
    """
    rng = RandomState(12345)

    n_clusters = 200

    # Generate model parameters
    lnM = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
    z = rng.uniform(0.1, 0.5, n_clusters)

    # Create a simple model: mu = a + b * lnM, sigma = constant
    a, b = 2.0, 0.8
    mu = a + b * (lnM - np.log(3e14))
    sigma = np.full(n_clusters, 0.3)

    # Generate observed data with truncation
    lnR_cut = np.log(10)  # Cut at richness = 10

    # Sample from truncated normal
    lnR_list = []
    for m, s in zip(mu, sigma):
        while True:
            val = rng.normal(m, s)
            if val >= lnR_cut:
                lnR_list.append(val)
                break

    lnR = np.array(lnR_list)

    # Add catalog uncertainties
    sigma_lnR = rng.uniform(0.05, 0.15, n_clusters)

    # Compute total scatter
    sigma_total = np.sqrt(sigma**2 + sigma_lnR**2)

    # Compute predicted truncated statistics
    mean_pred = mean_lnR_truncated(mu, sigma_total, lnR_cut)
    std_pred = std_lnR_truncated(mu, sigma_total, lnR_cut)

    return {
        "mu": mu,
        "sigma": sigma,
        "lnR": lnR,
        "sigma_lnR": sigma_lnR,
        "sigma_total": sigma_total,
        "mean_pred": mean_pred,
        "std_pred": std_pred,
        "lnR_cut": lnR_cut,
        "lnM": lnM,
        "z": z,
    }


@pytest.fixture(name="diagnostic_stats")
def fixture_diagnostic_stats(mock_diagnostic_data: dict[str, Any]) -> dict[str, Any]:
    """Compute diagnostic statistics from mock data."""
    data = mock_diagnostic_data
    stats = compute_binned_statistics(
        data["mean_pred"],
        data["std_pred"],
        data["lnR"],
        data["mu"],
        data["sigma"],
        data["sigma_lnR"],
        data["lnR_cut"],
        nbins=10,  # Use fewer bins for testing
        use_quantiles=False,
    )
    return stats


class TestComputeBinnedStatistics:
    """Tests for compute_binned_statistics function."""

    def test_basic_computation(self, mock_diagnostic_data: dict[str, Any]) -> None:
        """Test basic computation of binned statistics."""
        data = mock_diagnostic_data
        stats = compute_binned_statistics(
            data["mean_pred"],
            data["std_pred"],
            data["lnR"],
            data["mu"],
            data["sigma"],
            data["sigma_lnR"],
            data["lnR_cut"],
            nbins=10,
        )

        # Check that all expected keys are present
        expected_keys = [
            "sample_mean",
            "sample_std",
            "bin_mean_pred",
            "bin_std_pred",
            "mu_bins",
            "mu_bin_centers",
            "mu_count",
            "bin_mu",
            "bin_sigma",
            "bin_sigma_lnR",
            "bin_sigma_total",
            "sample_mu",
            "sample_sigma",
        ]
        for key in expected_keys:
            assert key in stats, f"Missing key: {key}"

        # Check array shapes
        assert len(stats["mu_bins"]) == 11  # nbins + 1
        assert len(stats["sample_mean"]) == 10
        assert len(stats["mu_bin_centers"]) == 10

    def test_quantile_binning(self, mock_diagnostic_data: dict[str, Any]) -> None:
        """Test binning with quantiles."""
        data = mock_diagnostic_data
        stats = compute_binned_statistics(
            data["mean_pred"],
            data["std_pred"],
            data["lnR"],
            data["mu"],
            data["sigma"],
            data["sigma_lnR"],
            data["lnR_cut"],
            nbins=10,
            use_quantiles=True,
        )

        # Check that bins were computed
        assert "mu_bins" in stats
        assert len(stats["mu_bins"]) == 11

        # Check that bins span the data range
        assert stats["mu_bins"][0] <= data["mu"].min()
        assert stats["mu_bins"][-1] >= data["mu"].max()

    def test_different_bin_counts(self, mock_diagnostic_data: dict[str, Any]) -> None:
        """Test with different numbers of bins."""
        data = mock_diagnostic_data

        for nbins in [5, 10, 20]:
            stats = compute_binned_statistics(
                data["mean_pred"],
                data["std_pred"],
                data["lnR"],
                data["mu"],
                data["sigma"],
                data["sigma_lnR"],
                data["lnR_cut"],
                nbins=nbins,
            )

            assert len(stats["mu_bins"]) == nbins + 1
            assert len(stats["sample_mean"]) == nbins

    def test_statistics_values_reasonable(
        self, mock_diagnostic_data: dict[str, Any]
    ) -> None:
        """Test that computed statistics have reasonable values."""
        data = mock_diagnostic_data
        stats = compute_binned_statistics(
            data["mean_pred"],
            data["std_pred"],
            data["lnR"],
            data["mu"],
            data["sigma"],
            data["sigma_lnR"],
            data["lnR_cut"],
            nbins=10,
        )

        # Sample mean should be above the cut
        valid_means = stats["sample_mean"][~np.isnan(stats["sample_mean"])]
        assert np.all(valid_means >= data["lnR_cut"])

        # Sample std should be positive
        valid_stds = stats["sample_std"][~np.isnan(stats["sample_std"])]
        assert np.all(valid_stds > 0)

        # Counts should be non-negative
        assert np.all(stats["mu_count"] >= 0)


class TestPlotMuRecovery:
    """Tests for plot_mu_recovery function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic mu recovery plot."""
        fig, ax = plot_mu_recovery(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""
        assert ax.get_title() != ""

        plt.close(fig)

    def test_with_existing_axes(self, diagnostic_stats):
        """Test plotting on existing axes."""
        fig, ax = plt.subplots()
        fig_out, ax_out = plot_mu_recovery(diagnostic_stats, ax=ax, show=False)

        assert fig_out is fig
        assert ax_out is ax

        plt.close(fig)

    def test_has_legend(self, diagnostic_stats: dict[str, Any]) -> None:
        """Test that plot has a legend."""
        fig, ax = plot_mu_recovery(diagnostic_stats, show=False)

        legend = ax.get_legend()
        assert legend is not None
        assert len(legend.get_texts()) > 0

        plt.close(fig)


class TestPlotSigmaRecovery:
    """Tests for plot_sigma_recovery function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic sigma recovery plot."""
        fig, ax = plot_sigma_recovery(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""
        assert ax.get_title() != ""

        plt.close(fig)

    def test_with_existing_axes(self, diagnostic_stats):
        """Test plotting on existing axes."""
        fig, ax = plt.subplots()
        fig_out, ax_out = plot_sigma_recovery(diagnostic_stats, ax=ax, show=False)

        assert fig_out is fig
        assert ax_out is ax

        plt.close(fig)


class TestPlotBinCounts:
    """Tests for plot_bin_counts function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic bin counts plot."""
        fig, ax = plot_bin_counts(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""
        assert ax.get_title() != ""

        plt.close(fig)

    def test_log_scale(self, diagnostic_stats: dict[str, Any]) -> None:
        """Test that y-axis uses log scale."""
        fig, ax = plot_bin_counts(diagnostic_stats, show=False)

        assert ax.get_yscale() == "log"

        plt.close(fig)


class TestPlotMeanLnR:
    """Tests for plot_mean_lnR function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic mean lnR plot."""
        fig, ax = plot_mean_lnR(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None

        plt.close(fig)

    def test_with_errors(self, diagnostic_stats: dict[str, Any]) -> None:
        """Test plot with error bars."""
        fig, ax = plot_mean_lnR(diagnostic_stats, with_errors=True, show=False)

        assert fig is not None
        assert ax is not None

        plt.close(fig)

    def test_without_errors(self, diagnostic_stats: dict[str, Any]) -> None:
        """Test plot without error bars."""
        fig, ax = plot_mean_lnR(diagnostic_stats, with_errors=False, show=False)

        assert fig is not None
        assert ax is not None

        plt.close(fig)


class TestPlotEmpiricalVsModelSigma:
    """Tests for plot_empirical_vs_model_sigma function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic empirical vs model sigma plot."""
        fig, ax = plot_empirical_vs_model_sigma(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)

    def test_with_existing_axes(self, diagnostic_stats):
        """Test plotting on existing axes."""
        fig, ax = plt.subplots()
        fig_out, ax_out = plot_empirical_vs_model_sigma(
            diagnostic_stats, ax=ax, show=False
        )

        assert fig_out is fig
        assert ax_out is ax

        plt.close(fig)


class TestPlotSigmaResiduals:
    """Tests for plot_sigma_residuals function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic sigma residuals plot."""
        fig, ax = plot_sigma_residuals(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)


class TestPlotScatterComponents:
    """Tests for plot_scatter_components function."""

    def test_basic_plot(self, diagnostic_stats):
        """Test basic scatter components plot."""
        fig, ax = plot_scatter_components(diagnostic_stats, show=False)

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)

    def test_has_multiple_lines(self, diagnostic_stats: dict[str, Any]) -> None:
        """Test that plot shows multiple scatter components."""
        fig, ax = plot_scatter_components(diagnostic_stats, show=False)

        # Should have multiple scatter components (uses stairs, not lines)
        # Check that legend has multiple entries
        legend = ax.get_legend()
        assert legend is not None
        assert len(legend.get_texts()) >= 3

        plt.close(fig)


class TestPlotStandardizedResidualsQQ:
    """Tests for plot_standardized_residuals_qq function."""

    def test_basic_plot(self, mock_diagnostic_data):
        """Test basic Q-Q plot of standardized residuals."""
        data = mock_diagnostic_data
        fig, ax = plot_standardized_residuals_qq(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)

    def test_with_existing_axes(self, mock_diagnostic_data):
        """Test Q-Q plot on existing axes."""
        data = mock_diagnostic_data
        fig, ax = plt.subplots()
        fig_out, ax_out = plot_standardized_residuals_qq(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            ax=ax,
            show=False,
        )

        assert fig_out is fig
        assert ax_out is ax

        plt.close(fig)


class TestPlotStandardizedResidualsHistogram:
    """Tests for plot_standardized_residuals_histogram function."""

    def test_basic_plot(self, mock_diagnostic_data):
        """Test basic histogram of standardized residuals."""
        data = mock_diagnostic_data
        fig, ax = plot_standardized_residuals_histogram(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)

    def test_plot_has_histogram(self, mock_diagnostic_data: dict[str, Any]) -> None:
        """Test that plot contains a histogram."""
        data = mock_diagnostic_data
        fig, ax = plot_standardized_residuals_histogram(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        assert fig is not None
        assert ax is not None
        # Check that there's a legend with histogram and theoretical distribution
        legend = ax.get_legend()
        assert legend is not None

        plt.close(fig)


class TestPlotStandardizedResidualsVsPredicted:
    """Tests for plot_standardized_residuals_vs_predicted function."""

    def test_basic_plot(self, mock_diagnostic_data):
        """Test basic scatter plot of residuals vs predicted."""
        data = mock_diagnostic_data
        fig, ax = plot_standardized_residuals_vs_predicted(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        assert fig is not None
        assert ax is not None
        assert ax.get_xlabel() != ""
        assert ax.get_ylabel() != ""

        plt.close(fig)

    def test_with_existing_axes(self, mock_diagnostic_data):
        """Test scatter plot on existing axes."""
        data = mock_diagnostic_data
        fig, ax = plt.subplots()
        fig_out, ax_out = plot_standardized_residuals_vs_predicted(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            ax=ax,
            show=False,
        )

        assert fig_out is fig
        assert ax_out is ax

        plt.close(fig)


class TestPlotResidualsSummary:
    """Tests for plot_residuals_summary function."""

    def test_basic_plot(self, mock_diagnostic_data):
        """Test basic residuals summary with 3 subplots."""
        data = mock_diagnostic_data
        fig, axes = plot_residuals_summary(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        assert fig is not None
        assert len(axes) == 3
        assert all(ax is not None for ax in axes)

        plt.close(fig)

    def test_all_subplots_have_content(
        self, mock_diagnostic_data: dict[str, Any]
    ) -> None:
        """Test that all subplots have content."""
        data = mock_diagnostic_data
        fig, axes = plot_residuals_summary(
            data["lnR"],
            data["mu"],
            data["sigma_total"],
            data["lnR_cut"],
            show=False,
        )

        # Check each subplot has content
        for ax in axes:
            assert ax.get_xlabel() != ""
            assert ax.get_ylabel() != ""

        plt.close(fig)


class TestPlotDiagnosticSummary:
    """Tests for plot_diagnostic_summary function."""

    def test_basic_plot(
        self, diagnostic_stats: dict[str, Any], mock_diagnostic_data: dict[str, Any]
    ) -> None:
        """Test basic diagnostic summary with 9 subplots."""
        data = mock_diagnostic_data
        fig, axes = plot_diagnostic_summary(
            diagnostic_stats, data["lnR_cut"], show=False
        )

        assert fig is not None
        # Should have 3x3 grid = 9 subplots
        assert len(axes) == 3
        assert len(axes[0]) == 3
        assert len(axes[1]) == 3
        assert len(axes[2]) == 3

        plt.close(fig)

    def test_all_subplots_populated(
        self, diagnostic_stats: dict[str, Any], mock_diagnostic_data: dict[str, Any]
    ) -> None:
        """Test that all subplots are populated."""
        data = mock_diagnostic_data
        fig, axes = plot_diagnostic_summary(
            diagnostic_stats, data["lnR_cut"], show=False
        )

        # Check each subplot has content
        for row in axes:
            for ax in row:
                # Each axes should have some artists
                assert len(ax.get_children()) > 0

        plt.close(fig)


class TestIntegration:
    """Integration tests for diagnostic workflow."""

    def test_full_diagnostic_pipeline(
        self, mock_diagnostic_data: dict[str, Any]
    ) -> None:
        """Test complete diagnostic pipeline."""
        data = mock_diagnostic_data

        # Step 1: Compute binned statistics
        stats = compute_binned_statistics(
            data["mean_pred"],
            data["std_pred"],
            data["lnR"],
            data["mu"],
            data["sigma"],
            data["sigma_lnR"],
            data["lnR_cut"],
            nbins=10,
        )

        # Step 2: Create individual diagnostic plots
        fig1, _ = plot_mu_recovery(stats, show=False)
        fig2, _ = plot_sigma_recovery(stats, show=False)
        fig3, _ = plot_mean_lnR(stats, show=False)
        fig4, _ = plot_empirical_vs_model_sigma(stats, show=False)

        # All plots should be created successfully
        assert all(f is not None for f in [fig1, fig2, fig3, fig4])

        # Clean up
        for fig in [fig1, fig2, fig3, fig4]:
            plt.close(fig)

    def test_summary_plots(self, mock_diagnostic_data: dict[str, Any]) -> None:
        """Test summary plot functions."""
        data = mock_diagnostic_data

        # Compute statistics
        stats = compute_binned_statistics(
            data["mean_pred"],
            data["std_pred"],
            data["lnR"],
            data["mu"],
            data["sigma"],
            data["sigma_lnR"],
            data["lnR_cut"],
            nbins=10,
        )

        # Create summary plots
        fig1, _ = plot_diagnostic_summary(stats, data["lnR_cut"], show=False)
        fig2, _ = plot_residuals_summary(
            data["lnR"], data["mu"], data["sigma_total"], data["lnR_cut"], show=False
        )

        assert fig1 is not None
        assert fig2 is not None

        plt.close(fig1)
        plt.close(fig2)


class TestEdgeCases:
    """Test edge cases and error handling."""

    def test_small_sample(self) -> None:
        """Test with very small sample size."""
        rng = RandomState(999)

        n = 10
        mu = rng.normal(3.0, 0.5, n)
        sigma = np.full(n, 0.3)
        lnR_cut = np.log(5)

        # Sample truncated data
        lnR_list = []
        for m, s in zip(mu, sigma):
            while True:
                val = rng.normal(m, s)
                if val >= lnR_cut:
                    lnR_list.append(val)
                    break
        lnR = np.array(lnR_list)

        sigma_lnR = np.full(n, 0.1)
        sigma_total = np.sqrt(sigma**2 + sigma_lnR**2)
        mean_pred = mean_lnR_truncated(mu, sigma_total, lnR_cut)
        std_pred = std_lnR_truncated(mu, sigma_total, lnR_cut)

        assert isinstance(mean_pred, np.ndarray)
        assert isinstance(std_pred, np.ndarray)
        # Should still work with small sample
        stats = compute_binned_statistics(
            mean_pred,
            std_pred,
            lnR,
            mu,
            sigma,
            sigma_lnR,
            lnR_cut,
            nbins=3,  # Few bins for small sample
        )

        assert "sample_mean" in stats

    def test_with_nans_in_data(self) -> None:
        """Test handling of NaN values in statistics."""
        rng = RandomState(888)

        n = 50
        mu = rng.normal(3.0, 0.5, n)
        sigma = np.full(n, 0.3)
        lnR_cut = np.log(5)

        lnR_list = []
        for m, s in zip(mu, sigma):
            while True:
                val = rng.normal(m, s)
                if val >= lnR_cut:
                    lnR_list.append(val)
                    break
        lnR = np.array(lnR_list)

        sigma_lnR = np.full(n, 0.1)
        sigma_total = np.sqrt(sigma**2 + sigma_lnR**2)
        mean_pred = mean_lnR_truncated(mu, sigma_total, lnR_cut)
        std_pred = std_lnR_truncated(mu, sigma_total, lnR_cut)

        assert isinstance(mean_pred, np.ndarray)
        assert isinstance(std_pred, np.ndarray)

        stats = compute_binned_statistics(
            mean_pred,
            std_pred,
            lnR,
            mu,
            sigma,
            sigma_lnR,
            lnR_cut,
            nbins=10,
        )

        # Some bins may have NaN (empty bins), but function should handle it
        # and not crash
        fig, _ax = plot_mu_recovery(stats, show=False)
        assert fig is not None

        plt.close(fig)
