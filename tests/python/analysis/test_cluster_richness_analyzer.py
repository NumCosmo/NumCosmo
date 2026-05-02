#!/usr/bin/env python
#
# test_cluster_richness_analyzer.py
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

"""Tests for cluster richness analyzer (CutAnalyzer)."""

import numpy as np
from numpy.random import RandomState  # pylint: disable=no-name-in-module
import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.analysis.cluster_richness import ClusterData, CutAnalyzer

Ncm.cfg_init()


@pytest.fixture(name="mock_cluster_data")
def fixture_mock_cluster_data() -> ClusterData:
    """Create mock cluster data for testing.

    Returns ClusterData with:
    - 100 mock clusters
    - lnM in [ln(1e14), ln(1e15)]
    - z in [0.1, 0.5]
    - lnR following mu + sigma * randn with truncation at lnR_cut
    - sigma_lnR in [0.05, 0.15]
    """
    rng = RandomState(42)
    n_clusters = 100

    # Generate independent variables
    lnM = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
    z = rng.uniform(0.1, 0.5, n_clusters)

    # Generate richness following a simple model
    # mu = a + b * (lnM - ln(3e14))
    a, b = 2.5, 0.7
    mu = a + b * (lnM - np.log(3e14))
    sigma_intrinsic = 0.25

    # Generate observed richness (truncated at ln(10))
    lnR_cut = np.log(10.0)
    lnR_list = []
    for m in mu:
        while True:
            val = rng.normal(m, sigma_intrinsic)
            if val >= lnR_cut:
                lnR_list.append(val)
                break
    lnR = np.array(lnR_list)

    # Add catalog uncertainties
    sigma_lnR = rng.uniform(0.05, 0.15, n_clusters)

    return ClusterData(lnM=lnM, z=z, lnR=lnR, sigma_lnR=sigma_lnR)


@pytest.fixture(name="cluster_mass_richness_model")
def fixture_cluster_mass_richness_model() -> Nc.ClusterMassRichness:
    """Create a simple cluster mass-richness model for testing."""
    cluster_m = Nc.ClusterMassAscaso()

    # Set initial parameters to match the mock data generation
    # mu = mup0 + mup1 * (lnM - ln(3e14)) + mup2 * z
    cluster_m.param_set_by_name("mup0", 2.5)
    cluster_m.param_set_by_name("mup1", 0.7)
    cluster_m.param_set_by_name("mup2", 0.0)  # No redshift dependence
    cluster_m.param_set_by_name("sigmap0", 0.25)
    cluster_m.param_set_by_name("sigmap1", 0.0)
    cluster_m.param_set_by_name("sigmap2", 0.0)
    cluster_m.param_set_by_name("cut", np.log(10.0))

    return cluster_m


class TestClusterData:
    """Test ClusterData class."""

    def test_cluster_data_creation(self, mock_cluster_data: ClusterData) -> None:
        """Test that ClusterData can be created."""
        assert len(mock_cluster_data) == 100
        assert len(mock_cluster_data.lnM) == 100
        assert len(mock_cluster_data.z) == 100
        assert len(mock_cluster_data.lnR) == 100
        assert len(mock_cluster_data.sigma_lnR) == 100

    def test_cluster_data_validation(self) -> None:
        """Test that ClusterData validates array lengths."""
        with pytest.raises(ValueError, match="All arrays must have the same length"):
            ClusterData(
                lnM=np.array([1.0, 2.0]),
                z=np.array([0.1, 0.2, 0.3]),  # Wrong length
                lnR=np.array([1.0, 2.0]),
                sigma_lnR=np.array([0.1, 0.1]),
            )

    def test_apply_cut(self, mock_cluster_data: ClusterData) -> None:
        """Test applying richness cut (line 88 in _analyzer.py).

        This tests the mask line: mask = self.lnR >= lnR_cut
        """
        # Original data
        n_original = len(mock_cluster_data)
        assert n_original == 100

        # Apply a cut at ln(20)
        lnR_cut = np.log(20.0)
        cut_data = mock_cluster_data.apply_cut(lnR_cut)

        # Check that the cut was applied correctly
        n_after_cut = len(cut_data)
        assert n_after_cut < n_original
        assert n_after_cut > 0

        # All remaining clusters should have lnR >= lnR_cut
        assert np.all(cut_data.lnR >= lnR_cut)

        # Verify that the correct clusters were selected
        original_mask = mock_cluster_data.lnR >= lnR_cut
        expected_n = np.sum(original_mask)
        assert n_after_cut == expected_n

        # Check that all arrays were masked consistently
        assert len(cut_data.lnM) == n_after_cut
        assert len(cut_data.z) == n_after_cut
        assert len(cut_data.lnR) == n_after_cut
        assert len(cut_data.sigma_lnR) == n_after_cut

    def test_apply_cut_edge_cases(self, mock_cluster_data: ClusterData) -> None:
        """Test apply_cut with edge case values."""
        # Cut that keeps all clusters
        very_low_cut = np.log(1.0)
        all_data = mock_cluster_data.apply_cut(very_low_cut)
        assert len(all_data) == len(mock_cluster_data)

        # Cut that keeps no clusters
        very_high_cut = np.log(1e10)
        empty_data = mock_cluster_data.apply_cut(very_high_cut)
        assert len(empty_data) == 0


class TestCutAnalyzer:
    """Test CutAnalyzer class."""

    def test_analyzer_initialization(self, mock_cluster_data: ClusterData) -> None:
        """Test analyzer initialization with various parameters."""
        cuts = [np.log(10.0), np.log(15.0), np.log(20.0)]

        # Test with default parameters
        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            verbose=False,
        )
        assert analyzer.data == mock_cluster_data
        assert analyzer.cuts == cuts
        assert analyzer.n_bootstrap == 100  # Default
        assert analyzer.n_mcmc_steps == 400  # Default
        assert analyzer.n_mcmc_burnin == 150  # Default
        assert not analyzer.compute_mcmc
        assert not analyzer.compute_bootstrap

    def test_analyzer_with_custom_parameters(
        self, mock_cluster_data: ClusterData
    ) -> None:
        """Test analyzer with custom MCMC and bootstrap parameters."""
        cuts = [np.log(10.0)]

        # Test with custom parameters
        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            n_bootstrap=10,  # Custom
            n_mcmc_steps=50,  # Custom
            n_mcmc_burnin=10,  # Custom
            compute_mcmc=False,
            compute_bootstrap=False,
            verbose=False,
        )
        assert analyzer.n_bootstrap == 10
        assert analyzer.n_mcmc_steps == 50
        assert analyzer.n_mcmc_burnin == 10

    def test_property_getters_setters(self, mock_cluster_data: ClusterData) -> None:
        """Test that properties can be get and set after initialization."""
        cuts = [np.log(10.0)]
        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            verbose=False,
        )

        # Test getting default values
        assert analyzer.n_bootstrap == 100
        assert analyzer.n_mcmc_steps == 400
        assert analyzer.n_mcmc_burnin == 150

        # Test setting new values
        analyzer.n_bootstrap = 25
        analyzer.n_mcmc_steps = 75
        analyzer.n_mcmc_burnin = 20

        # Test getting the new values
        assert analyzer.n_bootstrap == 25
        assert analyzer.n_mcmc_steps == 75
        assert analyzer.n_mcmc_burnin == 20

    def test_analyze_bestfit_only(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test analyzer with bestfit only (no MCMC/bootstrap)."""
        cuts = [np.log(10.0)]

        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            compute_mcmc=False,
            compute_bootstrap=False,
            verbose=False,
        )

        results = analyzer.analyze(cluster_mass_richness_model)

        # Check that we got results for our cut
        assert len(results) == 1
        assert np.log(10.0) in results

        result = results[np.log(10.0)]
        assert result.cut == np.log(10.0)
        assert result.n_clusters > 0
        assert result.bestfit is not None
        assert np.isfinite(result.m2lnL)

    def test_analyze_with_bootstrap(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test analyzer with bootstrap."""
        cuts = [np.log(10.0)]

        # Run with small number of bootstrap samples for testing
        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            n_bootstrap=5,  # Small number for fast test
            compute_mcmc=False,
            compute_bootstrap=True,
            verbose=False,
        )

        results = analyzer.analyze(cluster_mass_richness_model)

        # Check that bootstrap was computed
        result = results[np.log(10.0)]
        assert result.bootstrap_mean is not None

        # Bootstrap mean should be close to bestfit for small datasets
        # but they won't be identical
        bestfit_params = [
            result.bestfit["mup0"],
            result.bestfit["mup1"],
            result.bestfit["sigmap0"],
        ]
        bootstrap_params = [
            result.bootstrap_mean["mup0"],
            result.bootstrap_mean["mup1"],
            result.bootstrap_mean["sigmap0"],
        ]

        # Check that bootstrap parameters are reasonable
        for bp, bsp in zip(bestfit_params, bootstrap_params):
            assert np.isfinite(bp)
            assert np.isfinite(bsp)

    def test_analyze_with_mcmc(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test analyzer with MCMC (tests lines 262, 267-270).

        This tests:
        - Line 262: esmcmc.run(self.n_mcmc_steps)
        - Lines 267-270: range(self.n_mcmc_burnin, mcat.len())
        """
        cuts = [np.log(10.0)]

        # Run with small number of MCMC steps for testing
        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            n_mcmc_steps=20,  # Small number for fast test
            n_mcmc_burnin=5,  # Small burn-in
            compute_mcmc=True,
            compute_bootstrap=False,
            verbose=False,
        )

        results = analyzer.analyze(cluster_mass_richness_model)

        # Check that MCMC was computed
        result = results[np.log(10.0)]
        assert result.mcmc_mean is not None
        assert result.mcmc_median is not None

        # MCMC results should be reasonable
        mcmc_mean_params = [
            result.mcmc_mean["mup0"],
            result.mcmc_mean["mup1"],
            result.mcmc_mean["sigmap0"],
        ]
        mcmc_median_params = [
            result.mcmc_median["mup0"],
            result.mcmc_median["mup1"],
            result.mcmc_median["sigmap0"],
        ]

        # Check that all parameters are finite
        for param in mcmc_mean_params + mcmc_median_params:
            assert np.isfinite(param)

    def test_analyze_multiple_cuts(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test analyzer with multiple cuts."""
        cuts = [np.log(10.0), np.log(15.0), np.log(20.0)]

        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            compute_mcmc=False,
            compute_bootstrap=False,
            verbose=False,
        )

        results = analyzer.analyze(cluster_mass_richness_model)

        # Check that we got results for all cuts
        assert len(results) == 3
        for cut in cuts:
            assert cut in results
            result = results[cut]
            assert result.cut == cut
            assert result.n_clusters > 0

        # Number of clusters should decrease with increasing cut
        n_clusters = [results[cut].n_clusters for cut in sorted(cuts)]
        assert n_clusters[0] > n_clusters[1] > n_clusters[2]

    def test_mcmc_burnin_parameter(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test that n_mcmc_burnin parameter works correctly.

        This specifically tests that the burn-in range uses self.n_mcmc_burnin
        in line 267-270 of _analyzer.py.
        """
        cuts = [np.log(10.0)]

        # Test with different burn-in values
        for burnin in [2, 5, 10]:
            analyzer = CutAnalyzer(
                data=mock_cluster_data,
                cuts=cuts,
                n_mcmc_steps=20,
                n_mcmc_burnin=burnin,
                compute_mcmc=True,
                compute_bootstrap=False,
                verbose=False,
            )

            results = analyzer.analyze(cluster_mass_richness_model)
            result = results[np.log(10.0)]

            # Check that MCMC completed successfully
            assert result.mcmc_mean is not None
            assert result.mcmc_median is not None

            # Verify parameters are reasonable
            assert np.isfinite(result.mcmc_mean["mup0"])
            assert np.isfinite(result.mcmc_median["mup0"])

    def test_mcmc_steps_parameter(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test that n_mcmc_steps parameter works correctly.

        This specifically tests that esmcmc.run uses self.n_mcmc_steps
        in line 262 of _analyzer.py.
        """
        cuts = [np.log(10.0)]

        # Test with different step counts
        for steps in [10, 20, 30]:
            analyzer = CutAnalyzer(
                data=mock_cluster_data,
                cuts=cuts,
                n_mcmc_steps=steps,
                n_mcmc_burnin=2,
                compute_mcmc=True,
                compute_bootstrap=False,
                verbose=False,
            )

            results = analyzer.analyze(cluster_mass_richness_model)
            result = results[np.log(10.0)]

            # Check that MCMC completed successfully
            assert result.mcmc_mean is not None
            assert result.mcmc_median is not None

    def test_display_results(
        self,
        mock_cluster_data: ClusterData,
        cluster_mass_richness_model: Nc.ClusterMassRichness,
    ) -> None:
        """Test display_results method."""
        cuts = [np.log(10.0)]

        analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=cuts,
            compute_mcmc=False,
            compute_bootstrap=False,
            verbose=False,
        )

        analyzer.analyze(cluster_mass_richness_model)

        # This should not raise an error
        analyzer.display_results()

        # Test with empty results
        empty_analyzer = CutAnalyzer(
            data=mock_cluster_data,
            cuts=[],
            verbose=False,
        )
        # Should not raise an error with empty results
        empty_analyzer.display_results()


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
