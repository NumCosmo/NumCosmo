#!/usr/bin/env python
#
# test_cluster_richness_app.py
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

"""Tests for cluster richness CLI commands."""

from pathlib import Path
from typing import Tuple

import numpy as np
import pytest

pytest.importorskip("getdist")
pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table
from typer.testing import CliRunner

from numcosmo_py import Ncm
from numcosmo_py.app import app
from numcosmo_py.analysis.cluster_richness import RichnessModelType

pytestmark = pytest.mark.app
runner = CliRunner()

Ncm.cfg_init()


@pytest.fixture(name="mock_cluster_data")
def fixture_mock_cluster_data(tmp_path: Path) -> Tuple[Path, dict]:
    """Create a mock FITS file with cluster data.

    Creates a realistic but small dataset with:
    - 100 clusters
    - Masses between 1e14 and 1e15 M_sun
    - Redshifts between 0.1 and 0.5
    - Richnesses following a mass-richness relation with scatter
    """
    rng = np.random.RandomState(42)  # pylint: disable=no-member

    # Generate mock data
    n_clusters = 100
    log_mass = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
    redshift = rng.uniform(0.1, 0.5, n_clusters)

    # Mass-richness relation: ln(R) = ln(R_0) + alpha * ln(M/M_0)
    ln_R_0 = 3.0  # corresponds to R_0 ~ 20
    alpha = 1.2
    M_0 = 3e14
    sigma_lnR = 0.3

    ln_richness = (
        ln_R_0 + alpha * (log_mass - np.log(M_0)) + rng.normal(0, sigma_lnR, n_clusters)
    )

    # Ensure richness is reasonable (> 5)
    ln_richness = np.maximum(ln_richness, np.log(5))

    mass = np.exp(log_mass)
    richness = np.exp(ln_richness)
    richness_err = richness * 0.1  # 10% error

    # Create FITS table
    table = Table()
    table["halo_mass"] = mass
    table["redshift"] = redshift
    table["richness"] = richness
    table["richness_err"] = richness_err

    fits_file = tmp_path / "mock_clusters.fits"
    table.write(fits_file, format="fits", overwrite=True)

    metadata = {
        "n_clusters": n_clusters,
        "mass_range": (mass.min(), mass.max()),
        "z_range": (redshift.min(), redshift.max()),
        "richness_range": (richness.min(), richness.max()),
    }

    return fits_file, metadata


@pytest.fixture(name="mock_cluster_data_custom_columns")
def fixture_mock_cluster_data_custom_columns(tmp_path: Path) -> Path:
    """Create a mock FITS file with custom column names."""
    rng = np.random.RandomState(123)  # pylint: disable=no-member

    n_clusters = 50
    log_mass = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
    redshift = rng.uniform(0.1, 0.5, n_clusters)

    # Simple mass-richness relation
    ln_richness = 3.0 + 1.2 * (log_mass - np.log(3e14)) + rng.normal(0, 0.3, n_clusters)
    ln_richness = np.maximum(ln_richness, np.log(5))

    mass = np.exp(log_mass)
    richness = np.exp(ln_richness)
    richness_err = richness * 0.15

    # Use custom column names
    table = Table()
    table["M500c"] = mass
    table["z_spec"] = redshift
    table["lambda"] = richness
    table["lambda_sigma"] = richness_err

    fits_file = tmp_path / "mock_clusters_custom.fits"
    table.write(fits_file, format="fits", overwrite=True)

    return fits_file


@pytest.fixture(name="model_type", params=[e.value for e in RichnessModelType])
def fixture_model_type(request: pytest.FixtureRequest) -> str:
    """Return model type."""
    return request.param


class TestBasicDataLoading:
    """Test basic data loading and summary display."""

    def test_load_default_columns(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test loading data with default column names."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app, ["analysis", "cluster-richness", fits_file.as_posix()]
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Data Summary" in result.output
        assert "Total clusters:" in result.output
        assert "Mass range:" in result.output
        assert "Redshift range:" in result.output
        assert "Richness range:" in result.output

    def test_load_custom_columns(self, mock_cluster_data_custom_columns: Path) -> None:
        """Test loading data with custom column names."""
        fits_file = mock_cluster_data_custom_columns
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--mass-col",
                "M500c",
                "--redshift-col",
                "z_spec",
                "--richness-col",
                "lambda",
                "--sigma-lnR-col",
                "lambda_sigma",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Data Summary" in result.output
        assert "Total clusters:" in result.output

    def test_missing_column_error(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test error handling for missing column."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--mass-col",
                "nonexistent_column",
            ],
        )

        assert result.exit_code == 1, "Should fail with missing column"
        assert "not found" in result.output.lower() or "error" in result.output.lower()


class TestModelTypes:
    """Test different richness model types."""

    def test_ascaso_model(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with ASCASO model."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--model-type",
                "ascaso",
                "--n-bootstrap",
                "0",
                "--cuts",
                "10",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 1: Real Data Analysis" in result.output

    def test_ext_model(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with EXT model."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--model-type",
                "ext",
                "--n-bootstrap",
                "0",
                "--cuts",
                "10",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 1: Real Data Analysis" in result.output

    def test_all_model_types(self, mock_cluster_data: Tuple[Path, dict], model_type: str) -> None:
        """Test all available model types."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--model-type",
                model_type,
                "--n-bootstrap",
                "0",
                "--cuts",
                "10",
            ],
        )

        assert (
            result.exit_code == 0
        ), f"Command failed with {model_type}: {result.output}"


class TestCutAnalysis:
    """Test richness cut analysis."""

    def test_single_cut(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with a single richness cut."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--cuts",
                "20",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "λ ≥ 20.0" in result.output

    def test_multiple_cuts(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with multiple richness cuts."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--cuts",
                "5,10,15,20",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Check that all cuts appear in output
        for cut_val in [5.0, 10.0, 15.0, 20.0]:
            assert f"λ ≥ {cut_val}" in result.output

    def test_default_cuts(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with default cuts."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Default cuts are 5,10,15,20,30
        for cut_val in [5.0, 10.0, 15.0, 20.0, 30.0]:
            assert f"λ ≥ {cut_val}" in result.output


class TestAnalysisOptions:
    """Test various analysis options."""

    def test_ignore_noise(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis ignoring richness uncertainties."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--ignore-noise",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 1: Real Data Analysis" in result.output

    def test_custom_seed(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test analysis with custom random seed."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--seed",
                "12345",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

    def test_custom_output_prefix(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test analysis with custom output prefix."""
        fits_file, _ = mock_cluster_data

        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--output-prefix",
                "test_output",
                "--output-dir",
                tmp_path.as_posix(),
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

    def test_custom_hdu(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test loading from specific HDU."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--hdu",
                "1",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

    def test_bootstrap_disabled(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test with bootstrap disabled."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--no-bootstrap",
                "--cuts",
                "10",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

    def test_bootstrap_enabled_small(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test with small number of bootstrap samples."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--bootstrap",
                "--n-bootstrap",
                "10",
                "--cuts",
                "10",
                "--output-dir",
                tmp_path.as_posix(),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"


class TestMockStudy:
    """Test mock study functionality."""

    def test_mock_study_minimal(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test mock study with minimal number of mocks."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-mocks",
                "--n-mocks",
                "10",  # Very small for speed
                "--n-bootstrap",
                "0",
                "--cuts",
                "10",
                "--output-dir",
                tmp_path.as_posix(),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 2: Mock Study" in result.output
        assert "PHASE 3: Goodness-of-Fit Analysis" in result.output

    def test_mock_study_zero_mocks(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test that zero mocks skips mock study."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-mocks",
                "--n-mocks",
                "0",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
                "--output-dir",
                tmp_path.as_posix(),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Skipping mock study" in result.output

    def test_mock_study_requires_analysis(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test that mock study runs only when analysis runs."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-mocks",
                "--n-mocks",
                "5",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
                "--output-dir",
                tmp_path.as_posix(),
            ],
        )

        # Should run analysis and mocks
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 1: Real Data Analysis" in result.output
        assert "PHASE 2: Mock Study" in result.output


class TestDiagnostics:
    """Test diagnostic analysis functionality."""

    def test_diagnostics_basic(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test basic diagnostic analysis."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-diagnostics",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Diagnostic Analysis" in result.output

    def test_diagnostics_custom_bins(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test diagnostics with custom number of bins."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-diagnostics",
                "--n-bins",
                "10",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Diagnostic Analysis" in result.output

    def test_diagnostics_with_multiple_cuts(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test diagnostics with multiple cuts."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-diagnostics",
                "--cuts",
                "5,10,15",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Check diagnostics ran for multiple cuts
        assert "Diagnostic Analysis" in result.output


class TestIntegration:
    """Integration tests combining multiple features."""

    def test_full_pipeline_minimal(self, mock_cluster_data: Tuple[Path, dict], tmp_path: Path) -> None:
        """Test full pipeline with minimal settings."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--run-mocks",
                "--n-mocks",
                "5",
                "--run-diagnostics",
                "--cuts",
                "10,15",
                "--n-bootstrap",
                "0",
                "--output-dir",
                tmp_path.as_posix(),
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "PHASE 1: Real Data Analysis" in result.output
        assert "PHASE 2: Mock Study" in result.output
        assert "PHASE 3: Goodness-of-Fit Analysis" in result.output
        assert "Diagnostic Analysis" in result.output

    def test_reproducibility(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test that same seed produces same results."""
        fits_file, _ = mock_cluster_data

        # Run twice with same seed
        results = []
        for _ in range(2):
            result = runner.invoke(
                app,
                [
                    "analysis",
                    "cluster-richness",
                    fits_file.as_posix(),
                    "--seed",
                    "424242",
                    "--cuts",
                    "10",
                    "--n-bootstrap",
                    "0",
                ],
            )
            results.append(result)

        # Both should succeed
        assert all(r.exit_code == 0 for r in results)
        # Outputs should be identical (parameter values)
        assert results[0].output == results[1].output

    def test_multiple_model_types_same_data(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test that different models can analyze same data."""
        fits_file, _ = mock_cluster_data

        for model_type in ["ascaso", "ext"]:
            result = runner.invoke(
                app,
                [
                    "analysis",
                    "cluster-richness",
                    fits_file.as_posix(),
                    "--model-type",
                    model_type,
                    "--cuts",
                    "10",
                    "--n-bootstrap",
                    "0",
                ],
            )

            assert result.exit_code == 0, f"Failed for {model_type}: {result.output}"
            assert "Best-fit" in result.output


class TestEdgeCases:
    """Test edge cases and boundary conditions."""

    def test_very_small_dataset(self, tmp_path: Path) -> None:
        """Test with minimal dataset (10 clusters)."""
        rng = np.random.RandomState(999)  # pylint: disable=no-member

        n_clusters = 10
        log_mass = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
        mass = np.exp(log_mass)
        redshift = rng.uniform(0.2, 0.4, n_clusters)
        richness = 20.0 + rng.normal(0, 5, n_clusters)
        richness = np.maximum(richness, 5.0)
        richness_err = richness * 0.1

        table = Table()
        table["halo_mass"] = mass
        table["redshift"] = redshift
        table["richness"] = richness
        table["richness_err"] = richness_err

        fits_file = tmp_path / "small_dataset.fits"
        table.write(fits_file, format="fits", overwrite=True)

        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        # Should still work with small dataset
        assert result.exit_code == 0, f"Command failed: {result.output}"
        assert "Total clusters: 10" in result.output

    def test_high_richness_cut(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test with high richness cut that filters most data."""
        fits_file, metadata = mock_cluster_data
        max_richness = metadata["richness_range"][1]

        # Use cut close to maximum richness
        high_cut = max_richness * 0.9

        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--cuts",
                f"{high_cut:.0f}",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"

    def test_consistent_column_naming(self, tmp_path: Path) -> None:
        """Test that column name variations work correctly."""
        rng = np.random.RandomState(555)  # pylint: disable=no-member

        n_clusters = 50
        log_mass = rng.uniform(np.log(1e14), np.log(1e15), n_clusters)
        mass = np.exp(log_mass)
        redshift = rng.uniform(0.2, 0.4, n_clusters)
        richness = 20.0 + rng.normal(0, 5, n_clusters)
        richness = np.maximum(richness, 5.0)

        table = Table()
        table["M"] = mass
        table["z"] = redshift
        table["R"] = richness
        table["R_err"] = richness * 0.1

        fits_file = tmp_path / "simple_columns.fits"
        table.write(fits_file, format="fits", overwrite=True)

        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--mass-col",
                "M",
                "--redshift-col",
                "z",
                "--richness-col",
                "R",
                "--sigma-lnR-col",
                "R_err",
                "--cuts",
                "10",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"


class TestDataValidation:
    """Test data validation and quality checks."""

    def test_cluster_counts_per_cut(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test that cluster counts are reported correctly."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--cuts",
                "5,10,20",
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Should show cluster counts for each cut
        assert "Clusters per Richness Cut" in result.output
        assert "N clusters" in result.output
        assert "Fraction" in result.output

    def test_data_summary_ranges(self, mock_cluster_data: Tuple[Path, dict]) -> None:
        """Test that data summary includes all range information."""
        fits_file, _ = mock_cluster_data
        result = runner.invoke(
            app,
            [
                "analysis",
                "cluster-richness",
                fits_file.as_posix(),
                "--n-bootstrap",
                "0",
            ],
        )

        assert result.exit_code == 0, f"Command failed: {result.output}"
        # Check all expected ranges are reported
        assert "Mass range:" in result.output
        assert "Redshift range:" in result.output
        assert "Richness range:" in result.output
        assert "σ(ln R) range:" in result.output
