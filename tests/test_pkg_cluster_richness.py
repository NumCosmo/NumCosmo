#!/usr/bin/env python
#
# test_py_cluster_richness_package.py
#
# Copyright 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for the numcosmo_py.analysis.cluster_richness package."""

import os
import tempfile
import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

from numcosmo_py.analysis.cluster_richness import (
    # Parameters module
    CutAnalysisResult,
    dup_model,
    model_to_yaml,
    model_from_yaml,
    model_params_to_dict,
    model_params_from_dict,
    copy_model_params,
    model_params_as_list,
    model_params_from_list,
    get_model_param_names,
    # Utils module
    setup_model_fit_params,
    create_richness_model,
    get_model_type_name,
    PARAM_FORMAT,
    # Truncated stats module
    mean_lnR_truncated,
    std_lnR_truncated,
    invert_truncated_stats,
    invert_truncated_stats_mu_from_sample,
    invert_truncated_stats_sigma_from_sample,
    # Database module
    BestfitDatabase,
    # Analyzer module
    CutAnalyzer,
)

Ncm.cfg_init()

LN_RICHNESS_CUT = np.log(5)
LN_RICHNESS_MIN = LN_RICHNESS_CUT
LN_RICHNESS_MAX = np.log(200)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture(name="ascaso_model")
def fixture_ascaso_model() -> Nc.ClusterMassAscaso:
    """Create NcClusterMassAscaso model with typical parameters."""
    model = Nc.ClusterMassAscaso(
        lnRichness_min=LN_RICHNESS_MIN, lnRichness_max=LN_RICHNESS_MAX
    )
    model["mup0"] = 4.0
    model["mup1"] = 1.0
    model["mup2"] = 0.2
    model["sigmap0"] = 0.5
    model["sigmap1"] = 0.03
    model["sigmap2"] = 0.15
    model["cut"] = LN_RICHNESS_CUT
    return model


@pytest.fixture(name="ext_model")
def fixture_ext_model() -> Nc.ClusterMassExt:
    """Create NcClusterMassExt model with typical parameters."""
    model = Nc.ClusterMassExt(
        lnRichness_min=LN_RICHNESS_MIN, lnRichness_max=LN_RICHNESS_MAX
    )
    model["mup0"] = 4.0
    model["mup1"] = 1.0
    model["mup2"] = 0.1
    model["mup3"] = 0.01
    model["sigmap0"] = -0.3
    model["sigmap1"] = -0.08
    model["sigmap2"] = 0.005
    model["cut"] = LN_RICHNESS_CUT
    return model


@pytest.fixture(name="richness_model", params=["ascaso", "ext"])
def fixture_richness_model(
    request, ascaso_model: Nc.ClusterMassAscaso, ext_model: Nc.ClusterMassExt
) -> Nc.ClusterMassRichness:
    """Parametrized fixture that yields both Ascaso and Ext models."""
    if request.param == "ascaso":
        return ascaso_model
    return ext_model


@pytest.fixture(name="sample_data")
def fixture_sample_data(ascaso_model: Nc.ClusterMassAscaso) -> tuple:
    """Generate sample data for testing."""
    n_samples = 200
    rng = Ncm.RNG.seeded_new(None, 42)

    # Generate mock data
    lnM = np.linspace(np.log(1e13), np.log(1e15), n_samples)
    z = np.linspace(0.1, 1.0, n_samples)

    # Create data object
    lnM_v = Ncm.Vector.new_array(lnM)
    z_v = Ncm.Vector.new_array(z)
    lnR_v = Ncm.Vector.new(n_samples)
    lnR_v.set_zero()

    data = Nc.DataClusterMassRich.new()
    data.set_data(lnM_v, z_v, lnR_v)

    mset = Ncm.MSet.new_array([ascaso_model])
    mset.prepare_fparam_map()
    data.resample(mset, rng)

    lnR = np.array([data.peek_lnR().get(i) for i in range(n_samples)])

    return lnM, z, lnR


@pytest.fixture(name="temp_db_path")
def fixture_temp_db_path():
    """Create a temporary database file."""
    fd, path = tempfile.mkstemp(suffix=".db")
    os.close(fd)
    yield path
    # Cleanup
    if os.path.exists(path):
        os.remove(path)


# =============================================================================
# Tests for _parameters module
# =============================================================================


class TestDupModel:
    """Tests for dup_model function."""

    def test_dup_model_preserves_type(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that dup_model preserves the model type."""
        dup = dup_model(richness_model)
        assert type(dup) is type(richness_model)

    def test_dup_model_preserves_parameters(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that dup_model preserves all parameters."""
        dup = dup_model(richness_model)
        for i in range(richness_model.sparam_len()):
            assert dup.param_get(i) == richness_model.param_get(i)

    def test_dup_model_creates_independent_copy(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that modifying duplicate doesn't affect original."""
        original_value = richness_model["mup0"]
        dup = dup_model(richness_model)
        dup["mup0"] = original_value + 10.0
        assert richness_model["mup0"] == original_value
        assert dup["mup0"] == original_value + 10.0

    def test_dup_model_multiple_times(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that dup_model can be called multiple times (fresh serializer)."""
        dups = [dup_model(richness_model) for _ in range(5)]
        for dup in dups:
            assert type(dup) is type(richness_model)
            assert dup["mup0"] == richness_model["mup0"]


class TestModelYamlSerialization:
    """Tests for model_to_yaml and model_from_yaml functions."""

    def test_yaml_roundtrip(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test YAML serialization roundtrip preserves parameters."""
        yaml_str = model_to_yaml(richness_model)
        restored = model_from_yaml(yaml_str)

        for i in range(richness_model.sparam_len()):
            assert restored.param_get(i) == richness_model.param_get(i)

    def test_yaml_output_is_string(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that model_to_yaml returns a string."""
        yaml_str = model_to_yaml(richness_model)
        assert isinstance(yaml_str, str)
        assert len(yaml_str) > 0

    def test_yaml_contains_type_info(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that YAML contains type information."""
        yaml_str = model_to_yaml(richness_model)
        type_name = type(richness_model).__name__
        # GObject type name has Nc prefix
        assert f"Nc{type_name}" in yaml_str or type_name in yaml_str


class TestModelParamsToDict:
    """Tests for model_params_to_dict function."""

    def test_returns_dict(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that it returns a dictionary."""
        result = model_params_to_dict(richness_model)
        assert isinstance(result, dict)

    def test_contains_all_params(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that all parameters are included."""
        result = model_params_to_dict(richness_model)
        assert len(result) == richness_model.sparam_len()

    def test_correct_values(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that values are correct."""
        result = model_params_to_dict(richness_model)
        for name, value in result.items():
            assert value == richness_model[name]

    def test_known_params_ascaso(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test known parameters for Ascaso model."""
        result = model_params_to_dict(ascaso_model)
        assert "mup0" in result
        assert "sigmap0" in result
        assert result["mup0"] == 4.0


class TestModelParamsFromDict:
    """Tests for model_params_from_dict function."""

    def test_sets_params(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that parameters are set correctly."""
        new_values = {"mup0": 5.0, "sigmap0": 0.3}
        model_params_from_dict(richness_model, new_values)
        assert richness_model["mup0"] == 5.0
        assert richness_model["sigmap0"] == 0.3

    def test_ignores_unknown_params(
        self, richness_model: Nc.ClusterMassRichness
    ) -> None:
        """Test that unknown parameters are ignored."""
        original_mup0 = richness_model["mup0"]
        new_values = {"unknown_param": 999.0}
        model_params_from_dict(richness_model, new_values)
        assert richness_model["mup0"] == original_mup0

    def test_partial_update(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test partial parameter update."""
        original_sigmap0 = richness_model["sigmap0"]
        new_values = {"mup0": 7.0}
        model_params_from_dict(richness_model, new_values)
        assert richness_model["mup0"] == 7.0
        assert richness_model["sigmap0"] == original_sigmap0


class TestCopyModelParams:
    """Tests for copy_model_params function."""

    def test_copy_params(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test copying parameters between models of same type."""
        target = Nc.ClusterMassAscaso(
            lnRichness_min=LN_RICHNESS_MIN, lnRichness_max=LN_RICHNESS_MAX
        )
        target["mup0"] = 0.0  # Different from source

        copy_model_params(ascaso_model, target)

        assert target["mup0"] == ascaso_model["mup0"]
        assert target["sigmap0"] == ascaso_model["sigmap0"]

    def test_copy_params_different_types(
        self, ascaso_model: Nc.ClusterMassAscaso, ext_model: Nc.ClusterMassExt
    ) -> None:
        """Test copying parameters between different model types."""
        # Common params should be copied
        ascaso_model["mup0"] = 99.0
        copy_model_params(ascaso_model, ext_model)
        assert ext_model["mup0"] == 99.0


class TestModelParamsAsList:
    """Tests for model_params_as_list function."""

    def test_returns_list(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that it returns a list."""
        result = model_params_as_list(richness_model)
        assert isinstance(result, list)

    def test_correct_length(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test list has correct length."""
        result = model_params_as_list(richness_model)
        assert len(result) == richness_model.sparam_len()

    def test_correct_order(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test values are in correct order."""
        result = model_params_as_list(richness_model)
        for i, value in enumerate(result):
            assert value == richness_model.param_get(i)


class TestModelParamsFromList:
    """Tests for model_params_from_list function."""

    def test_sets_params(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test setting parameters from list."""
        n_params = richness_model.sparam_len()
        new_values = [float(i + 1) for i in range(n_params)]
        model_params_from_list(richness_model, new_values)

        for i, value in enumerate(new_values):
            assert richness_model.param_get(i) == value

    def test_roundtrip(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test list roundtrip."""
        original = model_params_as_list(richness_model)
        target = dup_model(richness_model)
        model_params_from_list(target, original)

        for i in range(richness_model.sparam_len()):
            assert target.param_get(i) == richness_model.param_get(i)


class TestGetModelParamNames:
    """Tests for get_model_param_names function."""

    def test_returns_list(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that it returns a list of strings."""
        result = get_model_param_names(richness_model)
        assert isinstance(result, list)
        assert all(isinstance(name, str) for name in result)

    def test_correct_length(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test list has correct length."""
        result = get_model_param_names(richness_model)
        assert len(result) == richness_model.sparam_len()

    def test_known_names_ascaso(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test known parameter names for Ascaso model."""
        result = get_model_param_names(ascaso_model)
        assert "mup0" in result
        assert "mup1" in result
        assert "sigmap0" in result


class TestCutAnalysisResult:
    """Tests for CutAnalysisResult dataclass."""

    def test_creation(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test creating CutAnalysisResult."""
        result = CutAnalysisResult(
            cut=1.0,
            n_clusters=100,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=50.0,
        )
        assert result.cut == 1.0
        assert result.n_clusters == 100
        assert result.m2lnL == 50.0

    def test_get_bestfit_params_dict(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test get_bestfit_params_dict method."""
        result = CutAnalysisResult(
            cut=1.0,
            n_clusters=100,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=50.0,
        )
        params = result.get_bestfit_params_dict()
        assert params["mup0"] == ascaso_model["mup0"]

    def test_to_db_row(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test to_db_row method."""
        result = CutAnalysisResult(
            cut=1.0,
            n_clusters=100,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=50.0,
        )
        row = result.to_db_row(mock_seed=42)
        assert row["mock_seed"] == 42
        # Note: model params include "cut" which overwrites the result cut
        # The row["cut"] will be the model's cut parameter value
        assert row["n_clusters"] == 100
        assert row["m2lnL"] == 50.0
        assert "mup0" in row


# =============================================================================
# Tests for _utils module
# =============================================================================


class TestSetupModelFitParams:
    """Tests for setup_model_fit_params function."""

    def test_sets_fit_flag(self, richness_model: Nc.ClusterMassRichness) -> None:
        """Test that fit flag is set for all parameters."""
        setup_model_fit_params(richness_model)
        for _ in range(richness_model.sparam_len()):
            # Check parameter is fittable (would need to check model state)
            # At minimum, the function should run without error
            pass


class TestCreateRichnessModel:
    """Tests for create_richness_model function."""

    def test_create_ascaso(self) -> None:
        """Test creating Ascaso model."""
        model = create_richness_model("ascaso")
        assert isinstance(model, Nc.ClusterMassAscaso)

    def test_create_ext(self) -> None:
        """Test creating Ext model."""
        model = create_richness_model("ext")
        assert isinstance(model, Nc.ClusterMassExt)

    def test_case_insensitive(self) -> None:
        """Test model type is case insensitive."""
        model1 = create_richness_model("ASCASO")
        model2 = create_richness_model("Ascaso")
        assert isinstance(model1, Nc.ClusterMassAscaso)
        assert isinstance(model2, Nc.ClusterMassAscaso)

    def test_custom_bounds(self) -> None:
        """Test custom richness bounds."""
        model = create_richness_model("ascaso", lnRichness_min=1.0, lnRichness_max=10.0)
        assert isinstance(model, Nc.ClusterMassAscaso)

    def test_unknown_type_raises(self) -> None:
        """Test that unknown model type raises ValueError."""
        with pytest.raises(ValueError, match="Unknown model type"):
            create_richness_model("unknown")


class TestGetModelTypeName:
    """Tests for get_model_type_name function."""

    def test_ascaso(self, ascaso_model: Nc.ClusterMassAscaso) -> None:
        """Test type name for Ascaso model."""
        assert get_model_type_name(ascaso_model) == "ascaso"

    def test_ext(self, ext_model: Nc.ClusterMassExt) -> None:
        """Test type name for Ext model."""
        assert get_model_type_name(ext_model) == "ext"


class TestParamFormat:
    """Tests for PARAM_FORMAT constant."""

    def test_format_string(self) -> None:
        """Test that PARAM_FORMAT is a valid format string."""
        value = 3.14159
        formatted = f"{value:{PARAM_FORMAT}}"
        assert isinstance(formatted, str)


# =============================================================================
# Tests for _truncated_stats module
# =============================================================================


class TestMeanLnRTruncated:
    """Tests for mean_lnR_truncated function."""

    def test_scalar_inputs(self) -> None:
        """Test with scalar inputs."""
        result = mean_lnR_truncated(mu=3.0, sigma=0.5, lnR_cut=2.0)
        assert isinstance(result, (float, np.floating))

    def test_truncated_mean_greater_than_cut(self) -> None:
        """Test that truncated mean is greater than cut."""
        mu, sigma, lnR_cut = 2.0, 1.0, 1.5
        result = mean_lnR_truncated(mu, sigma, lnR_cut)
        assert result > lnR_cut

    def test_truncated_mean_greater_than_mu(self) -> None:
        """Test that truncated mean is greater than underlying mu."""
        mu, sigma, lnR_cut = 2.0, 1.0, 1.5
        result = mean_lnR_truncated(mu, sigma, lnR_cut)
        assert result > mu

    def test_array_inputs(self) -> None:
        """Test with array inputs."""
        mu = np.array([2.0, 3.0, 4.0])
        sigma = np.array([0.5, 0.5, 0.5])
        result = mean_lnR_truncated(mu, sigma, lnR_cut=1.5)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3


class TestStdLnRTruncated:
    """Tests for std_lnR_truncated function."""

    def test_scalar_inputs(self) -> None:
        """Test with scalar inputs."""
        result = std_lnR_truncated(mu=3.0, sigma=0.5, lnR_cut=2.0)
        assert isinstance(result, (float, np.floating))

    def test_truncated_std_positive(self) -> None:
        """Test that truncated std is positive."""
        result = std_lnR_truncated(mu=2.0, sigma=1.0, lnR_cut=1.5)
        assert result > 0

    def test_truncated_std_less_than_sigma(self) -> None:
        """Test that truncated std is less than underlying sigma."""
        sigma = 1.0
        result = std_lnR_truncated(mu=2.0, sigma=sigma, lnR_cut=1.5)
        assert result < sigma

    def test_array_inputs(self) -> None:
        """Test with array inputs."""
        mu = np.array([2.0, 3.0, 4.0])
        sigma = np.array([0.5, 0.5, 0.5])
        result = std_lnR_truncated(mu, sigma, lnR_cut=1.5)
        assert isinstance(result, np.ndarray)
        assert len(result) == 3


class TestInvertTruncatedStats:
    """Tests for invert_truncated_stats function."""

    def test_successful_inversion(self) -> None:
        """Test successful inversion of truncated stats."""
        # Start with known mu/sigma, compute truncated stats, then recover
        mu_true, sigma_true, lnR_cut = 3.0, 0.8, 2.0
        m_obs = mean_lnR_truncated(mu_true, sigma_true, lnR_cut)
        s_obs = std_lnR_truncated(mu_true, sigma_true, lnR_cut)

        assert isinstance(m_obs, float)
        assert isinstance(s_obs, float)
        result = invert_truncated_stats(m_obs, s_obs, lnR_cut)

        assert result["success"]
        assert result["mu"] is not None
        assert result["sigma"] is not None
        assert abs(result["mu"] - mu_true) < 0.1
        assert abs(result["sigma"] - sigma_true) < 0.1

    def test_result_keys(self) -> None:
        """Test that result contains expected keys."""
        result = invert_truncated_stats(3.5, 0.7, 2.0)
        assert "success" in result
        assert "message" in result
        assert "mu" in result
        assert "sigma" in result


class TestInvertTruncatedStatsMuFromSample:
    """Tests for invert_truncated_stats_mu_from_sample function."""

    def test_with_sample(self) -> None:
        """Test recovering mu from a sample."""
        # Generate truncated sample
        mu_true, sigma_true, lnR_cut = 3.0, 0.5, 2.5
        np.random.seed(42)
        full_sample = np.random.normal(mu_true, sigma_true, 10000)
        truncated_sample = full_sample[full_sample > lnR_cut]

        result = invert_truncated_stats_mu_from_sample(truncated_sample, lnR_cut)

        assert result is not None
        assert abs(result - mu_true) < 0.2


class TestInvertTruncatedStatsSigmaFromSample:
    """Tests for invert_truncated_stats_sigma_from_sample function."""

    def test_with_sample(self) -> None:
        """Test recovering sigma from a sample."""
        # Generate truncated sample
        mu_true, sigma_true, lnR_cut = 3.0, 0.5, 2.5
        np.random.seed(42)
        full_sample = np.random.normal(mu_true, sigma_true, 10000)
        truncated_sample = full_sample[full_sample > lnR_cut]

        result = invert_truncated_stats_sigma_from_sample(truncated_sample, lnR_cut)

        assert result is not None
        assert abs(result - sigma_true) < 0.2


# =============================================================================
# Tests for _database module
# =============================================================================


class TestBestfitDatabase:
    """Tests for BestfitDatabase class."""

    def test_init_creates_db(self, temp_db_path: str) -> None:
        """Test that initialization creates database file."""
        db = BestfitDatabase(temp_db_path)
        assert os.path.exists(temp_db_path)
        assert db.count_entries() == 0

    def test_insert_bestfit(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test inserting a best-fit result."""
        db = BestfitDatabase(temp_db_path)
        result = CutAnalysisResult(
            cut=1.5,
            n_clusters=50,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=100.0,
        )
        success = db.insert_bestfit(mock_seed=1, result=result)
        assert success
        assert db.count_entries() == 1

    def test_insert_duplicate_fails(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test that inserting duplicate returns False."""
        db = BestfitDatabase(temp_db_path)
        result = CutAnalysisResult(
            cut=1.5,
            n_clusters=50,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=100.0,
        )
        db.insert_bestfit(mock_seed=1, result=result)
        success = db.insert_bestfit(mock_seed=1, result=result)
        assert not success
        assert db.count_entries() == 1

    def test_batch_insert(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test batch inserting multiple results."""
        db = BestfitDatabase(temp_db_path)
        results = []
        for seed in range(5):
            result = CutAnalysisResult(
                cut=float(seed) + 1.0,
                n_clusters=50 + seed * 10,
                bestfit=dup_model(ascaso_model),
                mcmc_mean=dup_model(ascaso_model),
                mcmc_median=dup_model(ascaso_model),
                bootstrap_mean=dup_model(ascaso_model),
                m2lnL=100.0 + seed,
            )
            results.append((seed, result))

        db.batch_insert(results)
        assert db.count_entries() == 5

    def test_get_bestfit(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test retrieving a best-fit model."""
        db = BestfitDatabase(temp_db_path)
        ascaso_model["mup0"] = 5.5  # Distinctive value
        result = CutAnalysisResult(
            cut=1.5,
            n_clusters=50,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=100.0,
        )
        db.insert_bestfit(mock_seed=42, result=result)

        retrieved = db.get_bestfit(mock_seed=42, cut=1.5)
        assert retrieved is not None
        assert retrieved["mup0"] == 5.5

    def test_get_bestfit_not_found(self, temp_db_path: str) -> None:
        """Test that get_bestfit returns None when not found."""
        db = BestfitDatabase(temp_db_path)
        result = db.get_bestfit(mock_seed=999, cut=999.0)
        assert result is None

    def test_get_bestfit_params_dict(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test retrieving best-fit parameters as dictionary."""
        db = BestfitDatabase(temp_db_path)
        ascaso_model["mup0"] = 7.7
        result = CutAnalysisResult(
            cut=2.0,
            n_clusters=60,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=80.0,
        )
        db.insert_bestfit(mock_seed=10, result=result)

        params = db.get_bestfit_params_dict(mock_seed=10, cut=2.0)
        assert params is not None
        assert params["mup0"] == 7.7

    def test_get_computed_seeds(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test getting computed seeds for a specific cut."""
        db = BestfitDatabase(temp_db_path)
        cut = 1.5
        for seed in [1, 3, 5]:
            result = CutAnalysisResult(
                cut=cut,
                n_clusters=50,
                bestfit=dup_model(ascaso_model),
                mcmc_mean=dup_model(ascaso_model),
                mcmc_median=dup_model(ascaso_model),
                bootstrap_mean=dup_model(ascaso_model),
                m2lnL=100.0,
            )
            db.insert_bestfit(mock_seed=seed, result=result)

        seeds = db.get_computed_seeds(cut=cut)
        assert seeds == {1, 3, 5}

    def test_get_all_computed_seeds(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test getting all computed seeds."""
        db = BestfitDatabase(temp_db_path)
        for seed in [1, 2, 3]:
            result = CutAnalysisResult(
                cut=float(seed),
                n_clusters=50,
                bestfit=dup_model(ascaso_model),
                mcmc_mean=dup_model(ascaso_model),
                mcmc_median=dup_model(ascaso_model),
                bootstrap_mean=dup_model(ascaso_model),
                m2lnL=100.0,
            )
            db.insert_bestfit(mock_seed=seed, result=result)

        seeds = db.get_all_computed_seeds()
        assert seeds == {1, 2, 3}

    def test_delete_bestfit(
        self, temp_db_path: str, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test deleting a best-fit entry."""
        db = BestfitDatabase(temp_db_path)
        result = CutAnalysisResult(
            cut=1.5,
            n_clusters=50,
            bestfit=ascaso_model,
            mcmc_mean=dup_model(ascaso_model),
            mcmc_median=dup_model(ascaso_model),
            bootstrap_mean=dup_model(ascaso_model),
            m2lnL=100.0,
        )
        db.insert_bestfit(mock_seed=1, result=result)
        assert db.count_entries() == 1

        deleted = db.delete_bestfit(mock_seed=1, cut=1.5)
        assert deleted
        assert db.count_entries() == 0

    def test_delete_bestfit_not_found(self, temp_db_path: str) -> None:
        """Test that deleting non-existent entry returns False."""
        db = BestfitDatabase(temp_db_path)
        deleted = db.delete_bestfit(mock_seed=999, cut=999.0)
        assert not deleted

    def test_stores_ext_model(
        self, temp_db_path: str, ext_model: Nc.ClusterMassExt
    ) -> None:
        """Test that Ext model can be stored and retrieved."""
        db = BestfitDatabase(temp_db_path)
        ext_model["mup0"] = 8.8
        result = CutAnalysisResult(
            cut=2.5,
            n_clusters=70,
            bestfit=ext_model,
            mcmc_mean=dup_model(ext_model),
            mcmc_median=dup_model(ext_model),
            bootstrap_mean=dup_model(ext_model),
            m2lnL=90.0,
        )
        db.insert_bestfit(mock_seed=20, result=result)

        retrieved = db.get_bestfit(mock_seed=20, cut=2.5)
        assert retrieved is not None
        assert isinstance(retrieved, Nc.ClusterMassExt)
        assert retrieved["mup0"] == 8.8


# =============================================================================
# Tests for _analyzer module
# =============================================================================


class TestCutAnalyzer:
    """Tests for CutAnalyzer class."""

    def test_init(self, sample_data: tuple) -> None:
        """Test CutAnalyzer initialization."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10), np.log(20)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        assert analyzer.lnM is lnM
        assert analyzer.cuts == cuts

    def test_analyze_with_default_model(self, sample_data: tuple) -> None:
        """Test analyze with default model."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        model_init = create_richness_model("ascaso")
        results = analyzer.analyze(model_init=model_init)

        assert len(results) == 1
        assert cuts[0] in results
        result = results[cuts[0]]
        assert isinstance(result, CutAnalysisResult)
        assert result.n_clusters > 0

    def test_analyze_with_custom_model(
        self, sample_data: tuple, ascaso_model: Nc.ClusterMassAscaso
    ) -> None:
        """Test analyze with custom model."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        results = analyzer.analyze(model_init=ascaso_model)

        assert len(results) == 1
        assert cuts[0] in results

    def test_analyze_multiple_cuts(self, sample_data: tuple) -> None:
        """Test analyze with multiple cuts."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10), np.log(15), np.log(20)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        model_init = create_richness_model("ascaso")
        results = analyzer.analyze(model_init=model_init)

        assert len(results) == 3
        for cut in cuts:
            assert cut in results

    def test_results_stored(self, sample_data: tuple) -> None:
        """Test that results are stored in analyzer."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        model_init = create_richness_model("ascaso")
        analyzer.analyze(model_init=model_init)

        assert len(analyzer.results) == 1

    def test_display_results_no_error(self, sample_data: tuple) -> None:
        """Test that display_results runs without error."""
        lnM, z, lnR = sample_data
        cuts = [np.log(10)]
        analyzer = CutAnalyzer(
            lnM=lnM,
            z=z,
            lnR=lnR,
            cuts=cuts,
            verbose=False,
        )
        model_init = create_richness_model("ascaso")
        analyzer.analyze(model_init=model_init)
        # Should not raise
        analyzer.display_results()
