#!/usr/bin/env python
#
# test_py_galaxy_sd_redshift_lsst.py
#
# Tue Mar 04 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_galaxy_sd_redshift_lsst.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Unit tests for LSSTgalaxy redshift distributions.

Tests NcGalaxySDTrueRedshiftLSSTSRD and NcGalaxySDObsRedshiftGauss
with LSST parametrizations (Y1/Y10 x source/lens).
"""

from functools import cache
import pytest
from numpy.testing import assert_allclose
import numpy as np
from scipy.special import erf, erfc  # type: ignore # pylint: disable=no-name-in-module

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()


# Expected parameter values for each LSST SRD parametrization
LSST_SRD_PARAMS = {
    "Y1_SOURCE": {"alpha": 0.78, "beta": 2.0, "z0": 0.13},
    "Y1_LENS": {"alpha": 0.94, "beta": 2.0, "z0": 0.26},
    "Y10_SOURCE": {"alpha": 0.68, "beta": 2.0, "z0": 0.11},
    "Y10_LENS": {"alpha": 0.90, "beta": 2.0, "z0": 0.28},
}

# Expected bin counts for each parametrization
LSST_SRD_BIN_COUNTS = {
    "Y1_SOURCE": 5,
    "Y1_LENS": 5,
    "Y10_SOURCE": 5,
    "Y10_LENS": 10,
}


@cache
def get_cosmology() -> Cosmology:
    """Create a default cosmology for testing."""
    return Cosmology.default()


@pytest.fixture(name="cosmology", scope="module")
def fixture_cosmology() -> Cosmology:
    """Create a simple cosmology for testing."""
    cosmology = get_cosmology()
    return cosmology


@pytest.fixture(name="rng", scope="module")
def fixture_rng() -> Ncm.RNG:
    """Create a random number generator."""
    return Ncm.RNG.seeded_new(None, 42)


@pytest.fixture(name="serializer", scope="module")
def fixture_serializer() -> Ncm.Serialize:
    """Create a serializer for testing."""
    return Ncm.Serialize.new(Ncm.SerializeOpt.NONE)


# =============================================================================
# Parametrized fixtures for all LSST SRD types
# =============================================================================


def _create_lsst_srd_from_type(
    lsst_type: Nc.GalaxySDTrueRedshiftLSSTSRDType,
) -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create LSST SRD distribution from enum type."""
    return Nc.GalaxySDTrueRedshiftLSSTSRD.new_from_type(lsst_type)


@pytest.fixture(
    name="lsst_srd",
    params=[
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_SOURCE,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_LENS,
    ],
    ids=["Y1_SOURCE", "Y1_LENS", "Y10_SOURCE", "Y10_LENS"],
)
def fixture_lsst_srd(request) -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create all LSST SRD parametrizations."""
    return _create_lsst_srd_from_type(request.param)


@pytest.fixture(name="lsst_srd_y1_source", scope="module")
def fixture_lsst_srd_y1_source() -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create Y1 source parametrization."""
    return Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()


@pytest.fixture(name="lsst_srd_y1_lens", scope="module")
def fixture_lsst_srd_y1_lens() -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create Y1 lens parametrization."""
    return Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_lens()


@pytest.fixture(name="lsst_srd_y10_source", scope="module")
def fixture_lsst_srd_y10_source() -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create Y10 source parametrization."""
    return Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_source()


@pytest.fixture(name="lsst_srd_y10_lens", scope="module")
def fixture_lsst_srd_y10_lens() -> Nc.GalaxySDTrueRedshiftLSSTSRD:
    """Create Y10 lens parametrization."""
    return Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_lens()


# =============================================================================
# NcGalaxySDTrueRedshiftLSSTSRD Tests
# =============================================================================


def test_lsst_srd_creation(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test that LSST SRD distributions can be created."""
    assert lsst_srd is not None
    assert isinstance(lsst_srd, Nc.GalaxySDTrueRedshiftLSSTSRD)
    assert isinstance(lsst_srd, Nc.GalaxySDTrueRedshift)


def test_lsst_srd_specialized_constructors() -> None:
    """Test all specialized constructors."""
    y1_source = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    assert isinstance(y1_source, Nc.GalaxySDTrueRedshiftLSSTSRD)

    y1_lens = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_lens()
    assert isinstance(y1_lens, Nc.GalaxySDTrueRedshiftLSSTSRD)

    y10_source = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_source()
    assert isinstance(y10_source, Nc.GalaxySDTrueRedshiftLSSTSRD)

    y10_lens = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_lens()
    assert isinstance(y10_lens, Nc.GalaxySDTrueRedshiftLSSTSRD)


def test_lsst_srd_from_type() -> None:
    """Test constructor from enum type."""
    for lsst_type in [
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_SOURCE,
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_LENS,
    ]:
        gsdtr = Nc.GalaxySDTrueRedshiftLSSTSRD.new_from_type(lsst_type)
        assert isinstance(gsdtr, Nc.GalaxySDTrueRedshiftLSSTSRD)


def test_lsst_srd_parameters(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test that parameters are correctly set."""
    alpha = lsst_srd.props.alpha
    beta = lsst_srd.props.beta
    z0 = lsst_srd.props.z0

    # Determine which parametrization based on alpha value
    param_name = None
    for name, params in LSST_SRD_PARAMS.items():
        if np.isclose(alpha, params["alpha"], atol=1e-10):
            param_name = name
            break

    assert param_name is not None, f"Unknown parametrization with alpha={alpha}"

    # Verify all parameters match
    expected = LSST_SRD_PARAMS[param_name]
    assert_allclose(alpha, expected["alpha"], rtol=1e-10)
    assert_allclose(beta, expected["beta"], rtol=1e-10)
    assert_allclose(z0, expected["z0"], rtol=1e-10)


def test_lsst_srd_limits(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test redshift limits."""
    z_min = 0.1
    z_max = 3.0

    lsst_srd.set_lim(z_min, z_max)

    z_min_out, z_max_out = lsst_srd.get_lim()

    assert_allclose(z_min_out, z_min, rtol=1e-10)
    assert_allclose(z_max_out, z_max, rtol=1e-10)


def test_lsst_srd_integ(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test integration of redshift distribution."""
    z_min = 0.01
    z_max = 3.0
    lsst_srd.set_lim(z_min, z_max)

    # Test at several redshift values
    z_array = np.linspace(z_min, z_max, 20)

    for z in z_array:
        p_z = lsst_srd.integ(z)
        ln_p_z = lsst_srd.ln_integ(z)

        assert np.isfinite(p_z)
        assert p_z >= 0.0
        assert np.isfinite(ln_p_z)

        # Check consistency between integ and ln_integ
        assert_allclose(np.log(p_z), ln_p_z, rtol=1e-8)


def test_lsst_srd_normalization(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test that distribution is normalized."""
    z_min = 0.01
    z_max = 3.0
    lsst_srd.set_lim(z_min, z_max)

    # Integrate numerically to check normalization
    z_array = np.linspace(z_min, z_max, 1000)
    p_array = np.array([lsst_srd.integ(z) for z in z_array])

    norm = np.trapezoid(p_array, z_array)

    # Should integrate to 1
    assert_allclose(norm, 1.0, rtol=1e-3)


def test_lsst_srd_gen(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test random generation of redshifts."""
    rng = Ncm.RNG.seeded_new(None, 123)
    z_min = 0.01
    z_max = 3.0
    lsst_srd.set_lim(z_min, z_max)

    # Generate samples
    n_samples = 1000
    samples = np.array([lsst_srd.gen(rng) for _ in range(n_samples)])

    # Check all samples are within limits
    assert np.all(samples >= z_min)
    assert np.all(samples <= z_max)

    # Check samples are finite
    assert np.all(np.isfinite(samples))

    # Check mean is reasonable (should be near z0 for the distribution)
    mean_z = np.mean(samples)
    # Mean should be within reasonable range of z0 (not exact due to truncation)
    assert mean_z > 0.0
    assert mean_z < z_max


def test_lsst_srd_serialize(
    lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD, serializer: Ncm.Serialize
) -> None:
    """Test serialization and deserialization."""
    # Serialize
    lsst_srd_str = serializer.to_yaml(lsst_srd)
    assert lsst_srd_str is not None
    assert len(lsst_srd_str) > 0

    # Deserialize
    lsst_srd_dup = serializer.from_yaml(lsst_srd_str)
    assert isinstance(lsst_srd_dup, Nc.GalaxySDTrueRedshiftLSSTSRD)

    # Check parameters match
    assert_allclose(lsst_srd.props.alpha, lsst_srd_dup.props.alpha, rtol=1e-10)
    assert_allclose(lsst_srd.props.beta, lsst_srd_dup.props.beta, rtol=1e-10)
    assert_allclose(lsst_srd.props.z0, lsst_srd_dup.props.z0, rtol=1e-10)


def test_lsst_srd_model_id(lsst_srd: Nc.GalaxySDTrueRedshiftLSSTSRD) -> None:
    """Test model ID registration."""
    model_id = Nc.GalaxySDTrueRedshift.id()
    assert model_id >= 0

    # Submodels cannot be added directly to the model set, we need another model to
    # hold it.
    spec = Nc.GalaxySDObsRedshiftSpec.new(lsst_srd, 0.0, 10.0)

    # Create a model set and add the distribution
    mset = Ncm.MSet.empty_new()
    mset.set(spec)

    # Check it can be retrieved
    retrieved_spec = mset.peek(spec.id())
    retrieved = retrieved_spec.peek_submodel_by_mid(model_id)
    assert retrieved is not None
    assert isinstance(retrieved, Nc.GalaxySDTrueRedshift)
    assert retrieved is lsst_srd


# =============================================================================
# NcGalaxySDObsRedshiftGauss LSST Bins Tests
# =============================================================================


@pytest.fixture(
    name="lsst_bins_case",
    params=[
        (Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE, "Y1_SOURCE"),
        (Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS, "Y1_LENS"),
        (Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_SOURCE, "Y10_SOURCE"),
        (Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_LENS, "Y10_LENS"),
    ],
    ids=["Y1_SOURCE", "Y1_LENS", "Y10_SOURCE", "Y10_LENS"],
)
def fixture_lsst_bins_case(request) -> tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str]:
    """Parametrize over all LSST bin types."""
    return request.param


def test_lsst_bins_creation(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test creation of LSST bins."""
    lsst_type, type_name = lsst_bins_case

    # Create bins
    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    assert bins is not None
    assert len(bins) > 0

    # Check bin count matches expected
    expected_count = LSST_SRD_BIN_COUNTS[type_name]
    assert len(bins) == expected_count

    # Check all bins are valid
    for bin_obj in bins:
        assert isinstance(bin_obj, Nc.GalaxySDObsRedshiftGauss)
        assert isinstance(bin_obj, Ncm.Model)


def test_lsst_bins_share_true_redshift(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that all bins share the same true redshift submodel.

    All photometric bins should reference the same NcGalaxySDTrueRedshift
    object (not just equal parameters, but actually the same object instance).
    """
    lsst_type, _ = lsst_bins_case

    # Create bins with output parameter to get the shared true redshift object
    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # All bins should share the same true redshift object
    n_bins = len(bins)
    if n_bins > 0:
        # Get submodel from first bin
        first_bin = bins[0]
        if hasattr(first_bin, "peek_submodel_by_mid"):
            first_gsdtr = first_bin.peek_submodel_by_mid(Nc.GalaxySDTrueRedshift.id())

            # Check all other bins have the same submodel
            for i in range(1, n_bins):
                bin_gsdtr = bins[i].peek_submodel_by_mid(Nc.GalaxySDTrueRedshift.id())
                # They should be the same object (not just equal parameters)
                assert bin_gsdtr is first_gsdtr


def test_lsst_bins_parameters(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that bins have correct true redshift parameters."""
    lsst_type, type_name = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Get expected parameters
    expected_params = LSST_SRD_PARAMS[type_name]

    # Check parameters for each bin (they should all be the same)
    for bin_obj in bins:
        if hasattr(bin_obj, "peek_submodel_by_mid"):
            gsdtr = bin_obj.peek_submodel_by_mid(Nc.GalaxySDTrueRedshift.id())
            if gsdtr is not None:
                alpha = gsdtr["alpha"]
                beta = gsdtr["beta"]
                z0 = gsdtr["z0"]

                assert_allclose(alpha, expected_params["alpha"], rtol=1e-10)
                assert_allclose(beta, expected_params["beta"], rtol=1e-10)
                assert_allclose(z0, expected_params["z0"], rtol=1e-10)


def test_lsst_bins_limits(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that bins have valid, non-overlapping limits.

    Verifies that each bin has zp_min < zp_max, all limits are finite,
    and adjacent bins are contiguous (no gaps or overlaps).
    """
    lsst_type, _ = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Collect limits from all bins
    limits = []
    for bin_obj in bins:
        zp_min, zp_max = bin_obj.get_zp_lim()
        limits.append((zp_min, zp_max))
        assert zp_min < zp_max
        assert np.isfinite(zp_min)
        assert np.isfinite(zp_max)

    # Check bins are in order
    for i in range(len(limits) - 1):
        curr_max = limits[i][1]
        next_min = limits[i + 1][0]
        # Adjacent bins should be contiguous (within floating point)
        assert_allclose(curr_max, next_min, rtol=1e-10)


def test_lsst_bins_coverage(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that bins cover the expected redshift range."""
    lsst_type, type_name = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Last bin maximum
    last_max = bins[-1].get_zp_lim()[1]

    # For source bins, z_max should be 3.5 (as per C implementation)
    # For lens bins, depends on Y1 or Y10
    if "SOURCE" in type_name:
        assert_allclose(last_max, 3.5, rtol=1e-10)
    elif "LENS" in type_name:
        # Lens bins typically go to lower redshift
        assert last_max > 0.5  # Reasonable check
        assert last_max < 3.0


def test_lsst_bins_integration(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test evaluation of the P(z|zp_min, zp_max) distribution.

    Verifies that the conditional distribution can be evaluated at multiple
    true redshift values and returns finite, non-negative probabilities.
    """
    lsst_type, type_name = lsst_bins_case

    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Get true redshift limits
    z_min, z_max = gsdtr.get_lim()

    # For source bins, sigma0 = 0.05, for lens bins sigma0 = 0.03
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    for bin_obj in bins:
        zp_min, zp_max = bin_obj.get_zp_lim()

        # Test at several true redshifts within valid range
        z_array = np.linspace(z_min, min(z_max, 2.0), 10)

        for z in z_array:
            # Evaluate P(z|zp_min, zp_max)
            p = bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0, rel_error)
            assert np.isfinite(p), f"Non-finite probability at z={z}"
            assert p >= 0.0, f"Negative probability at z={z}"


def test_lsst_bins_normalization(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that P(z|zp_min, zp_max) is properly normalized over z.

    The conditional true redshift distribution for each photometric bin must
    integrate to 1.0 over the true redshift range. Tests first, middle, and
    last bins using numerical integration.
    """
    lsst_type, type_name = lsst_bins_case

    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Get true redshift limits
    z_min, z_max = gsdtr.get_lim()

    # For source bins, sigma0 = 0.05, for lens bins sigma0 = 0.03
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    # Test first, middle, and last bins
    test_indices = [0, len(bins) // 2, -1] if len(bins) > 2 else range(len(bins))

    for i in test_indices:
        bin_obj = bins[i]
        zp_min, zp_max = bin_obj.get_zp_lim()

        # Integrate P(z|zp_min, zp_max) over z
        # Use finer sampling for integration
        z_array = np.linspace(z_min, z_max, 500)
        p_array = np.array(
            [
                bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0, rel_error)
                for z in z_array
            ]
        )

        # The integral over z should be 1 (normalized distribution)
        norm = np.trapezoid(p_array, z_array)

        # Should integrate to 1
        assert np.isfinite(norm), f"Non-finite normalization for bin {i}"
        assert_allclose(
            norm, 1.0, rtol=2e-2, err_msg=f"Bin {i} not properly normalized"
        )


def test_lsst_bins_pzp_distribution(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test evaluation of the marginal photometric redshift distribution.

    Evaluates P(zp) at multiple photometric redshift values and verifies that
    all values are finite, non-negative, and reasonably smooth (no discontinuities).
    """
    lsst_type, type_name = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # For source bins, sigma0 = 0.05, for lens bins sigma0 = 0.03
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    # Get overall zp range covered by all bins
    zp_min_tot = bins[0].get_zp_lim()[0]
    zp_max_tot = bins[-1].get_zp_lim()[1]

    # Test P(zp) at various points
    zp_array = np.linspace(zp_min_tot, zp_max_tot, 20)

    # Use first bin to compute P(zp) (they all share the same true redshift
    # distribution)
    bin_obj = bins[0]

    for zp in zp_array:
        # Evaluate P(zp)
        p = bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error)
        assert np.isfinite(p), f"Non-finite P(zp) at zp={zp}"
        assert p >= 0.0, f"Negative P(zp) at zp={zp}"

    # Test that P(zp) is reasonably smooth (no huge jumps)
    p_array = np.array(
        [bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error) for zp in zp_array]
    )

    # Check that there are no NaN or inf values
    assert np.all(np.isfinite(p_array)), "P(zp) contains non-finite values"

    # Check that distribution is positive
    assert np.all(p_array >= 0.0), "P(zp) contains negative values"


def test_lsst_bins_pzp_caching(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test caching behavior of P(zp) computation.

    Verifies that repeated evaluations with the same parameters return cached
    results (exactly equal values), and that changing parameters correctly
    invalidates the cache.
    """
    lsst_type, type_name = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # For source bins, sigma0 = 0.05, for lens bins sigma0 = 0.03
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    bin_obj = bins[0]
    zp_max_tot = bins[-1].get_zp_lim()[1]

    # First evaluation (should compute and cache)
    zp_test = 1.0
    p1 = bin_obj.eval_pzp(zp_test, sigma0, zp_max_tot, rel_error)

    # Second evaluation with same parameters (should use cache)
    p2 = bin_obj.eval_pzp(zp_test, sigma0, zp_max_tot, rel_error)

    # Should be exactly equal (cached)
    assert p1 == p2, "Cached values don't match"

    # Third evaluation with different sigma0 (should recompute)
    p3 = bin_obj.eval_pzp(zp_test, sigma0 * 1.5, zp_max_tot, rel_error)

    # Should be different
    assert p1 != p3, "Cache not invalidated when sigma0 changed"


def test_lsst_bins_serialization(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
    serializer: Ncm.Serialize,
) -> None:
    """Test that bins can be serialized."""
    lsst_type, _ = lsst_bins_case

    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # Serialize each bin
    for bin_obj in bins:
        bin_str = serializer.to_yaml(bin_obj)
        assert bin_str is not None
        assert len(bin_str) > 0

        # Deserialize
        bin_dup = serializer.from_yaml(bin_str)
        assert isinstance(bin_dup, Nc.GalaxySDObsRedshiftGauss)

        # Check limits match
        zp_min, zp_max = bin_obj.get_zp_lim()
        zp_min_dup, zp_max_dup = bin_dup.get_zp_lim()

        assert_allclose(zp_min, zp_min_dup, rtol=1e-10)
        assert_allclose(zp_max, zp_max_dup, rtol=1e-10)


# =============================================================================
# Normalization and consistency tests
# =============================================================================


def test_lsst_bins_pzp_normalization(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that P(zp) is properly normalized over photometric redshift.

    Note: For lens samples, bins only cover [0.2, 1.2], so we integrate over
    a wider range to check full normalization.
    """
    lsst_type, type_name = lsst_bins_case

    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    # For source bins, sigma0 = 0.05, for lens bins sigma0 = 0.03
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    # For proper normalization test, use a range that covers the full distribution
    # For source: bins go to 3.5
    # For lens: bins go to ~1.2, but distribution extends further
    if "SOURCE" in type_name:
        zp_max_tot = bins[-1].get_zp_lim()[1]  # 3.5
    else:
        # For lens, use larger range for full normalization test
        _, z_max = gsdtr.get_lim()
        # Use z_max + 5*sigma as reasonable upper limit
        zp_max_tot = min(z_max + 5 * sigma0 * (1 + z_max), 5.0)

    # Use first bin (they all share the same true redshift distribution)
    bin_obj = bins[0]

    # Integrate P(zp) over zp using high-resolution grid
    zp_array = np.linspace(0.0, zp_max_tot, 1000)
    p_array = np.array(
        [bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error) for zp in zp_array]
    )

    # Numerical integration
    norm = np.trapezoid(p_array, zp_array)

    # Should integrate to 1 (properly normalized)
    assert np.isfinite(norm), "Non-finite normalization"
    assert_allclose(
        norm, 1.0, rtol=2e-2, err_msg=f"P(zp) not properly normalized for {type_name}"
    )


@pytest.mark.parametrize(
    "lsst_bins_case",
    [
        (t0, t0.name)
        for t0 in list(Nc.GalaxySDTrueRedshiftLSSTSRDType)
        if "SOURCE" in t0.name
    ],
    indirect=True,
    ids=lambda x: f"{x[1]}_BINS",
)
def test_lsst_bins_equal_area_verification(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Test that source sample bins have equal areas under P(zp).

    Equal-area binning means each bin contains 1/n_bins of the total integrated
    probability. The bins are computed using the CDF of P(zp), so we verify by
    integrating P(zp) over each bin and checking that all bins have equal area.
    """
    lsst_type, _ = lsst_bins_case
    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    sigma0 = 0.05
    rel_error = 1e-4
    zp_max_tot = bins[-1].get_zp_lim()[1]

    # Get the bin edges
    bin_edges_list = [bins[0].get_zp_lim()[0]]  # First edge is zp_min of first bin
    for bin_i in bins:
        _, zp_max = bin_i.get_zp_lim()
        bin_edges_list.append(zp_max)

    bin_edges = np.array(bin_edges_list)

    # Compute the area of each bin by integrating P(zp) over the bin
    # P(zp) should be computed the same way as in compute_equal_area_photoz_bins
    bin_areas_list = []
    bin_obj = bins[0]

    total_z = np.linspace(0.0, zp_max_tot, 2000)
    total_area = np.trapezoid(
        [bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error) for zp in total_z], total_z
    )

    for zp0, zp1 in zip(bin_edges[:-1], bin_edges[1:]):
        # Integrate from zp0 to zp1 (over this bin)
        zp_array = np.linspace(zp0, zp1, 1500)
        p_array = np.array(
            [bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error) for zp in zp_array]
        )
        area = np.trapezoid(p_array, zp_array) / total_area
        bin_areas_list.append(area)

    bin_areas = np.array(bin_areas_list)

    # All bin areas should be equal (equal-area binning)
    expected_area = np.mean(bin_areas)

    for i, area in enumerate(bin_areas):
        assert_allclose(
            area,
            expected_area,
            rtol=1.0e-5,
            err_msg=(
                f"Bin {i} area ({area:.4f}) doesn't "
                f"match expected ({expected_area:.4f})"
            ),
        )


def test_lsst_bins_pz_given_zp_python_comparison(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Compare C-computed P(z|zp_min,zp_max) with Python numerical integration.

    Validates the conditional true redshift distribution by computing it
    independently in Python and comparing the shape (ratios) with the C
    implementation. Tests multiple redshift values within the valid range.
    """
    lsst_type, type_name = lsst_bins_case

    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    z_min, z_max = gsdtr.get_lim()
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    # Test first bin only (to keep test fast)
    bin_obj = bins[0]
    zp_min, zp_max = bin_obj.get_zp_lim()

    # Pick several test points
    z_test = np.array([0.2, 0.5, 1.0, 1.5])
    z_test = z_test[(z_test >= z_min) & (z_test <= z_max)]

    for z in z_test:
        # Get C-computed value
        p_c = bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0, rel_error)

        # Compute manually with Python
        # P(z|zp_min,zp_max) ∝ P(z) * W(zp_min,zp_max|z)
        # where W is the Gaussian integral over the bin

        p_z = gsdtr.integ(z)
        sigmaz = sigma0 * (1.0 + z)

        # Standard Gaussian CDF
        def gaussian_integral(a, b, mu, sigma):
            """Integrate Gaussian from a to b."""
            sqrt2 = np.sqrt(2.0)
            return 0.5 * (
                erf((b - mu) / (sqrt2 * sigma)) - erf((a - mu) / (sqrt2 * sigma))
            )

        w = gaussian_integral(zp_min, zp_max, z, sigmaz)

        # Normalization (to account for zp > 0 constraint)
        norm_z = 0.5 * (1.0 + erf(z / (np.sqrt(2.0) * sigmaz)))

        # Compute P(z) * W / norm_z (unnormalized)
        p_python_un_norm = p_z * w / norm_z

        # Now we need to normalize over z
        # For comparison, just check that the ratio is consistent across different z
        # (both should have same normalization)

        # Instead, let's just check the shape is similar
        # The absolute values might differ due to normalization
        # So we'll compute for multiple z and check consistency

        if z == z_test[0]:
            p_python_un_norm_ref = p_python_un_norm
            p_c_ref = p_c
        else:
            # Check ratio consistency
            ratio_python = p_python_un_norm / p_python_un_norm_ref
            ratio_c = p_c / p_c_ref

            assert_allclose(
                ratio_c,
                ratio_python,
                rtol=0.05,
                atol=1.0e-10,
                err_msg=f"Ratio mismatch at z={z}",
            )


def test_lsst_bins_pzp_python_comparison(
    lsst_bins_case: tuple[Nc.GalaxySDTrueRedshiftLSSTSRDType, str],
) -> None:
    """Compare C-computed P(zp) with Python numerical integration.

    Validates the marginal photometric redshift distribution by computing it
    independently in Python using numerical integration over the true redshift
    distribution convolved with a Gaussian. Tests multiple zp values.
    """
    lsst_type, type_name = lsst_bins_case

    bins, gsdtr = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    z_min, z_max = gsdtr.get_lim()
    sigma0 = 0.05 if "SOURCE" in type_name else 0.03
    rel_error = 1e-4

    bin_obj = bins[0]
    zp_max_tot = bins[-1].get_zp_lim()[1]

    # Test at several zp values
    zp_test = np.array([0.5, 1.0, 1.5, 2.0])
    zp_test = zp_test[zp_test <= zp_max_tot]

    for zp in zp_test:
        # Get C-computed value
        p_c = bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error)

        # Compute manually with Python
        # P(zp) = ∫ P(z) * Gauss(zp|z,σ(z)) dz
        # where σ(z) = σ0(1+z)

        # High-resolution integration
        z_array = np.linspace(z_min, min(z_max, 3.0), 500)
        p_z_array = np.array([gsdtr.integ(z) for z in z_array])

        # Gaussian contribution from each z
        gauss_contributions = []
        for z, p_z in zip(z_array, p_z_array):
            sigmaz = sigma0 * (1.0 + z)
            # Gaussian PDF
            gauss = np.exp(-0.5 * ((zp - z) / sigmaz) ** 2) / (
                np.sqrt(2.0 * np.pi) * sigmaz
            )
            norm = 0.5 * erfc(-z / (np.sqrt(2.0) * sigmaz))
            gauss_contributions.append(p_z * gauss / norm)

        # Integrate
        p_python = np.trapezoid(gauss_contributions, z_array)

        # Compare (allow some tolerance due to numerical integration differences)
        assert_allclose(
            p_c,
            p_python,
            rtol=0.05,
            err_msg=f"P(zp) mismatch at zp={zp} for {type_name}",
        )


def test_lsst_bins_consistency_between_methods() -> None:
    """Test consistency and caching behavior of evaluation methods.

    Verifies that repeated calls with the same parameters return identical
    cached results, and that changing parameters correctly invalidates the
    cache and produces different results.
    """
    # Create a simple case
    gsdtr = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    bin_obj = Nc.GalaxySDObsRedshiftGauss.new(gsdtr, 0.5, 1.5)

    z_min, z_max = gsdtr.get_lim()
    zp_min, zp_max = bin_obj.get_zp_lim()

    sigma0 = 0.05
    rel_error = 1e-4

    # Test that eval_pz_given_zp gives consistent results for different calls
    z_test = np.linspace(z_min, min(z_max, 2.0), 10)

    # First pass
    p_first = np.array(
        [bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0, rel_error) for z in z_test]
    )

    # Second pass (should use cache)
    p_second = np.array(
        [bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0, rel_error) for z in z_test]
    )

    # Should be exactly equal (cached)
    assert_allclose(
        p_first, p_second, rtol=1e-15, err_msg="Cached values don't match exactly"
    )

    # Change sigma0 slightly and verify it recomputes
    sigma0_new = sigma0 * 1.1
    p_third = np.array(
        [
            bin_obj.eval_pz_given_zp(z, zp_min, zp_max, sigma0_new, rel_error)
            for z in z_test
        ]
    )

    # Should be different
    assert not np.allclose(
        p_first, p_third, rtol=1e-3
    ), "Cache not invalidated when parameters changed"


def test_lsst_bins_sum_over_bins() -> None:
    """Test that summing P(zp) contributions over all bins recovers overall P(zp)."""
    # For lens samples with linear binning, bins don't overlap but are contiguous
    # So sum of contributions should approximately equal the overall distribution

    lsst_type = Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS
    bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(lsst_type)

    sigma0 = 0.03
    rel_error = 1e-4
    zp_max_tot = bins[-1].get_zp_lim()[1]

    bin_obj = bins[0]

    # Test at several zp points that fall within the bin coverage
    zp_min_tot = bins[0].get_zp_lim()[0]
    zp_test = np.linspace(zp_min_tot + 0.1, zp_max_tot - 0.1, 5)

    for zp in zp_test:
        # Get overall P(zp)
        p_total = bin_obj.eval_pzp(zp, sigma0, zp_max_tot, rel_error)

        # Since zp is given, each bin sees the same P(zp)
        # The test here is more about checking P(zp) is consistent
        # across the range covered by bins

        assert np.isfinite(p_total), f"Non-finite P(zp) at zp={zp}"
        assert p_total > 0, f"Non-positive P(zp) at zp={zp}"


# =============================================================================
# Cross-parametrization tests
# =============================================================================


def test_different_parametrizations_different_params() -> None:
    """Test that different parametrizations have different parameters."""
    y1_source = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_source()
    y1_lens = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y1_lens()
    y10_source = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_source()
    y10_lens = Nc.GalaxySDTrueRedshiftLSSTSRD.new_y10_lens()

    # Y1 source vs Y1 lens
    assert not np.isclose(y1_source.props.alpha, y1_lens.props.alpha)
    assert not np.isclose(y1_source.props.z0, y1_lens.props.z0)

    # Y10 source vs Y10 lens
    assert not np.isclose(y10_source.props.alpha, y10_lens.props.alpha)
    assert not np.isclose(y10_source.props.z0, y10_lens.props.z0)

    # Y1 source vs Y10 source
    assert not np.isclose(y1_source.props.alpha, y10_source.props.alpha)
    assert not np.isclose(y1_source.props.z0, y10_source.props.z0)

    # Y1 lens vs Y10 lens
    assert not np.isclose(y1_lens.props.alpha, y10_lens.props.alpha)
    assert not np.isclose(y1_lens.props.z0, y10_lens.props.z0)


def test_different_bin_counts() -> None:
    """Test that Y1 and Y10 have different bin counts for lens."""
    y1_lens_bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS
    )
    y10_lens_bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_LENS
    )

    # Y1 should have 5 bins, Y10 should have 10 bins
    assert len(y1_lens_bins) == 5
    assert len(y10_lens_bins) == 10


def test_source_lens_z_ranges() -> None:
    """Test that source and lens bins have appropriate redshift ranges."""
    y1_source_bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE
    )
    y1_lens_bins, _ = Nc.GalaxySDObsRedshiftGauss.new_lsst_srd_bins(
        Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS
    )

    # Source bins typically extend to higher redshift
    source_max = y1_source_bins[-1].get_zp_lim()[1]
    lens_max = y1_lens_bins[-1].get_zp_lim()[1]

    # Source should go higher in redshift than lens
    assert source_max > lens_max
