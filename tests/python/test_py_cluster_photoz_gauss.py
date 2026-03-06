#!/usr/bin/env python
#
# test_py_cluster_photoz_gauss.py
#
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
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

"""Tests for NcClusterPhotozGauss model."""

import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()


@pytest.fixture(name="cluster_pz")
def fixture_cluster_pz() -> Nc.ClusterPhotozGauss:
    """Create cluster photometric redshift model."""
    return Nc.ClusterPhotozGauss(pz_min=0.1, pz_max=0.9)


def test_serialization_deserialization(cluster_pz: Nc.ClusterPhotozGauss) -> None:
    """Test serialization and deserialization."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    cluster_pz_dup = ser.dup_obj(cluster_pz)

    assert cluster_pz_dup is not None
    assert isinstance(cluster_pz_dup, Nc.ClusterPhotozGauss)
    assert cluster_pz_dup is not cluster_pz


def test_cluster_photoz_gauss_p_basic(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution function."""
    lnM = np.log(1e14)
    z = 0.5
    z_obs = [0.52]
    z_obs_params = [0.02, 0.03]  # bias, sigma

    p_val = cluster_pz.p(cosmo, lnM, z, z_obs, z_obs_params)
    assert p_val > 0.0

    # Test with different observed redshifts
    z_obs_close = [0.51]
    z_obs_far = [0.8]

    p_close = cluster_pz.p(cosmo, lnM, z, z_obs_close, z_obs_params)
    p_far = cluster_pz.p(cosmo, lnM, z, z_obs_far, z_obs_params)

    assert p_close > p_far


def test_cluster_photoz_gauss_intp_basic(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test integrated probability function."""
    lnM = np.log(1e14)
    z = 0.5

    intp_val = cluster_pz.intp(cosmo, lnM, z)
    assert 0.0 <= intp_val <= 1.0


def test_cluster_photoz_gauss_intp_bin_basic(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test binned integration of probability distribution."""
    lnM = np.log(1e14)
    z = 0.5
    z_obs_params = [0.02, 0.03]  # bias, sigma

    z_obs_lower = [0.4]
    z_obs_upper = [0.6]

    bin_val = cluster_pz.intp_bin(cosmo, lnM, z, z_obs_lower, z_obs_upper, z_obs_params)
    assert bin_val > 0.0
    assert bin_val <= 1.0


def test_cluster_photoz_gauss_resample(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test resampling function."""
    lnM = np.log(1e14)
    z = 0.5
    z_obs_params = [0.02, 0.03]  # bias, sigma
    rng = Ncm.RNG.new()

    success, z_obs = cluster_pz.resample(cosmo, lnM, z, z_obs_params, rng)

    assert isinstance(success, bool)
    assert z_obs > 0.0  # Should be positive


def test_cluster_photoz_gauss_p_limits(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test probability limits function."""
    z_obs = [0.5]
    z_obs_params = [0.02, 0.03]  # bias, sigma

    z_lower, z_upper = cluster_pz.p_limits(cosmo, z_obs, z_obs_params)

    assert z_lower >= 0.0
    assert z_upper > z_lower
    assert z_upper - z_lower > 10 * z_obs_params[1]  # Should span multiple sigmas


def test_cluster_photoz_gauss_p_bin_limits(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test probability bin limits function."""
    z_obs_lower = [0.4]
    z_obs_upper = [0.6]
    z_obs_params = [0.02, 0.03]  # bias, sigma

    z_lower, z_upper = cluster_pz.p_bin_limits(
        cosmo, z_obs_lower, z_obs_upper, z_obs_params
    )

    assert z_lower >= 0.0
    assert z_upper > z_lower


def test_cluster_photoz_gauss_n_limits(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test normalization limits function."""
    z_lower, z_upper = cluster_pz.n_limits(cosmo)

    assert z_lower >= 0.0
    assert z_upper > z_lower


def test_cluster_photoz_gauss_volume(cluster_pz: Nc.ClusterPhotozGauss) -> None:
    """Test volume calculation."""
    volume = cluster_pz.volume()
    assert volume >= 0.0


def test_cluster_photoz_gauss_gaussian_properties(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test Gaussian distribution properties."""
    lnM = np.log(1e14)
    z = 0.5
    bias = 0.02
    sigma = 0.03
    z_obs_params = [bias, sigma]

    # Test that probability is maximum at the mean
    z_mean = z + bias
    z_obs_mean = [z_mean]
    z_obs_offset = [z_mean + sigma]

    p_mean = cluster_pz.p(cosmo, lnM, z, z_obs_mean, z_obs_params)
    p_offset = cluster_pz.p(cosmo, lnM, z, z_obs_offset, z_obs_params)

    assert p_mean > p_offset


def test_cluster_photoz_gauss_bias_effect(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test effect of bias parameter."""
    lnM = np.log(1e14)
    z = 0.5
    sigma = 0.03
    z_obs = [0.52]

    # Test with different bias values
    z_obs_params_no_bias = [0.0, sigma]
    z_obs_params_with_bias = [0.02, sigma]

    p_no_bias = cluster_pz.p(cosmo, lnM, z, z_obs, z_obs_params_no_bias)
    p_with_bias = cluster_pz.p(cosmo, lnM, z, z_obs, z_obs_params_with_bias)

    # With positive bias, the distribution shifts, so probabilities should differ
    assert p_no_bias != p_with_bias


def test_cluster_photoz_gauss_sigma_effect(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test effect of sigma parameter."""
    lnM = np.log(1e14)
    z = 0.5
    bias = 0.02
    z_obs = [0.52]

    # Test with different sigma values
    z_obs_params_small_sigma = [bias, 0.01]
    z_obs_params_large_sigma = [bias, 0.05]

    p_small_sigma = cluster_pz.p(cosmo, lnM, z, z_obs, z_obs_params_small_sigma)
    p_large_sigma = cluster_pz.p(cosmo, lnM, z, z_obs, z_obs_params_large_sigma)

    # Smaller sigma should give higher peak probability
    assert p_small_sigma > p_large_sigma


def test_cluster_photoz_gauss_p_integral(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution integrates correctly using trapezoid."""
    lnM = np.log(1e14)
    z = 0.5
    bias = 0.02
    sigma = 0.03
    z_obs_params = [bias, sigma]

    # Create array of observed redshifts around the mean
    z_mean = z + bias
    z_obs_array = np.linspace(z_mean - 5 * sigma, z_mean + 5 * sigma, 200)
    p_array = np.array(
        [cluster_pz.p(cosmo, lnM, z, [z_obs], z_obs_params) for z_obs in z_obs_array]
    )
    assert all(p > 0 for p in p_array)

    integral = np.trapezoid(p_array, z_obs_array)
    intp_val = cluster_pz.intp(cosmo, lnM, z)

    assert np.isclose(integral, intp_val, rtol=1e-2)


def test_cluster_photoz_gauss_intp_bin_consistency(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test intp_bin matches manual integration of p using trapezoid."""
    lnM = np.log(1e14)
    z = 0.5
    bias = 0.02
    sigma = 0.03
    z_obs_params = [bias, sigma]

    z_obs_lower = [0.4]
    z_obs_upper = [0.6]

    # Manual integration using trapezoid
    z_obs_array = np.linspace(z_obs_lower[0], z_obs_upper[0], 200)
    p_array = np.array(
        [cluster_pz.p(cosmo, lnM, z, [z_obs], z_obs_params) for z_obs in z_obs_array]
    )
    manual_integral = np.trapezoid(p_array, z_obs_array)

    # Compare with intp_bin
    intp_bin_val = cluster_pz.intp_bin(
        cosmo, lnM, z, z_obs_lower, z_obs_upper, z_obs_params
    )

    assert np.isclose(manual_integral, intp_bin_val, rtol=1e-2)


def test_cluster_photoz_gauss_bin_sum_consistency(
    cluster_pz: Nc.ClusterPhotozGauss, cosmo: Nc.HICosmo
) -> None:
    """Test sum of bins matches total integral using trapezoid."""
    lnM = np.log(1e14)
    z = 0.5
    bias = 0.02
    sigma = 0.03
    z_obs_params = [bias, sigma]

    # Create bins spanning the distribution
    z_mean = z + bias
    z_edges = np.linspace(z_mean - 5 * sigma, z_mean + 5 * sigma, 11)

    bin_sum = sum(
        cluster_pz.intp_bin(cosmo, lnM, z, [z_edges[i]], [z_edges[i + 1]], z_obs_params)
        for i in range(len(z_edges) - 1)
    )

    total_integral = cluster_pz.intp(cosmo, lnM, z)

    assert np.isclose(bin_sum, total_integral, rtol=1e-2)
