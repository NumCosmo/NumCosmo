#!/usr/bin/env python
#
# test_py_cluster_mass_ascaso.py
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

"""Tests for NcClusterMassAscaso model."""

import math
import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

LN_RICHNESS_CUT = np.log(5)


@pytest.fixture(name="cluster_m")
def fixture_cluster_m() -> Nc.ClusterMassAscaso:
    """Create cluster mass-richness relation model."""
    cluster_m = Nc.ClusterMassAscaso(
        lnRichness_min=LN_RICHNESS_CUT, lnRichness_max=np.log(200)
    )
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473)
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048)
    cluster_m.param_set_by_name("cut", LN_RICHNESS_CUT)
    return cluster_m


def test_get_cut(cluster_m: Nc.ClusterMassAscaso) -> None:
    """Test richness cut parameter."""
    assert cluster_m.get_cut() == LN_RICHNESS_CUT


def test_cluster_mass_ascaso_mean_std(cluster_m: Nc.ClusterMassAscaso) -> None:
    """Test mean and standard deviation of truncated distribution."""
    nsize = 100
    lnM = np.linspace(np.log(1e13), np.log(1e16), nsize)
    z = np.linspace(0, 1.1, nsize)

    for i in range(nsize):
        assert cluster_m.get_mean(lnM[i], 0.5) >= cluster_m.mu(lnM[i], 0.5)
        assert cluster_m.get_mean(np.log(1e14), z[i]) >= cluster_m.mu(
            np.log(1e14), z[i]
        )
        assert cluster_m.get_std(lnM[i], 0.5) <= cluster_m.sigma(lnM[i], 0.5)
        assert cluster_m.get_std(np.log(1e14), z[i]) <= cluster_m.sigma(
            np.log(1e14), z[i]
        )


def test_serialization_deserialization(cluster_m: Nc.ClusterMassAscaso) -> None:
    """Test serialization and deserialization."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    cluster_m_dup = ser.dup_obj(cluster_m)

    assert cluster_m_dup is not None
    assert isinstance(cluster_m_dup, Nc.ClusterMassAscaso)
    assert cluster_m_dup is not cluster_m


def test_cluster_mass_ascaso_p_basic(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution function."""
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.mu(lnM, z)

    p_mean = cluster_m.p(cosmo, lnM, z, [mean_lnR], None)
    assert p_mean > 0.0

    p_above = cluster_m.p(cosmo, lnM, z, [mean_lnR + 1.0], None)
    p_below = cluster_m.p(cosmo, lnM, z, [mean_lnR - 1.0], None)
    assert p_above < p_mean
    assert p_below < p_mean

    p_far = cluster_m.p(cosmo, lnM, z, [mean_lnR + 5.0], None)
    assert p_far < p_mean * 0.01


def test_cluster_mass_ascaso_p_integral(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test probability distribution integrates correctly."""
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.mu(lnM, z)
    std_lnR = cluster_m.sigma(lnM, z)

    lnR_array = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 200)
    p_array = np.array([cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array])

    integral = np.trapezoid(p_array, lnR_array)
    intp_val = cluster_m.intp(cosmo, lnM, z)

    assert np.isclose(integral, intp_val, rtol=1e-3)


def test_cluster_mass_ascaso_intp_basic(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test cumulative distribution function."""
    lnM = np.log(1e14)
    z = 0.5

    intp_val = cluster_m.intp(cosmo, lnM, z)
    assert 0.0 <= intp_val <= 1.0


def test_cluster_mass_ascaso_intp_bin_basic(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test binned integration of probability distribution."""
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.mu(lnM, z)
    std_lnR = cluster_m.sigma(lnM, z)

    bin_val = cluster_m.intp_bin(cosmo, lnM, z, [mean_lnR], [mean_lnR + std_lnR], None)
    assert bin_val > 0.0

    lnR_edges = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 11)
    bin_sum = sum(
        cluster_m.intp_bin(cosmo, lnM, z, [lnR_edges[i]], [lnR_edges[i + 1]], None)
        for i in range(len(lnR_edges) - 1)
    )
    total_integral = cluster_m.intp(cosmo, lnM, z)

    assert np.isclose(bin_sum, total_integral, rtol=1e-3)


def test_cluster_mass_ascaso_intp_bin_consistency(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test intp_bin() matches manual integration of p()."""
    lnM = np.log(1e14)
    z = 0.5
    mean_lnR = cluster_m.mu(lnM, z)
    std_lnR = cluster_m.sigma(lnM, z)

    lnR_lower = mean_lnR - std_lnR
    lnR_upper = mean_lnR + std_lnR

    lnR_array = np.linspace(lnR_lower, lnR_upper, 200)
    p_array = np.array([cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array])
    manual_integral = np.trapezoid(p_array, lnR_array)

    intp_bin_val = cluster_m.intp_bin(cosmo, lnM, z, [lnR_lower], [lnR_upper], None)

    assert np.isclose(manual_integral, intp_bin_val, rtol=1e-3)


def test_cluster_mass_ascaso_truncated_gaussian(
    cluster_m: Nc.ClusterMassAscaso,
) -> None:
    """Test truncated Gaussian mean and std corrections."""
    lnM = np.log(1e14)
    z = 0.5
    cut = cluster_m.get_cut()

    mean_untrunc = cluster_m.mu(lnM, z)
    std_untrunc = cluster_m.sigma(lnM, z)

    mean_trunc = cluster_m.get_mean(lnM, z)
    std_trunc = cluster_m.get_std(lnM, z)

    assert mean_trunc > mean_untrunc
    assert std_trunc < std_untrunc

    alpha = (cut - mean_untrunc) / std_untrunc
    phi_alpha = np.exp(-0.5 * alpha**2) / np.sqrt(2 * np.pi)
    Phi_alpha = 0.5 * (1 + math.erf(alpha / np.sqrt(2)))

    expected_mean = mean_untrunc + std_untrunc * phi_alpha / (1 - Phi_alpha)
    assert np.isclose(mean_trunc, expected_mean, rtol=1e-10)

    lambda_alpha = phi_alpha / (1 - Phi_alpha)
    expected_std = std_untrunc * np.sqrt(1 + alpha * lambda_alpha - lambda_alpha**2)
    assert np.isclose(std_trunc, expected_std, rtol=1e-10)


def test_cluster_mass_ascaso_p_vec_z_lnMobs(
    cluster_m: Nc.ClusterMassAscaso, cosmo: Nc.HICosmo
) -> None:
    """Test vectorized probability distribution function.

    Verifies that p_vec_z_lnMobs() matches individual p() calls
    for multiple redshifts and richness values.
    """
    lnM = np.log(1e14)
    z_array = np.linspace(0.1, 1.1, 100)
    lnR_array = np.log(np.geomspace(3.0, 1e3, 100))

    # Create Ncm.Vector for z and Ncm.Matrix for lnR (one column)
    z_vec = Ncm.Vector.new_array(z_array)
    lnR_mat = Ncm.Matrix.new_array(lnR_array, 1)

    # Allocate result array
    res = Ncm.Vector.new(len(z_array))

    # Vectorized call
    cluster_m.p_vec_z_lnMobs(cosmo, lnM, z_vec, lnR_mat, None, res)

    # Compare with individual calls
    for i, (z, lnR) in enumerate(zip(z_array, lnR_array)):
        p_scalar = cluster_m.p(cosmo, lnM, z, [lnR], None)
        p_vec = res.get(i)
        assert np.isclose(p_scalar, p_vec, rtol=1e-10)
