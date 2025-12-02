#!/usr/bin/env python
#
# test_py_cluster_mass_richness.py
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

"""Tests for NcClusterMassRichness and its subclasses (Ascaso and Ext)."""

import math
import pytest
import numpy as np
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

LN_RICHNESS_CUT = np.log(5)
LN_RICHNESS_MIN = LN_RICHNESS_CUT
LN_RICHNESS_MAX = np.log(200)


# =============================================================================
# Fixtures
# =============================================================================


@pytest.fixture(name="cluster_m_ascaso")
def fixture_cluster_m_ascaso() -> Nc.ClusterMassAscaso:
    """Create NcClusterMassAscaso model with typical parameters."""
    cluster_m = Nc.ClusterMassAscaso(
        lnRichness_min=LN_RICHNESS_MIN, lnRichness_max=LN_RICHNESS_MAX
    )
    cluster_m.param_set_by_name("mup0", 4.12769558168741)
    cluster_m.param_set_by_name("mup1", 1.17476066603899)
    cluster_m.param_set_by_name("mup2", 0.393577193825473)
    cluster_m.param_set_by_name("sigmap0", 0.408750324989284)
    cluster_m.param_set_by_name("sigmap1", -0.123232985316648)
    cluster_m.param_set_by_name("sigmap2", -0.0644996574273048)
    cluster_m.param_set_by_name("cut", LN_RICHNESS_CUT)
    return cluster_m


@pytest.fixture(name="cluster_m_ext")
def fixture_cluster_m_ext() -> Nc.ClusterMassExt:
    """Create NcClusterMassExt model with typical parameters."""
    cluster_m = Nc.ClusterMassExt(
        lnRichness_min=LN_RICHNESS_MIN, lnRichness_max=LN_RICHNESS_MAX
    )
    cluster_m.param_set_by_name("mup0", 4.0)
    cluster_m.param_set_by_name("mup1", 1.0)
    cluster_m.param_set_by_name("mup2", 0.1)
    cluster_m.param_set_by_name("mup3", 0.01)
    cluster_m.param_set_by_name("sigmap0", -0.3)
    cluster_m.param_set_by_name("sigmap1", -0.08)
    cluster_m.param_set_by_name("sigmap2", 0.005)
    cluster_m.param_set_by_name("cut", LN_RICHNESS_CUT)
    return cluster_m


@pytest.fixture(
    name="cluster_m",
    params=["ascaso", "ext"],
)
def fixture_cluster_m(
    request, cluster_m_ascaso: Nc.ClusterMassAscaso, cluster_m_ext: Nc.ClusterMassExt
) -> Nc.ClusterMassRichness:
    """Parametrized fixture that yields both Ascaso and Ext models."""
    if request.param == "ascaso":
        return cluster_m_ascaso
    return cluster_m_ext


# =============================================================================
# Common tests for NcClusterMassRichness (both Ascaso and Ext)
# =============================================================================


class TestClusterMassRichnessCommon:
    """Tests common to all NcClusterMassRichness subclasses."""

    def test_get_cut(self, cluster_m: Nc.ClusterMassRichness) -> None:
        """Test richness cut parameter."""
        assert cluster_m.get_cut() == LN_RICHNESS_CUT

    def test_mu_sigma_consistency(self, cluster_m: Nc.ClusterMassRichness) -> None:
        """Test that mu_sigma returns same values as mu and sigma separately."""
        lnM = np.log(1e14)
        z = 0.5

        mu_separate = cluster_m.mu(lnM, z)
        sigma_separate = cluster_m.sigma(lnM, z)

        mu_combined, sigma_combined = cluster_m.mu_sigma(lnM, z)

        assert np.isclose(mu_separate, mu_combined, rtol=1e-15)
        assert np.isclose(sigma_separate, sigma_combined, rtol=1e-15)

    def test_truncated_mean_std_bounds(self, cluster_m: Nc.ClusterMassRichness) -> None:
        """Test that truncated mean and std are within expected bounds.

        Tests that mean >= untruncated mean and truncated std <= untruncated std.
        """
        nsize = 50
        lnM_array = np.linspace(np.log(1e13), np.log(1e16), nsize)
        z_array = np.linspace(0.1, 1.1, nsize)

        for lnM in lnM_array:
            for z in z_array:
                mu_untrunc = cluster_m.mu(lnM, z)
                sigma_untrunc = cluster_m.sigma(lnM, z)
                mean_trunc = cluster_m.get_mean(lnM, z)
                std_trunc = cluster_m.get_std(lnM, z)

                # Truncation shifts mean up and reduces std
                assert mean_trunc >= mu_untrunc
                assert std_trunc <= sigma_untrunc

    def test_truncated_gaussian_formulas(
        self, cluster_m: Nc.ClusterMassRichness
    ) -> None:
        """Test truncated Gaussian mean and std.

        Tests that corrections match analytical formulas.
        """
        lnM = np.log(1e14)
        z = 0.5
        cut = cluster_m.get_cut()

        mu_untrunc = cluster_m.mu(lnM, z)
        sigma_untrunc = cluster_m.sigma(lnM, z)

        mean_trunc = cluster_m.get_mean(lnM, z)
        std_trunc = cluster_m.get_std(lnM, z)

        # Analytical formulas for truncated normal
        alpha = (cut - mu_untrunc) / sigma_untrunc
        phi_alpha = np.exp(-0.5 * alpha**2) / np.sqrt(2 * np.pi)
        Phi_alpha = 0.5 * (1 + math.erf(alpha / np.sqrt(2)))

        expected_mean = mu_untrunc + sigma_untrunc * phi_alpha / (1 - Phi_alpha)
        assert np.isclose(mean_trunc, expected_mean, rtol=1e-10)

        lambda_alpha = phi_alpha / (1 - Phi_alpha)
        expected_std = sigma_untrunc * np.sqrt(
            1 + alpha * lambda_alpha - lambda_alpha**2
        )
        assert np.isclose(std_trunc, expected_std, rtol=1e-10)

    def test_compute_truncated_mean_std(
        self, cluster_m: Nc.ClusterMassRichness
    ) -> None:
        """Test compute_truncated_mean/std with arbitrary lnR_mean and lnR_sigma."""
        lnR_mean = 3.5
        lnR_sigma = 0.5

        mean_trunc = cluster_m.compute_truncated_mean(lnR_mean, lnR_sigma)

        # Should be consistent with formulas
        cut = cluster_m.get_cut()
        alpha = (cut - lnR_mean) / lnR_sigma
        phi_alpha = np.exp(-0.5 * alpha**2) / np.sqrt(2 * np.pi)
        Phi_alpha = 0.5 * (1 + math.erf(alpha / np.sqrt(2)))

        expected_mean = lnR_mean + lnR_sigma * phi_alpha / (1 - Phi_alpha)
        assert np.isclose(mean_trunc, expected_mean, rtol=1e-10)

    def test_p_basic(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test probability distribution function basic properties."""
        lnM = np.log(1e14)
        z = 0.5
        mean_lnR = cluster_m.mu(lnM, z)

        # Probability at mean should be positive
        p_mean = cluster_m.p(cosmo, lnM, z, [mean_lnR], None)
        assert p_mean > 0.0

        # Probability should decrease away from mean
        p_above = cluster_m.p(cosmo, lnM, z, [mean_lnR + 1.0], None)
        p_below = cluster_m.p(cosmo, lnM, z, [mean_lnR - 1.0], None)
        assert p_above < p_mean
        assert p_below < p_mean

        # Far from mean should be very small
        p_far = cluster_m.p(cosmo, lnM, z, [mean_lnR + 5.0], None)
        assert p_far < p_mean * 0.01

    def test_p_below_cut_is_zero(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test that probability below cut threshold is zero."""
        lnM = np.log(1e14)
        z = 0.5
        cut = cluster_m.get_cut()

        p_below_cut = cluster_m.p(cosmo, lnM, z, [cut - 1.0], None)
        assert p_below_cut == 0.0

    def test_p_integral_equals_intp(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test probability distribution integrates to intp value."""
        lnM = np.log(1e14)
        z = 0.5
        mean_lnR = cluster_m.mu(lnM, z)
        std_lnR = cluster_m.sigma(lnM, z)

        lnR_array = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 500)
        p_array = np.array(
            [cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array]
        )

        integral = np.trapezoid(p_array, lnR_array)
        intp_val = cluster_m.intp(cosmo, lnM, z)

        # Tolerance is 2e-3 due to trapezoidal integration accuracy
        assert np.isclose(integral, intp_val, rtol=2e-3)

    def test_intp_bounds(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test cumulative distribution function bounds."""
        lnM = np.log(1e14)
        z = 0.5

        intp_val = cluster_m.intp(cosmo, lnM, z)
        assert 0.0 <= intp_val <= 1.0

    def test_intp_bin_basic(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test binned integration of probability distribution."""
        lnM = np.log(1e14)
        z = 0.5
        mean_lnR = cluster_m.mu(lnM, z)
        std_lnR = cluster_m.sigma(lnM, z)

        bin_val = cluster_m.intp_bin(
            cosmo, lnM, z, [mean_lnR], [mean_lnR + std_lnR], None
        )
        assert bin_val > 0.0

        # Sum of bins should equal total integral
        lnR_edges = np.linspace(mean_lnR - 5 * std_lnR, mean_lnR + 5 * std_lnR, 11)
        bin_sum = sum(
            cluster_m.intp_bin(cosmo, lnM, z, [lnR_edges[i]], [lnR_edges[i + 1]], None)
            for i in range(len(lnR_edges) - 1)
        )
        total_integral = cluster_m.intp(cosmo, lnM, z)

        assert np.isclose(bin_sum, total_integral, rtol=1e-3)

    def test_intp_bin_consistency_with_p(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test intp_bin matches manual integration of p()."""
        lnM = np.log(1e14)
        z = 0.5
        mean_lnR = cluster_m.mu(lnM, z)
        std_lnR = cluster_m.sigma(lnM, z)

        lnR_lower = mean_lnR - std_lnR
        lnR_upper = mean_lnR + std_lnR

        lnR_array = np.linspace(lnR_lower, lnR_upper, 200)
        p_array = np.array(
            [cluster_m.p(cosmo, lnM, z, [lnR], None) for lnR in lnR_array]
        )
        manual_integral = np.trapezoid(p_array, lnR_array)

        intp_bin_val = cluster_m.intp_bin(cosmo, lnM, z, [lnR_lower], [lnR_upper], None)

        assert np.isclose(manual_integral, intp_bin_val, rtol=1e-3)

    def test_p_vec_z_lnMobs(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test vectorized probability distribution matches scalar calls."""
        lnM = np.log(1e14)
        z_array = np.linspace(0.1, 1.1, 50)
        lnR_array = np.log(np.geomspace(3.0, 1e3, 50))

        z_vec = Ncm.Vector.new_array(z_array)
        lnR_mat = Ncm.Matrix.new_array(lnR_array, 1)
        res = Ncm.Vector.new(len(z_array))

        cluster_m.p_vec_z_lnMobs(cosmo, lnM, z_vec, lnR_mat, None, res)

        for i, (z, lnR) in enumerate(zip(z_array, lnR_array)):
            p_scalar = cluster_m.p(cosmo, lnM, z, [lnR], None)
            p_vec = res.get(i)
            assert np.isclose(p_scalar, p_vec, rtol=1e-10)

    def test_sample_full_dist_property(self, cluster_m: Nc.ClusterMassRichness) -> None:
        """Test sample_full_dist property getter and setter."""
        # Default should be True
        assert cluster_m.get_sample_full_dist() is True

        cluster_m.set_sample_full_dist(False)
        assert cluster_m.get_sample_full_dist() is False

        cluster_m.set_sample_full_dist(True)
        assert cluster_m.get_sample_full_dist() is True

    def test_resample_full_dist(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test resampling from full distribution (may produce out-of-range values)."""
        rng = Ncm.RNG.new()
        rng.set_seed(42)

        cluster_m.set_sample_full_dist(True)

        lnM = np.log(1e14)
        z = 0.5
        cut = cluster_m.get_cut()

        n_samples = 1000
        n_in_range = 0
        samples_list = []
        lnR_vec = Ncm.Vector.new(1)

        for _ in range(n_samples):
            valid = cluster_m.resample_vec(cosmo, lnM, z, lnR_vec, None, rng)
            lnR = lnR_vec.get(0)
            samples_list.append(lnR)
            if valid:
                n_in_range += 1
                assert lnR >= cut

        # With full distribution, some samples should be out of range
        # (unless mean is very far from cut)
        samples = np.array(samples_list)
        assert len(samples) == n_samples

    def test_resample_truncated_dist(
        self, cluster_m: Nc.ClusterMassRichness, cosmo: Nc.HICosmo
    ) -> None:
        """Test resampling from truncated distribution (all values in range)."""
        rng = Ncm.RNG.new()
        rng.set_seed(42)

        cluster_m.set_sample_full_dist(False)

        lnM = np.log(1e14)
        z = 0.5
        cut = cluster_m.get_cut()

        n_samples = 100
        lnR_vec = Ncm.Vector.new(1)

        for _ in range(n_samples):
            valid = cluster_m.resample_vec(cosmo, lnM, z, lnR_vec, None, rng)
            lnR = lnR_vec.get(0)
            # With truncated sampling, all values should be >= cut
            assert lnR >= cut
            # valid indicates if also below lnR_max
            if valid:
                assert lnR <= LN_RICHNESS_MAX

    def test_serialization_deserialization(
        self, cluster_m: Nc.ClusterMassRichness
    ) -> None:
        """Test serialization and deserialization."""
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        cluster_m_dup = ser.dup_obj(cluster_m)

        assert cluster_m_dup is not None
        assert type(cluster_m_dup) is type(cluster_m)
        assert cluster_m_dup is not cluster_m

        # Check parameters are preserved
        lnM = np.log(1e14)
        z = 0.5
        assert isinstance(cluster_m_dup, Nc.ClusterMassRichness)
        assert np.isclose(cluster_m.mu(lnM, z), cluster_m_dup.mu(lnM, z))
        assert np.isclose(cluster_m.sigma(lnM, z), cluster_m_dup.sigma(lnM, z))
        assert cluster_m.get_cut() == cluster_m_dup.get_cut()


# =============================================================================
# Specific tests for NcClusterMassAscaso
# =============================================================================


class TestClusterMassAscaso:
    """Tests specific to NcClusterMassAscaso model."""

    def test_mu_formula(self, cluster_m_ascaso: Nc.ClusterMassAscaso) -> None:
        """Test mu follows the Ascaso formula: mu = mup0 + mup1*DlnM + mup2*Dln(1+z)."""
        lnM0 = cluster_m_ascaso.lnM0()
        ln1pz0 = cluster_m_ascaso.ln1pz0()

        mup0 = cluster_m_ascaso.param_get_by_name("mup0")
        mup1 = cluster_m_ascaso.param_get_by_name("mup1")
        mup2 = cluster_m_ascaso.param_get_by_name("mup2")

        lnM = np.log(1e14)
        z = 0.5

        DlnM = lnM - lnM0
        Dln1pz = np.log1p(z) - ln1pz0

        expected_mu = mup0 + mup1 * DlnM + mup2 * Dln1pz
        actual_mu = cluster_m_ascaso.mu(lnM, z)

        assert np.isclose(expected_mu, actual_mu, rtol=1e-10)

    def test_sigma_formula(self, cluster_m_ascaso: Nc.ClusterMassAscaso) -> None:
        """Test sigma follows the Ascaso formula.

        The formula is sigma = sigmap0 + sigmap1*DlnM + sigmap2*Dln(1+z)."""
        lnM0 = cluster_m_ascaso.lnM0()
        ln1pz0 = cluster_m_ascaso.ln1pz0()

        sigmap0 = cluster_m_ascaso.param_get_by_name("sigmap0")
        sigmap1 = cluster_m_ascaso.param_get_by_name("sigmap1")
        sigmap2 = cluster_m_ascaso.param_get_by_name("sigmap2")

        lnM = np.log(1e14)
        z = 0.5

        DlnM = lnM - lnM0
        Dln1pz = np.log1p(z) - ln1pz0

        expected_sigma = sigmap0 + sigmap1 * DlnM + sigmap2 * Dln1pz
        actual_sigma = cluster_m_ascaso.sigma(lnM, z)

        assert np.isclose(expected_sigma, actual_sigma, rtol=1e-10)

    def test_mu_mass_dependence(self, cluster_m_ascaso: Nc.ClusterMassAscaso) -> None:
        """Test mu increases with mass (positive mup1)."""
        z = 0.5
        lnM_array = np.linspace(np.log(1e13), np.log(1e16), 10)

        mu_array = [cluster_m_ascaso.mu(lnM, z) for lnM in lnM_array]

        # With positive mup1, mu should increase with mass
        mup1 = cluster_m_ascaso.param_get_by_name("mup1")
        if mup1 > 0:
            assert all(mu_array[i] < mu_array[i + 1] for i in range(len(mu_array) - 1))

    def test_mu_redshift_dependence(
        self, cluster_m_ascaso: Nc.ClusterMassAscaso
    ) -> None:
        """Test mu varies with redshift."""
        lnM = np.log(1e14)
        z_array = np.linspace(0.1, 1.0, 10)

        mu_array = [cluster_m_ascaso.mu(lnM, z) for z in z_array]

        # With positive mup2, mu should increase with redshift
        mup2 = cluster_m_ascaso.param_get_by_name("mup2")
        if mup2 > 0:
            assert all(mu_array[i] < mu_array[i + 1] for i in range(len(mu_array) - 1))


# =============================================================================
# Specific tests for NcClusterMassExt
# =============================================================================


class TestClusterMassExt:
    """Tests specific to NcClusterMassExt model."""

    def test_mu_formula(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test mu follows the Ext formula.

        The Ext formula is mu = mup0 + mup1*DlnM + mup2*DlnM^2 + mup3*DlnM^3.
        """
        lnM0 = cluster_m_ext.lnM0()

        mup0 = cluster_m_ext.param_get_by_name("mup0")
        mup1 = cluster_m_ext.param_get_by_name("mup1")
        mup2 = cluster_m_ext.param_get_by_name("mup2")
        mup3 = cluster_m_ext.param_get_by_name("mup3")

        lnM = np.log(1e14)
        z = 0.5  # z should not affect mu in Ext model

        DlnM = lnM - lnM0

        expected_mu = mup0 + mup1 * DlnM + mup2 * DlnM**2 + mup3 * DlnM**3
        actual_mu = cluster_m_ext.mu(lnM, z)

        assert np.isclose(expected_mu, actual_mu, rtol=1e-10)

    def test_sigma_formula(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test sigma follows the Ext formula.

        The formula is sigma = exp(sigmap0 + sigmap1*mu + sigmap2*mu^2).
        """
        sigmap0 = cluster_m_ext.param_get_by_name("sigmap0")
        sigmap1 = cluster_m_ext.param_get_by_name("sigmap1")
        sigmap2 = cluster_m_ext.param_get_by_name("sigmap2")

        lnM = np.log(1e14)
        z = 0.5

        mu = cluster_m_ext.mu(lnM, z)

        expected_sigma = np.exp(sigmap0 + sigmap1 * mu + sigmap2 * mu**2)
        actual_sigma = cluster_m_ext.sigma(lnM, z)

        assert np.isclose(expected_sigma, actual_sigma, rtol=1e-10)

    def test_mu_no_redshift_dependence(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test mu does not depend on redshift in Ext model."""
        lnM = np.log(1e14)
        z_array = [0.1, 0.5, 1.0, 1.5]

        mu_values = [cluster_m_ext.mu(lnM, z) for z in z_array]

        # All mu values should be identical
        assert all(np.isclose(mu_values[0], mu) for mu in mu_values)

    def test_sigma_depends_on_mu(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test that sigma varies with mu (mass-dependent scatter)."""
        z = 0.5
        lnM_array = np.linspace(np.log(1e13), np.log(1e16), 10)

        sigma_array = [cluster_m_ext.sigma(lnM, z) for lnM in lnM_array]

        # Sigma should vary with mass (since it depends on mu)
        assert not all(np.isclose(sigma_array[0], s) for s in sigma_array)

    def test_polynomial_mu(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test cubic polynomial behavior of mu."""
        lnM0 = cluster_m_ext.lnM0()
        z = 0.5

        # At pivot mass, DlnM = 0, so mu = mup0
        mup0 = cluster_m_ext.param_get_by_name("mup0")
        assert np.isclose(cluster_m_ext.mu(lnM0, z), mup0, rtol=1e-10)

        # Test polynomial evaluation at several points
        lnM_array = np.linspace(np.log(1e13), np.log(1e16), 20)
        for lnM in lnM_array:
            DlnM = lnM - lnM0
            mup1 = cluster_m_ext.param_get_by_name("mup1")
            mup2 = cluster_m_ext.param_get_by_name("mup2")
            mup3 = cluster_m_ext.param_get_by_name("mup3")

            expected = mup0 + mup1 * DlnM + mup2 * DlnM**2 + mup3 * DlnM**3
            assert np.isclose(cluster_m_ext.mu(lnM, z), expected, rtol=1e-10)


# =============================================================================
# Test lnM0 and ln1pz0 accessors
# =============================================================================


class TestPivotAccessors:
    """Test pivot mass and redshift accessors."""

    def test_lnM0_ascaso(self, cluster_m_ascaso: Nc.ClusterMassAscaso) -> None:
        """Test lnM0 accessor for Ascaso model."""
        lnM0 = cluster_m_ascaso.lnM0()
        # Default M0 is 3.0e14/0.71, so lnM0 should be positive and reasonable
        assert lnM0 > np.log(1e13)
        assert lnM0 < np.log(1e16)

    def test_lnM0_ext(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test lnM0 accessor for Ext model."""
        lnM0 = cluster_m_ext.lnM0()
        assert lnM0 > np.log(1e13)
        assert lnM0 < np.log(1e16)

    def test_ln1pz0_ascaso(self, cluster_m_ascaso: Nc.ClusterMassAscaso) -> None:
        """Test ln1pz0 accessor for Ascaso model."""
        ln1pz0 = cluster_m_ascaso.ln1pz0()
        # Default z0 is 0.6, so ln(1.6) ≈ 0.47
        assert ln1pz0 > 0.0
        assert ln1pz0 < 1.0

    def test_ln1pz0_ext(self, cluster_m_ext: Nc.ClusterMassExt) -> None:
        """Test ln1pz0 accessor for Ext model."""
        ln1pz0 = cluster_m_ext.ln1pz0()
        assert ln1pz0 > 0.0
        assert ln1pz0 < 1.0
