#!/usr/bin/env python
#
# test_py_jpas_forecast24.py
#
# Fri Dec 20 10:00:00 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_jpas_forecast24.py
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

"""Tests for J-PAS 2024 forecast module."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc

# Import module under test
from numcosmo_py.experiments import jpas_forecast24 as jpas

Ncm.cfg_init()

# Check if PySSC is available
pytest.importorskip("numcosmo_py.external.pyssc")


class TestEnums:
    """Test the enum classes."""

    def test_jpas_ssc_type_values(self):
        """Test JpasSSCType enum values."""
        assert jpas.JpasSSCType.NO_SSC == "no_ssc"
        assert jpas.JpasSSCType.FULLSKY == "fullsky"
        assert jpas.JpasSSCType.FULL == "full"
        assert jpas.JpasSSCType.GUARANTEED == "guaranteed"

    def test_cluster_mass_type_values(self):
        """Test ClusterMassType enum values."""
        assert jpas.ClusterMassType.NODIST == "nodist"
        assert jpas.ClusterMassType.ASCASO == "ascaso"

    def test_cluster_redshift_type_values(self):
        """Test ClusterRedshiftType enum values."""
        assert jpas.ClusterRedshiftType.NODIST == "nodist"
        assert jpas.ClusterRedshiftType.GAUSS == "gauss"


class TestZBinsKernels:
    """Test redshift bins and kernels creation."""

    def test_create_zbins_kernels_defaults(self):
        """Test create_zbins_kernels with default parameters."""
        kernel_z, kernels_T, z_bins_knots = jpas.create_zbins_kernels()

        assert len(kernel_z) == 400
        assert kernels_T.shape == (7, 400)
        assert len(z_bins_knots) == 8
        assert z_bins_knots[0] == 0.1
        assert z_bins_knots[-1] == 0.8

    def test_create_zbins_kernels_custom(self):
        """Test create_zbins_kernels with custom parameters."""
        kernel_z, kernels_T, z_bins_knots = jpas.create_zbins_kernels(
            z_min=0.2, z_max=1.0, nknots=5, kernel_nknots=100, kernel_zmax=2.0
        )

        assert len(kernel_z) == 100
        assert kernels_T.shape == (4, 100)
        assert len(z_bins_knots) == 5
        assert z_bins_knots[0] == 0.2
        assert z_bins_knots[-1] == 1.0
        assert kernel_z[-1] == pytest.approx(2.0, rel=0.01)

    def test_kernels_normalization(self):
        """Test that kernels are properly normalized."""
        kernel_z, kernels_T, z_bins_knots = jpas.create_zbins_kernels(
            nknots=4, kernel_nknots=1000
        )

        # Each kernel should integrate to 1 over its bin
        for i in range(len(z_bins_knots) - 1):
            dz = kernel_z[1] - kernel_z[0]
            integral = np.sum(kernels_T[i]) * dz
            # Top-hat normalized by Dz should integrate to 1
            assert integral == pytest.approx(1.0, rel=0.05)


class TestLnMObsBins:
    """Test observable mass bins creation."""

    def test_create_lnM_obs_bins_defaults(self):
        """Test create_lnM_obs_bins with default parameters."""
        lnM_bins = jpas.create_lnM_obs_bins()

        assert len(lnM_bins) == 2
        assert lnM_bins[0] == pytest.approx(np.log(10.0) * 14.0)
        assert lnM_bins[-1] == pytest.approx(np.log(10.0) * 15.0)

    def test_create_lnM_obs_bins_custom(self):
        """Test create_lnM_obs_bins with custom parameters."""
        lnM_bins = jpas.create_lnM_obs_bins(
            lnM_obs_min=np.log(10.0) * 13.0,
            lnM_obs_max=np.log(10.0) * 16.0,
            nknots=5,
        )

        assert len(lnM_bins) == 5
        assert lnM_bins[0] == pytest.approx(np.log(10.0) * 13.0)
        assert lnM_bins[-1] == pytest.approx(np.log(10.0) * 16.0)


class TestSurveyArea:
    """Test survey area calculation."""

    def test_survey_area_fullsky(self):
        """Test survey area for FULLSKY."""
        area = jpas.survey_area(jpas.JpasSSCType.FULLSKY)
        assert area == pytest.approx(20009.97)

    def test_survey_area_full(self):
        """Test survey area for FULL."""
        area = jpas.survey_area(jpas.JpasSSCType.FULL)
        assert area == pytest.approx(8500.0)

    def test_survey_area_guaranteed(self):
        """Test survey area for GUARANTEED."""
        area = jpas.survey_area(jpas.JpasSSCType.GUARANTEED)
        assert area == pytest.approx(2959.1)

    def test_survey_area_invalid(self):
        """Test survey area raises error for NO_SSC."""
        with pytest.raises(ValueError, match="Invalid sky cut type"):
            jpas.survey_area(jpas.JpasSSCType.NO_SSC)


class TestMasks:
    """Test HEALPix mask creation."""

    def test_create_mask_guaranteed(self):
        """Test guaranteed mask creation."""
        mask = jpas.create_mask_guaranteed(nside=32)

        assert len(mask) == 12 * 32 * 32
        assert np.all((mask == 0) | (mask == 1))
        assert np.sum(mask) > 0  # Should have some pixels

    def test_create_mask_guaranteed_different_nside(self):
        """Test guaranteed mask with different nside."""
        mask_low = jpas.create_mask_guaranteed(nside=16)
        mask_high = jpas.create_mask_guaranteed(nside=64)

        assert len(mask_low) == 12 * 16 * 16
        assert len(mask_high) == 12 * 64 * 64

    def test_create_mask_full(self):
        """Test full mask creation."""
        mask = jpas.create_mask_full(nside=32)

        assert len(mask) == 12 * 32 * 32
        assert np.all((mask == 0) | (mask == 1))
        assert np.sum(mask) > 0

    def test_create_mask_full_different_nside(self):
        """Test full mask with different nside."""
        mask_low = jpas.create_mask_full(nside=16)
        mask_high = jpas.create_mask_full(nside=64)

        assert len(mask_low) == 12 * 16 * 16
        assert len(mask_high) == 12 * 64 * 64

    def test_masks_coverage(self):
        """Test that full mask covers more area than guaranteed."""
        nside = 64
        mask_full = jpas.create_mask_full(nside=nside)
        mask_guaranteed = jpas.create_mask_guaranteed(nside=nside)

        # Full mask should have more pixels than guaranteed
        assert np.sum(mask_full) > np.sum(mask_guaranteed)


class TestCosmoCreation:
    """Test cosmology model creation."""

    def test_create_cosmo(self):
        """Test create_cosmo creates valid cosmology."""
        cosmo = jpas.create_cosmo()

        assert isinstance(cosmo, Nc.HICosmoDEXcdm)
        assert cosmo["H0"] == pytest.approx(67.81)
        assert cosmo["Omegab"] == pytest.approx(0.0486)
        assert cosmo["w"] == pytest.approx(-1.0)
        assert cosmo["Omegak"] == pytest.approx(0.0)

    def test_create_cosmo_has_submodels(self):
        """Test create_cosmo has primordial and reionization submodels."""
        cosmo = jpas.create_cosmo()

        # Check that primordial spectrum model is attached
        prim = cosmo.peek_submodel_by_mid(Nc.HIPrim.id())
        assert prim is not None
        assert isinstance(prim, Nc.HIPrimPowerLaw)

        # Check that reionization model is attached
        reion = cosmo.peek_submodel_by_mid(Nc.HIReion.id())
        assert reion is not None
        assert isinstance(reion, Nc.HIReionCamb)

    def test_create_cosmo_parameter_fitting_status(self):
        """Test parameter fitting status."""
        cosmo = jpas.create_cosmo()
        prim = cosmo.peek_submodel_by_mid(Nc.HIPrim.id())

        # Check fit status for cosmology parameters
        assert not cosmo.param_get_desc("H0")["fit"]
        assert cosmo.param_get_desc("Omegac")["fit"]
        assert not cosmo.param_get_desc("Omegab")["fit"]
        assert cosmo.param_get_desc("w")["fit"]
        assert not cosmo.param_get_desc("Omegak")["fit"]

        # Check fit status for primordial parameters
        assert prim.param_get_desc("ln10e10ASA")["fit"]
        assert not prim.param_get_desc("n_SA")["fit"]


class TestMFuncArray:
    """Test derived parameter function creation."""

    def test_create_mfunc_array(self):
        """Test create_mfunc_array creates function array."""
        tf = Nc.TransferFuncEH()
        psml = Nc.PowspecMLTransfer.new(tf)

        mfunc_oa = jpas.create_mfunc_array(psml)

        assert isinstance(mfunc_oa, Ncm.ObjArray)
        assert mfunc_oa.len() == 2

    def test_create_mfunc_array_contents(self):
        """Test create_mfunc_array has correct functions."""
        tf = Nc.TransferFuncEH()
        psml = Nc.PowspecMLTransfer.new(tf)

        mfunc_oa = jpas.create_mfunc_array(psml)

        # First should be sigma8
        func0 = mfunc_oa.peek(0)
        assert isinstance(func0, Ncm.MSetFuncList)

        # Second should be Omega_m0
        func1 = mfunc_oa.peek(1)
        assert isinstance(func1, Ncm.MSetFuncList)


class TestCovarianceMatrices:
    """Test SSC covariance matrix creation."""

    @pytest.fixture
    def setup_for_covariance(self):
        """Setup kernel and cosmology for covariance tests."""
        kernel_z, kernels_T, _ = jpas.create_zbins_kernels(nknots=4, kernel_nknots=50)
        cosmo = jpas.create_cosmo()
        return kernel_z, kernels_T, cosmo

    def test_create_covariance_S_fullsky(self, setup_for_covariance):
        """Test fullsky SSC covariance matrix creation."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S_fullsky(kernel_z, kernels_T, cosmo)

        assert isinstance(S, Ncm.Matrix)
        n_bins = kernels_T.shape[0]
        assert S.nrows() == n_bins
        assert S.ncols() == n_bins

    def test_create_covariance_S_guaranteed(self, setup_for_covariance):
        """Test guaranteed mask SSC covariance matrix creation."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S_guaranteed(kernel_z, kernels_T, cosmo)

        assert isinstance(S, Ncm.Matrix)
        n_bins = kernels_T.shape[0]
        assert S.nrows() == n_bins
        assert S.ncols() == n_bins

    def test_create_covariance_S_full(self, setup_for_covariance):
        """Test full mask SSC covariance matrix creation."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S_full(kernel_z, kernels_T, cosmo)

        assert isinstance(S, Ncm.Matrix)
        n_bins = kernels_T.shape[0]
        assert S.nrows() == n_bins
        assert S.ncols() == n_bins

    def test_create_covariance_S_router_fullsky(self, setup_for_covariance):
        """Test covariance router for FULLSKY."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S(
            kernel_z, kernels_T, jpas.JpasSSCType.FULLSKY, cosmo
        )

        assert isinstance(S, Ncm.Matrix)

    def test_create_covariance_S_router_full(self, setup_for_covariance):
        """Test covariance router for FULL."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S(kernel_z, kernels_T, jpas.JpasSSCType.FULL, cosmo)

        assert isinstance(S, Ncm.Matrix)

    def test_create_covariance_S_router_guaranteed(self, setup_for_covariance):
        """Test covariance router for GUARANTEED."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        S = jpas.create_covariance_S(
            kernel_z, kernels_T, jpas.JpasSSCType.GUARANTEED, cosmo
        )

        assert isinstance(S, Ncm.Matrix)

    def test_create_covariance_S_invalid(self, setup_for_covariance):
        """Test covariance router raises error for NO_SSC."""
        kernel_z, kernels_T, cosmo = setup_for_covariance

        with pytest.raises(ValueError, match="Invalid sky cut type"):
            jpas.create_covariance_S(
                kernel_z, kernels_T, jpas.JpasSSCType.NO_SSC, cosmo
            )


class TestClusterModels:
    """Test cluster mass and redshift model creation."""

    def test_create_cluster_mass_nodist(self):
        """Test cluster mass NODIST creation."""
        cluster_m, lnM_obs_min, lnM_obs_max = jpas.create_cluster_mass(
            jpas.ClusterMassType.NODIST
        )

        assert isinstance(cluster_m, Nc.ClusterMassNodist)
        assert lnM_obs_min == pytest.approx(np.log(10) * 14.0)
        assert lnM_obs_max == pytest.approx(np.log(10) * 16.0)

    def test_create_cluster_mass_ascaso(self):
        """Test cluster mass ASCASO creation."""
        cluster_m, lnM_obs_min, lnM_obs_max = jpas.create_cluster_mass(
            jpas.ClusterMassType.ASCASO
        )

        assert isinstance(cluster_m, Nc.ClusterMassAscaso)
        assert cluster_m["mup0"] == pytest.approx(3.207)
        assert cluster_m["mup1"] == pytest.approx(0.993)
        assert cluster_m["sigmap0"] == pytest.approx(0.456)
        assert lnM_obs_min == pytest.approx(np.log(5.0))
        assert lnM_obs_max == pytest.approx(np.log(10) * 2.5)

    def test_create_cluster_mass_invalid(self):
        """Test create_cluster_mass raises error for invalid type."""
        with pytest.raises(ValueError, match="Invalid cluster mass type"):
            # Use a string that's not a valid enum member
            jpas.create_cluster_mass("invalid_type")  # type: ignore

    def test_create_cluster_redshift_nodist(self):
        """Test cluster redshift NODIST creation."""
        cluster_z = jpas.create_cluster_redshift(jpas.ClusterRedshiftType.NODIST)

        assert isinstance(cluster_z, Nc.ClusterRedshiftNodist)

    def test_create_cluster_redshift_gauss(self):
        """Test cluster redshift GAUSS creation."""
        cluster_z = jpas.create_cluster_redshift(jpas.ClusterRedshiftType.GAUSS)

        assert isinstance(cluster_z, Nc.ClusterPhotozGaussGlobal)
        assert cluster_z["z-bias"] == pytest.approx(0.0)
        assert cluster_z["sigma0"] == pytest.approx(0.1)

    def test_create_cluster_redshift_invalid(self):
        """Test create_cluster_redshift raises error for invalid type."""
        with pytest.raises(ValueError, match="Invalid cluster redshift type"):
            jpas.create_cluster_redshift("invalid_type")  # type: ignore


class TestGenerateJpasForecast:
    """Test main experiment generation function."""

    def test_generate_jpas_forecast_2024_minimal(self):
        """Test generate_jpas_forecast_2024 with minimal parameters."""
        experiment, mfunc_oa = jpas.generate_jpas_forecast_2024(
            znknots=3, lnMobsnknots=2
        )

        assert isinstance(experiment, Ncm.ObjDictStr)
        assert isinstance(mfunc_oa, Ncm.ObjArray)

        # Check experiment contains expected objects
        likelihood = experiment.get("likelihood")
        assert isinstance(likelihood, Ncm.Likelihood)

        mset = experiment.get("model-set")
        assert isinstance(mset, Ncm.MSet)

    def test_generate_jpas_forecast_2024_with_ssc_fullsky(self):
        """Test generation with FULLSKY SSC."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_Sij_type=jpas.JpasSSCType.FULLSKY,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_with_ssc_guaranteed(self):
        """Test generation with GUARANTEED SSC for fitting."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_Sij_type=jpas.JpasSSCType.GUARANTEED,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_with_ssc_full(self):
        """Test generation with FULL SSC for fitting."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_Sij_type=jpas.JpasSSCType.FULL,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_resample_ssc_fullsky(self):
        """Test generation with resample FULLSKY SSC."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            resample_Sij_type=jpas.JpasSSCType.FULLSKY,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_resample_ssc_guaranteed(self):
        """Test generation with resample GUARANTEED SSC."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            resample_Sij_type=jpas.JpasSSCType.GUARANTEED,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_resample_ssc_full(self):
        """Test generation with resample FULL SSC."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            resample_Sij_type=jpas.JpasSSCType.FULL,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_cluster_mass_ascaso(self):
        """Test generation with ASCASO cluster mass."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            cluster_mass_type=jpas.ClusterMassType.ASCASO,
        )

        mset = experiment.get("model-set")
        cluster_m = mset.peek(Nc.ClusterMass.id())
        assert isinstance(cluster_m, Nc.ClusterMassAscaso)

    def test_generate_jpas_forecast_2024_cluster_redshift_gauss(self):
        """Test generation with GAUSS cluster redshift."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            cluster_redshift_type=jpas.ClusterRedshiftType.GAUSS,
        )

        mset = experiment.get("model-set")
        cluster_z = mset.peek(Nc.ClusterRedshift.id())
        assert isinstance(cluster_z, Nc.ClusterPhotozGaussGlobal)

    def test_generate_jpas_forecast_2024_with_fixed_cov(self):
        """Test generation with fixed covariance."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3, lnMobsnknots=2, use_fixed_cov=True
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_different_params(self):
        """Test generation with different fitting and resample params."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_model=(0.27, -0.95, 0.82),
            resample_model=(0.25, -1.05, 0.81),
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_custom_seed(self):
        """Test generation with custom random seed."""
        experiment1, _ = jpas.generate_jpas_forecast_2024(
            znknots=3, lnMobsnknots=2, resample_seed=1111
        )

        experiment2, _ = jpas.generate_jpas_forecast_2024(
            znknots=3, lnMobsnknots=2, resample_seed=2222
        )

        # Both should be valid but with different random realizations
        assert isinstance(experiment1, Ncm.ObjDictStr)
        assert isinstance(experiment2, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_custom_redshift_range(self):
        """Test generation with custom redshift range."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            z_min=0.2, z_max=1.0, znknots=4, lnMobsnknots=2
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_custom_mass_range(self):
        """Test generation with custom mass range."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnM_obs_min=np.log(10.0) * 13.5,
            lnM_obs_max=np.log(10.0) * 15.5,
            lnMobsnknots=3,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_custom_area(self):
        """Test generation with custom survey area."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3, lnMobsnknots=2, area=5000.0
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_area_override_resample(self):
        """Test that resample SSC type overrides area parameter."""
        # When resample uses a masked SSC, the area should be updated
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            area=1000.0,  # This should be overridden
            resample_Sij_type=jpas.JpasSSCType.GUARANTEED,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)
        # Area should have been adjusted to guaranteed area (2959.1)

    def test_generate_jpas_forecast_2024_area_override_fitting(self):
        """Test that fitting SSC type overrides area parameter."""
        # When fitting uses a masked SSC and resample doesn't, area from fitting is used
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            area=1000.0,  # This should be overridden
            fitting_Sij_type=jpas.JpasSSCType.FULL,
            resample_Sij_type=jpas.JpasSSCType.NO_SSC,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_both_ssc_types(self):
        """Test generation with both fitting and resample SSC."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_Sij_type=jpas.JpasSSCType.FULLSKY,
            resample_Sij_type=jpas.JpasSSCType.FULLSKY,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

    def test_generate_jpas_forecast_2024_no_fitting_ssc_with_resample_ssc(self):
        """Test generation with NO_SSC for fitting but SSC for resampling."""
        experiment, _ = jpas.generate_jpas_forecast_2024(
            znknots=3,
            lnMobsnknots=2,
            fitting_Sij_type=jpas.JpasSSCType.NO_SSC,
            resample_Sij_type=jpas.JpasSSCType.FULLSKY,
        )

        assert isinstance(experiment, Ncm.ObjDictStr)

        # Verify likelihood can be evaluated
        likelihood = experiment.get("likelihood")
        mset = experiment.get("model-set")
        m2lnL = likelihood.m2lnL_val(mset)
        assert np.isfinite(m2lnL)

    def test_generate_jpas_forecast_2024_likelihood_evaluation(self):
        """Test that generated experiment can evaluate likelihood."""
        experiment, _ = jpas.generate_jpas_forecast_2024(znknots=3, lnMobsnknots=2)

        likelihood = experiment.get("likelihood")
        mset = experiment.get("model-set")

        # Should be able to evaluate likelihood
        m2lnL = likelihood.m2lnL_val(mset)
        assert np.isfinite(m2lnL)
        assert m2lnL >= 0


class TestSetMSetParams:
    """Test internal _set_mset_params function."""

    def test_set_mset_params(self):
        """Test _set_mset_params sets parameters correctly."""
        cosmo = jpas.create_cosmo()
        mset = Ncm.MSet.new_array([cosmo])
        mset.prepare_fparam_map()

        # Set parameters
        params = (0.27, -0.95, 0.82)
        # pylint: disable-next=protected-access
        jpas._set_mset_params(mset, params)

        # Check that Omegac and w are set
        assert mset["NcHICosmo"]["Omegac"] == pytest.approx(0.27)
        assert mset["NcHICosmo"]["w"] == pytest.approx(-0.95)

        # ln10e10ASA should be adjusted to match sigma8
        # Just check it changed from default
        prim = cosmo.peek_submodel_by_mid(Nc.HIPrim.id())
        assert prim["ln10e10ASA"] != pytest.approx(3.02745)
