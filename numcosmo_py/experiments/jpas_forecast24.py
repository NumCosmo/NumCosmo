#
# jpas_forecast24.py
#
# Mon Feb 20 22:31:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# jpas_forecast24.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Factory functions for J-Pas 2024 forcasting.

This module provides factory functions to create the J-Pas 2024 forecast
experiment dictionary and the extra functions for the forecast.
"""

from typing import cast
from enum import Enum
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.helper import npa_to_seq
from numcosmo_py.external.pyssc import pyssc as PySSC


class JpasSSCType(str, Enum):
    """J-Pas 2024 Super Sample Covariance types."""

    NO_SSC = "no_ssc"
    FULLSKY = "fullsky"
    FULL = "full"
    GUARANTEED = "guaranteed"


class ClusterMassType(str, Enum):
    """Mass-observable relation types."""

    NODIST = "nodist"
    ASCASO = "ascaso"


class ClusterRedshiftType(str, Enum):
    """Photoz types."""

    NODIST = "nodist"
    GAUSS = "gauss"


def create_zbins_kernels(
    z_min: float = 0.1,
    z_max: float = 0.8,
    nknots: int = 8,
    kernel_nknots: int = 400,
    kernel_zmax: float = 1.9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create the redshift bins and kernels for the J-Pas 2024 forecast."""
    # Redshift bins and kernels
    z_bins_len = nknots - 1
    z_bins_knots = np.linspace(z_min, z_max, num=nknots)

    kernel_z = np.linspace(0.0, kernel_zmax, num=kernel_nknots + 1)[1:]

    kernels_T = np.zeros((z_bins_len, kernel_nknots))

    for i, (zminbin, zmaxbin) in enumerate(zip(z_bins_knots[:-1], z_bins_knots[1:])):
        Dz = zmaxbin - zminbin

        kernel = np.zeros_like(kernel_z)
        kernel[(kernel_z >= zminbin) & (kernel_z <= zmaxbin)] = 1.0
        kernels_T[i] = kernel / Dz

    return kernel_z, kernels_T, z_bins_knots


def create_lnM_bins(
    lnM_min: float = np.log(10.0) * 14.0,
    lnM_max: float = np.log(10.0) * 15.0,
    nknots: int = 2,
) -> np.ndarray:
    """Create the mass bins for the J-Pas 2024 forecast."""
    lnM_bins_knots = np.linspace(lnM_min, lnM_max, nknots)

    return lnM_bins_knots


def survey_area(sky_cut: JpasSSCType):
    """Return the survey area for the J-Pas 2024 forecast."""
    match sky_cut:
        case JpasSSCType.FULLSKY:
            survey_area0 = 20009.97
        case JpasSSCType.GUARANTEED:
            survey_area0 = 2959.1
        case _:
            raise ValueError(f"Invalid sky cut type: {sky_cut}")

    return survey_area0


def create_mask_guaranteed(nside: int = 512) -> np.ndarray:
    """Create the mask for the J-Pas 2024 forecast."""
    smap = Ncm.SphereMap.new(nside)

    npix = smap.get_npix()

    angles = np.array([smap.pix2ang_ring(i) for i in range(npix)])

    mask1_guaranteed = np.zeros(npix)
    mask2_guaranteed = np.zeros(npix)
    mask3_guaranteed = np.zeros(npix)

    pix_theta_ecl = angles[:, 0]
    pix_phi_ecl = angles[:, 1]

    mask1_guaranteed_condition = (
        (pix_phi_ecl > 3.0 * np.pi / 4.0)
        & (pix_phi_ecl < 13.0 * np.pi / 12.0)
        & (pix_theta_ecl < np.pi / 2.0 - 30.0 * np.pi / 180)
        & (pix_theta_ecl > np.pi / 2.0 - 80.0 * np.pi / 180)
    )
    mask2_guaranteed_condition = (
        (pix_phi_ecl > 11.0 * np.pi / 6.0)
        & (pix_theta_ecl < np.pi / 2.0 - 20.0 * np.pi / 180)
        & (pix_theta_ecl > np.pi / 2.0 - 40.0 * np.pi / 180)
    )
    mask3_guaranteed_condition = (
        (pix_phi_ecl < np.pi / 4.0)
        & (pix_theta_ecl < np.pi / 2.0 - 20.0 * np.pi / 180)
        & (pix_theta_ecl > np.pi / 2.0 - 40.0 * np.pi / 180)
    )

    mask1_guaranteed[mask1_guaranteed_condition] = 1
    mask2_guaranteed[mask2_guaranteed_condition] = 1
    mask3_guaranteed[mask3_guaranteed_condition] = 1

    mask_guaranteed = mask1_guaranteed + mask2_guaranteed + mask3_guaranteed

    return mask_guaranteed


def create_mask_full(nside: int = 512) -> np.ndarray:
    """Create the mask for the J-Pas 2024 forecast."""
    smap = Ncm.SphereMap.new(nside)

    npix = smap.get_npix()

    angles = np.array([smap.pix2ang_ring(i) for i in range(npix)])
    pix_theta_ecl = angles[:, 0]
    pix_phi_ecl = angles[:, 1]

    # Total mask
    mask1_full = np.zeros(npix)
    mask2_full = np.zeros(npix)
    mask3_full = np.zeros(npix)

    mask1_full_condition = (
        (pix_phi_ecl > 2.0 * np.pi / 3.0)
        & (pix_phi_ecl < 3.0 * np.pi / 2.0)
        & (pix_theta_ecl < np.pi / 2.0 - 10.0 * np.pi / 180)
        & (pix_theta_ecl > np.pi / 2.0 - 80.0 * np.pi / 180)
    )
    mask2_full_condition = (
        (pix_phi_ecl > 11.0 * np.pi / 6.0)
        & (pix_theta_ecl < np.pi / 2.0)
        & (pix_theta_ecl > np.pi / 2.0 - 45.0 * np.pi / 180)
    )
    mask3_full_condition = (
        (pix_phi_ecl < np.pi / 4.0)
        & (pix_theta_ecl < np.pi / 2.0)
        & (pix_theta_ecl > np.pi / 2.0 - 45.0 * np.pi / 180)
    )

    mask1_full[mask1_full_condition] = 1
    mask2_full[mask2_full_condition] = 1
    mask3_full[mask3_full_condition] = 1

    mask_full = mask1_full + mask2_full + mask3_full

    return mask_full


def create_cosmo() -> Nc.HICosmo:
    """Create a cosmology for J-Pas 2024 forecast."""
    cosmo = Nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.omega_x2omega_k()
    cosmo["H0"] = 67.81
    cosmo["Omegab"] = 0.0486
    cosmo["Omegac"] = 0.2612
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.00

    cosmo.param_set_desc("H0", {"fit": False})
    cosmo.param_set_desc("Omegac", {"fit": True})
    cosmo.param_set_desc("Omegab", {"fit": False})
    cosmo.param_set_desc("w", {"fit": True})
    cosmo.param_set_desc("Omegak", {"fit": False})

    prim = Nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745
    prim["n_SA"] = 0.9660

    prim.param_set_desc("ln10e10ASA", {"fit": True})
    prim.param_set_desc("n_SA", {"fit": False})

    reion = Nc.HIReionCamb.new()
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def create_mfunc_array(psml: Nc.PowspecML) -> Ncm.ObjArray:
    """Create a list of extra functions for J-Pas 2024 forecast."""
    mfunc_oa = Ncm.ObjArray.new()

    psml.require_kmin(1.0e-5)
    psml.require_kmax(1.0e1)
    psml.require_zi(0.0)
    psml.require_zf(1.0)

    psf_cbe = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf_cbe.set_best_lnr0()

    mfunc_sigma8 = Ncm.MSetFuncList.new("NcHICosmo:sigma8", psf_cbe)
    mfunc_oa.add(mfunc_sigma8)

    mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
    mfunc_oa.add(mfunc_Omegam)

    return mfunc_oa


def create_covariance_S_fullsky(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base covariance matrix S_ij  based on the full sky."""
    S_fullsky_array = PySSC.Sij(kernel_z, kernels_T, cosmo)

    S_fullsky = Ncm.Matrix.new_array(
        S_fullsky_array.flatten(), S_fullsky_array.shape[1]
    )

    return S_fullsky


def create_covariance_S_guaranteed(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base covariance matrix S_ij for the guaranteed mask."""
    mask = create_mask_guaranteed()

    S_guaranteed_array = PySSC.Sij_psky(kernel_z, kernels_T, cosmo, mask=mask)

    S_guaranteed = Ncm.Matrix.new_array(
        S_guaranteed_array.flatten(), S_guaranteed_array.shape[1]
    )

    return S_guaranteed


def create_covariance_S_full(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base covariance matrix S_ij for the full mask."""
    mask = create_mask_full()

    S_full_array = PySSC.Sij_psky(kernel_z, kernels_T, cosmo, mask=mask)

    S_full = Ncm.Matrix.new_array(S_full_array.flatten(), S_full_array.shape[1])

    return S_full


def create_covariance_S(
    kernel_z: np.ndarray,
    kernels_T: np.ndarray,
    sky_cut: JpasSSCType,
    cosmo: Nc.HICosmo,
) -> Ncm.Matrix:
    """Create the base covariance matrix S_ij."""
    if sky_cut == JpasSSCType.FULLSKY:
        S = create_covariance_S_fullsky(kernel_z, kernels_T, cosmo)
    elif sky_cut == JpasSSCType.FULL:
        S = create_covariance_S_full(kernel_z, kernels_T, cosmo)
    elif sky_cut == JpasSSCType.GUARANTEED:
        S = create_covariance_S_guaranteed(kernel_z, kernels_T, cosmo)
    else:
        raise ValueError(f"Invalid sky cut type: {sky_cut}")

    return S


def create_cluster_mass(
    cluster_mass_type: ClusterMassType,
) -> Nc.ClusterMass:
    """Create the cluster mass-observable relation."""
    if cluster_mass_type == ClusterMassType.NODIST:
        cluster_m = Nc.ClusterMassNodist(
            lnM_min=np.log(10) * 14.0, lnM_max=np.log(10) * 16.0
        )
    elif cluster_mass_type == ClusterMassType.ASCASO:
        cluster_m = Nc.ClusterMassAscaso(
            lnRichness_min=0.0, lnRichness_max=np.log(10) * 2.5, z0=0
        )
        cluster_m["mup0"] = 3.207
        cluster_m["mup1"] = 0.993
        cluster_m["mup2"] = 0
        cluster_m["sigmap0"] = 0.456
        cluster_m["sigmap1"] = -0.169  # valor do murata e 0 valor do SRD
        cluster_m["sigmap2"] = 0
        # [2.996 , 5.393 , 5]
    else:
        raise ValueError(f"Invalid cluster mass type: {cluster_mass_type}")

    return cluster_m


def create_cluster_redshift(
    cluster_redshift_type: ClusterRedshiftType,
) -> Nc.ClusterRedshift:
    """Create the cluster photoz relation."""
    if cluster_redshift_type == ClusterRedshiftType.NODIST:
        cluster_z = Nc.ClusterRedshiftNodist(z_min=0.0, z_max=2.0)
    elif cluster_redshift_type == ClusterRedshiftType.GAUSS:
        cluster_z = Nc.ClusterPhotozGaussGlobal(pz_min=0.0, pz_max=1.0)
        cluster_z["z-bias"] = 0
        cluster_z["sigma0"] = 0.1  # year1 0.003 year10

    else:
        raise ValueError(f"Invalid cluster redshift type: {cluster_redshift_type}")

    return cluster_z


def _set_mset_params(mset: Ncm.MSet, params: tuple[float, float, float]) -> None:
<<<<<<< HEAD
    """Set the parameters for the mass model."""
<<<<<<< HEAD
    param_names = ["NcHICosmo:Omegac", "NcHICosmo:w", "NcHIPrim:ln10e10ASA"]
    
=======
>>>>>>> 2de4259916cf3226ca3cc5cf0d35e0962ac87ac5
=======
    """Set the parameters for the cosmology model.

    :param mset: The mass model set.
    :param params: The parameters for the cosmology model, (Omegac, w, sigma8).
    """
>>>>>>> f81b5976cdb4356db9de13d8c1da699e3120f0b7
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)
<<<<<<< HEAD
    
    for i, param_name in enumerate(param_names):
        pi = mset.fparam_get_pi_by_name(param_name)
        
        
        
        if i ==2:
            
            cosmo = mset.peek(mset.get_id_by_ns("NcHICosmo"))
            A_s = np.exp(mset.param_get(pi.mid, pi.pid)) * 1.0e-10
            
            fact = (params[2] / psml.sigma_tophat_R(cosmo, 1.0e-7, 0.0, 8.0 / 0.6774)) ** 2
            mset.param_set(pi.mid, pi.pid,  np.log(1.0e10 * A_s * fact))
        else:
            mset.param_set(pi.mid, pi.pid, params[i])
        
=======
    Omegac, w, sigma8 = params

    cosmo: Nc.HICosmo = cast(Nc.HICosmo, mset["NcHICosmo"])
    prim: Nc.HIPrim = cast(Nc.HIPrim, mset["NcHIPrim"])
    mset["NcHICosmo"]["Omegac"] = Omegac
    mset["NcHICosmo"]["w"] = w

    A_s = np.exp(prim["ln10e10ASA"]) * 1.0e-10
    fact = (sigma8 / psml.sigma_tophat_R(cosmo, 1.0e-7, 0.0, 8.0 / cosmo.h())) ** 2
    prim["ln10e10ASA"] = np.log(1.0e10 * A_s * fact)

>>>>>>> 2de4259916cf3226ca3cc5cf0d35e0962ac87ac5

def generate_jpas_forecast_2024(
    area: float = 2959.1,
    z_min: float = 0.1,
    z_max: float = 0.8,
    znknots: int = 8,
    cluster_redshift_type: ClusterRedshiftType = ClusterRedshiftType.NODIST,
    lnM_min: float = np.log(10.0) * 14.0,
    lnM_max: float = np.log(10.0) * 15.0,
    lnMnknots: int = 2,
    cluster_mass_type: ClusterMassType = ClusterMassType.NODIST,
    use_fixed_cov: bool = False,
    fitting_model: tuple[float, float, float] = (0.2612, -1.0, 0.8159),
    resample_model: tuple[float, float, float] = (0.2612, -1.0, 0.8159),
    resample_seed: int = 1234,
    fitting_Sij_type: JpasSSCType = JpasSSCType.FULLSKY,
    resample_Sij_type: JpasSSCType = JpasSSCType.NO_SSC,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate J-Pas forecast 2024 experiment dictionary."""
    # Computation tools
    if (
        resample_Sij_type == JpasSSCType.FULL
        or resample_Sij_type == JpasSSCType.GUARANTEED
    ):
        area = survey_area(fitting_Sij_type)
        print(
            f"sky cut is of type {fitting_Sij_type} adjusting to "
            f"corresponding survey area = {area:.2f} sqd. "
        )

    dist = Nc.Distance.new(2.0)

    tf = Nc.TransferFuncEH()

    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)

    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)

    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hbias_Tinker = Nc.HaloBiasTinker.new(hmf)
    cad = Nc.ClusterAbundance.new(hmf, hbias_Tinker)
    cad.set_area(area * (np.pi / 180) ** 2)
    # Models

    cluster_m = create_cluster_mass(cluster_mass_type)
    cluster_z = create_cluster_redshift(cluster_redshift_type)
    cosmo = create_cosmo()

    mset = Ncm.MSet.new_array([cosmo, cluster_m, cluster_z])
    mset.prepare_fparam_map()

    # Likelihood
    #   Bins and kernels

    kernel_z, kernels_T, z_bins_knots = create_zbins_kernels(
        z_min=z_min, z_max=z_max, nknots=znknots
    )
    lnM_bins_knots = create_lnM_bins(lnM_min=lnM_min, lnM_max=lnM_max, nknots=lnMnknots)

    z_bins_vec = Ncm.Vector.new_array(npa_to_seq(z_bins_knots))
    lnM_bins_vec = Ncm.Vector.new_array(npa_to_seq(lnM_bins_knots))

    #   NCountsGauss

    ncounts_gauss = Nc.DataClusterNCountsGauss.new(cad)
    ncounts_gauss.set_size((z_bins_vec.len() - 1) * (lnM_bins_vec.len() - 1))
    ncounts_gauss.set_init(True)
    ncounts_gauss.use_norma(True)
    ncounts_gauss.set_z_obs(z_bins_vec)
    ncounts_gauss.set_lnM_obs(lnM_bins_vec)

    if fitting_Sij_type != JpasSSCType.NO_SSC:
        _set_mset_params(mset, fitting_model)
        fitting_S_ij = create_covariance_S(kernel_z, kernels_T, fitting_Sij_type, cosmo)
        ncounts_gauss.set_s_matrix(fitting_S_ij)

    if resample_Sij_type != JpasSSCType.NO_SSC:
        _set_mset_params(mset, resample_model)
        resample_S_ij = create_covariance_S(
            kernel_z, kernels_T, resample_Sij_type, cosmo
        )
        ncounts_gauss.set_resample_s_matrix(resample_S_ij)

    dset = Ncm.Dataset.new_array([ncounts_gauss])
    likelihood = Ncm.Likelihood.new(dset)

    # Extra functions

    mfunc_oa = create_mfunc_array(psml)

    # Generate experiment mock data

    rng = Ncm.RNG.seeded_new(None, resample_seed)

    _set_mset_params(mset, resample_model)
    # If resampling uses SSC, then we need to add SSC to resample
    if resample_Sij_type != JpasSSCType.NO_SSC:
        ncounts_gauss.set_has_ssc(True)

    ncounts_gauss.resample(mset, rng)

    if fitting_S_ij == JpasSSCType.NO_SSC:
        ncounts_gauss.set_has_ssc(False)
    else:
        ncounts_gauss.set_has_ssc(True)

    if use_fixed_cov:
        print("\n")
        _set_mset_params(mset, fitting_model)
        cov, updated = ncounts_gauss.compute_cov(mset)

        assert updated
        ncounts_gauss.set_cov(cov)

        ncounts_gauss.set_fix_cov(True)
        print("\n")
    # Set the model back to the resample model
    _set_mset_params(mset, resample_model)
    # Save experiment
    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    return experiment, mfunc_oa