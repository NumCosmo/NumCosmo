#
# jpas_forecast24.py
#
# Mon Feb 20 22:31:10 2024
# Copyright 2024 Sandro Dias Pinto Vitenti
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""Factory functions for J-Pas 2024 cluster abundance forecasting.

This module provides functions to construct a complete cluster number counts
forecast experiment using the NumCosmo and PySSC libraries, tailored for
the J-PAS survey.

It includes:
1. Cosmological model setup (flat wCDM).
2. Astrophysical nuisance models (mass-observable relation and photoz scatter).
3. Binning, kernels, and HEALPix masks for the J-PAS sky area.
4. Calculation and integration of Super Sample Covariance (SSC).
5. Generation of mock data for forecast analysis.
"""

from typing import cast
from enum import StrEnum, auto
import time
import numpy as np

from numcosmo_py import Ncm, Nc
from numcosmo_py.external.pyssc import pyssc as PySSC


class JpasSSCType(StrEnum):
    """J-Pas 2024 Super Sample Covariance (SSC) types.

    Defines the types of sky cuts and corresponding SSC matrices to use for
    the covariance calculation.
    """

    NO_SSC = auto()
    """No Super Sample Covariance is included in the covariance matrix."""
    FULLSKY = auto()
    """Uses a full-sky approximation for the SSC matrix $S_{ij}$."""
    FULL = auto()
    """Uses the HEALPix mask for the full J-PAS sky coverage for SSC."""
    GUARANTEED = auto()
    """Uses the HEALPix mask for the guaranteed J-PAS sky coverage for SSC."""


class ClusterMassType(StrEnum):
    """Mass-observable relation types.

    Defines the model used to relate the cluster's observable property (e.g., richness,
    $N_{gal}$) to its underlying true mass ($lnM$ or $ln M_{200c}$).
    """

    NODIST = auto()
    """No distribution

    assumes a perfect mass-observable relation (no intrinsic scatter).
    """
    ASCASO = auto()
    """Uses the mass-richness relation parameters from Ascaso et al. (2015)."""


class ClusterRedshiftType(StrEnum):
    """Photoz types.

    Defines the model for the cluster's photometric redshift (photoz) distribution.
    """

    NODIST = auto()
    """No distribution; assumes perfect spectroscopic redshift (no photoz error)."""
    GAUSS = auto()
    """Uses a Gaussian distribution for the photometric redshift uncertainty."""


def create_zbins_kernels(
    z_min: float = 0.1,
    z_max: float = 0.8,
    nknots: int = 8,
    kernel_nknots: int = 400,
    kernel_zmax: float = 1.9,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Create the redshift bins and the SSC projection kernels.

    The kernels are top-hat functions in redshift used to project the 3D power
    spectrum when computing the Super Sample Covariance (SSC) matrix $S_{ij}$.

    :param z_min: Minimum redshift of the cluster bin range.
    :param z_max: Maximum redshift of the cluster bin range.
    :param nknots: Number of knots (boundaries) for the redshift bins (e.g., 8 knots
        give 7 bins).
    :param kernel_nknots: Number of points for the kernel's internal redshift array.
    :param kernel_zmax: Maximum redshift for the kernel evaluation range.
    :return: A tuple containing:
        - kernel_z (np.ndarray): Redshift values for kernel evaluation ($z$).
        - kernels_T (np.ndarray): The kernel matrix $T_i(z)$ where $T_i(z) = 1/Delta
          z_i$ inside bin $i$ and 0 otherwise.
        - z_bins_knots (np.ndarray): Redshift bin boundaries.
    """
    # Redshift bins and knots (boundaries)
    z_bins_len = nknots - 1
    z_bins_knots = np.linspace(z_min, z_max, num=nknots)

    # Kernel redshift axis (used for numerical integration of the power spectrum)
    kernel_z = np.linspace(0.0, kernel_zmax, num=kernel_nknots + 1)[1:]

    # Initialize the kernel matrix (number of bins x number of kernel points)
    kernels_T = np.zeros((z_bins_len, kernel_nknots))

    # Create the top-hat kernel for each redshift bin
    for i, (zminbin, zmaxbin) in enumerate(zip(z_bins_knots[:-1], z_bins_knots[1:])):
        Dz = zmaxbin - zminbin  # Bin width

        kernel = np.zeros_like(kernel_z)
        # Check which kernel points fall into the current bin
        kernel[(kernel_z >= zminbin) & (kernel_z <= zmaxbin)] = 1.0
        # Normalize the kernel by the bin width Dz
        kernels_T[i] = kernel / Dz

    return kernel_z, kernels_T, z_bins_knots


def create_lnM_obs_bins(
    lnM_obs_min: float = np.log(10.0) * 14.0,
    lnM_obs_max: float = np.log(10.0) * 15.0,
    nknots: int = 2,
) -> np.ndarray:
    """Create the log-observable mass bins for the J-Pas 2024 forecast.

    The mass variable is typically $ln M$ (log-base-e of $M_{200c}$) in units of
    $h^{-1} M_\\odot$. These bins are defined in terms of the observed proxy (e.g.,
    richness).

    :param lnM_obs_min: Minimum log-mass-observable knot boundary.
    :param lnM_obs_max: Maximum log-mass-observable knot boundary.
    :param nknots: Number of knots (boundaries) for the log-mass bins.
    :return: An array of log-mass bin boundaries (knots).
    """
    lnM_bins_obs_knots = np.linspace(lnM_obs_min, lnM_obs_max, nknots)

    return lnM_bins_obs_knots


def survey_area(sky_cut: JpasSSCType) -> float:
    """Survey area in square degrees ($\text{sqd}$) for the J-Pas 2024 forecast.

    :param sky_cut: The type of sky coverage assumed.
    :raises ValueError: If the sky cut type is not FULLSKY or a masked type.
    :return: The survey area in sq. degrees.
    """
    match sky_cut:
        case JpasSSCType.FULLSKY:
            # Area in sq. degrees for full-sky approximation
            survey_area0 = 20009.97
        case JpasSSCType.FULL:
            # Area in sq. degrees for full J-PAS survey region
            survey_area0 = 8500.0  # Fictional value, replace with actual
        case JpasSSCType.GUARANTEED:
            # Area in sq. degrees for guaranteed J-PAS survey region
            survey_area0 = 2959.1
        case _:
            raise ValueError(f"Invalid sky cut type for area calculation: {sky_cut}")

    return survey_area0


def create_mask_guaranteed(nside: int = 512) -> np.ndarray:
    """Create the HEALPix mask for the J-Pas 2024 guaranteed sky coverage.

    The mask is defined in ecliptic coordinates.

    :param nside: The HEALPix resolution parameter.
    :return: A numpy array representing the pixel mask (1=in survey, 0=out).
    """
    smap = Ncm.SphereMap.new(nside)
    npix = smap.get_npix()
    angles = np.array([smap.pix2ang_ring(i) for i in range(npix)])

    mask1_guaranteed = np.zeros(npix)
    mask2_guaranteed = np.zeros(npix)
    mask3_guaranteed = np.zeros(npix)

    # HEALPix coordinates (theta, phi) in ecliptic system
    pix_theta_ecl = angles[:, 0]
    pix_phi_ecl = angles[:, 1]

    # Define the three distinct guaranteed regions based on ecliptic coordinates
    # These definitions map specific sky areas to the mask.
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

    # Combine the masks
    mask_guaranteed = mask1_guaranteed + mask2_guaranteed + mask3_guaranteed

    return mask_guaranteed


def create_mask_full(nside: int = 512) -> np.ndarray:
    """Create the HEALPix mask for the J-Pas 2024 full sky coverage.

    The mask is defined in ecliptic coordinates.

    :param nside: The HEALPix resolution parameter.
    :return: A numpy array representing the pixel mask (1=in survey, 0=out).
    """
    smap = Ncm.SphereMap.new(nside)
    npix = smap.get_npix()
    angles = np.array([smap.pix2ang_ring(i) for i in range(npix)])
    pix_theta_ecl = angles[:, 0]
    pix_phi_ecl = angles[:, 1]

    # Total mask initialization
    mask1_full = np.zeros(npix)
    mask2_full = np.zeros(npix)
    mask3_full = np.zeros(npix)

    # Define the three distinct full regions based on ecliptic coordinates
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

    # Combine the masks
    mask_full = mask1_full + mask2_full + mask3_full

    return mask_full


def create_cosmo() -> Nc.HICosmo:
    """Create a fiducial flat wCDM cosmology model for J-Pas 2024 forecast.

    Uses an $text{Nc.HICosmoDEXcdm}$ model (Flat $Lambda$CDM extension).
    Sets fiducial values and specifies which parameters are free to be fitted.

    :return: An initialized NumCosmo HICosmoDEXcdm model.
    """
    cosmo = Nc.HICosmoDEXcdm()

    # Set default fitting types (typically linear scale, fixed tolerance)
    cosmo.params_set_default_ftype()
    # Ensure $Omega_x$ (Dark Energy) is calculated from $Omega_k$ (curvature) and
    # $Omega_m$
    cosmo.omega_x2omega_k()

    # --- Set Fiducial Values (from current best-fit cosmology) ---
    cosmo["H0"] = 67.81
    cosmo["Omegab"] = 0.0486
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.00  # Flat cosmology

    # --- Define Fitting Status ---
    cosmo.param_set_desc("H0", {"fit": False})
    # $Omega_c$ is a free parameter
    cosmo.param_set_desc(
        "Omegac",
        {
            "lower-bound": 0.1,
            "upper-bound": 0.9,
            "scale": 1.0e-2,
            "abstol": 1.0e-50,
            "fit": True,
            "value": 0.2612,
        },
    )
    cosmo.param_set_desc("Omegab", {"fit": False})
    # $w$ is a free parameter (Dark Energy EoS)
    cosmo.param_set_desc("w", {"fit": True})
    cosmo.param_set_desc("Omegak", {"fit": False})

    # --- Primordial Power Spectrum Model ---
    prim = Nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745  # Amplitude equivalent to $A_s$
    prim["n_SA"] = 0.9660  # Scalar spectral index

    prim.param_set_desc(
        "ln10e10ASA", {"fit": True}
    )  # $A_s$ amplitude is a fit parameter
    prim.param_set_desc("n_SA", {"fit": False})

    # --- Reionization Model (fixed) ---
    reion = Nc.HIReionCamb.new()
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def create_mfunc_array(psml: Nc.PowspecML) -> Ncm.ObjArray:
    """Create a list of extra functions (derived parameters) for J-Pas 2024 forecast.

    These functions compute parameters like $sigma_8$ and $Omega_{m0}$ from the primary
    cosmological parameters, often required for visualization or as nuisance
    parameters.

    :param psml: The power spectrum model to be used for calculations.
    :return: An array of Ncm.MSetFuncList objects (derived parameters).
    """
    mfunc_oa = Ncm.ObjArray.new()

    # Set k-range for power spectrum calculation
    psml.require_kmin(1.0e-5)
    psml.require_kmax(1.0e1)
    psml.require_zi(0.0)  # Initial redshift
    psml.require_zf(1.0)  # Final redshift (for some internal calculations)

    # Setup filter for $sigma_8$ calculation (Top-Hat in real space)
    psf_cbe = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf_cbe.set_best_lnr0()  # Set $R=8 text{ Mpc}/h$ for $sigma_8$

    # $sigma_8$ function
    mfunc_sigma8 = Ncm.MSetFuncList.new("NcHICosmo:sigma8", psf_cbe)
    mfunc_oa.add(mfunc_sigma8)

    # $Omega_{m0}$ function (current total matter density parameter)
    mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
    mfunc_oa.add(mfunc_Omegam)

    return mfunc_oa


def create_covariance_S_fullsky(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base SSC covariance matrix $S_{ij}$ based on the full sky (isotropic).

    Uses the PySSC function $text{Sij}$, which does not require a mask.

    :param kernel_z: Redshift values for kernel evaluation.
    :param kernels_T: The kernel matrix $T_i(z)$.
    :param cosmo: The fiducial cosmology model.
    :return: The $S_{ij}$ matrix as a NumCosmo.Matrix object.
    """
    S_fullsky_array = PySSC.Sij(kernel_z, kernels_T, cosmo)

    S_fullsky = Ncm.Matrix.new_array(
        S_fullsky_array.flatten(), S_fullsky_array.shape[1]
    )

    return S_fullsky


def create_covariance_S_guaranteed(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base SSC covariance matrix $S_{ij}$ for the guaranteed mask.

    Uses the PySSC function $text{Sij_psky}$ with the generated HEALPix mask.

    :param kernel_z: Redshift values for kernel evaluation.
    :param kernels_T: The kernel matrix $T_i(z)$.
    :param cosmo: The fiducial cosmology model.
    :return: The $S_{ij}$ matrix (NumCosmo.Matrix).
    """
    mask = create_mask_guaranteed()

    S_guaranteed_array = PySSC.Sij_psky(kernel_z, kernels_T, cosmo, mask=mask)

    S_guaranteed = Ncm.Matrix.new_array(
        S_guaranteed_array.flatten(), S_guaranteed_array.shape[1]
    )

    return S_guaranteed


def create_covariance_S_full(
    kernel_z: np.ndarray, kernels_T: np.ndarray, cosmo: Nc.HICosmo
) -> Ncm.Matrix:
    """Create the base SSC covariance matrix $S_{ij}$ for the full mask.

    Uses the PySSC function $text{Sij_psky}$ with the generated HEALPix mask.

    :param kernel_z: Redshift values for kernel evaluation.
    :param kernels_T: The kernel matrix $T_i(z)$.
    :param cosmo: The fiducial cosmology model.
    :return: The $S_{ij}$ matrix (NumCosmo.Matrix).
    """
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
    """Create the base SSC covariance matrix $S_{ij}$ based on the sky cut type.

    This acts as a router to the appropriate $S_{ij}$ calculation function.

    :param kernel_z: Redshift values for kernel evaluation.
    :param kernels_T: The kernel matrix $T_i(z)$.
    :param sky_cut: The type of sky coverage to use for $S_{ij}$ calculation.
    :param cosmo: The fiducial cosmology model.
    :raises ValueError: If the sky cut type is invalid (e.g., NO_SSC is passed).
    :return: The $S_{ij}$ matrix (NumCosmo.Matrix).
    """
    if sky_cut == JpasSSCType.FULLSKY:
        S = create_covariance_S_fullsky(kernel_z, kernels_T, cosmo)
    elif sky_cut == JpasSSCType.FULL:
        S = create_covariance_S_full(kernel_z, kernels_T, cosmo)
    elif sky_cut == JpasSSCType.GUARANTEED:
        S = create_covariance_S_guaranteed(kernel_z, kernels_T, cosmo)
    else:
        raise ValueError(f"Invalid sky cut type for Sij calculation: {sky_cut}")

    return S


def create_cluster_mass(
    cluster_mass_type: ClusterMassType,
) -> tuple[Nc.ClusterMass, float, float]:
    """Create the cluster mass-observable relation model.

    :param cluster_mass_type: The type of mass-observable relation to use.
    :raises ValueError: If the cluster mass type is invalid.
    :return: A tuple containing:
        - cluster_m: An initialized NumCosmo ClusterMass model.
        - lnM_obs_min: Recommended minimum log-observable for this model.
        - lnM_obs_max: Recommended maximum log-observable for this model.
    """
    cluster_m: Nc.ClusterMassAscaso | Nc.ClusterMassNodist
    if cluster_mass_type == ClusterMassType.NODIST:
        # Assumes perfect relation: $M_{obs} = M_{true}$
        lnM_obs_min = np.log(10) * 14.0
        lnM_obs_max = np.log(10) * 16.0
        cluster_m = Nc.ClusterMassNodist(lnM_min=lnM_obs_min, lnM_max=lnM_obs_max)
    elif cluster_mass_type == ClusterMassType.ASCASO:
        # Ascaso et al. (2015) relation for richness
        lnM_obs_min = np.log(5.0)  # ln(richness_min) = ln(5)
        lnM_obs_max = np.log(10) * 2.5  # ln(richness_max) = ln(10^2.5)
        cluster_m = Nc.ClusterMassAscaso(
            lnRichness_min=lnM_obs_min, lnRichness_max=lnM_obs_max, z0=0
        )
        # Set fiducial values for the mass-richness relation parameters
        cluster_m["mup0"] = 3.207  # $mu_{M|lambda}$ normalization
        cluster_m["mup1"] = 0.993  # $mu_{M|lambda}$ mass dependence
        cluster_m["mup2"] = 0  # $mu_{M|lambda}$ redshift dependence
        cluster_m["sigmap0"] = 0.456  # $sigma_{M|lambda}$ normalization
        cluster_m["sigmap1"] = -0.169  # $sigma_{M|lambda}$ mass dependence
        cluster_m["sigmap2"] = 0  # $sigma_{M|lambda}$ redshift dependence
    else:
        raise ValueError(f"Invalid cluster mass type: {cluster_mass_type}")

    return cluster_m, lnM_obs_min, lnM_obs_max


def create_cluster_redshift(
    cluster_redshift_type: ClusterRedshiftType,
) -> Nc.ClusterRedshift:
    """Create the cluster photometric redshift (photoz) relation model.

    :param cluster_redshift_type: The type of photoz relation to use.
    :param cluster_redshift_type: The type of photoz relation to use.
    :raises ValueError: If the cluster redshift type is invalid.
    :return: An initialized NumCosmo ClusterRedshift model.
    """
    cluster_z: Nc.ClusterRedshiftNodist | Nc.ClusterPhotozGaussGlobal
    if cluster_redshift_type == ClusterRedshiftType.NODIST:
        # Assumes perfect redshift: $z_{obs} = z_{true}$
        cluster_z = Nc.ClusterRedshiftNodist(z_min=0.0, z_max=2.0)
    elif cluster_redshift_type == ClusterRedshiftType.GAUSS:
        # Gaussian photoz error distribution
        cluster_z = Nc.ClusterPhotozGaussGlobal(pz_min=0.0, pz_max=1.0)
        cluster_z["z-bias"] = 0  # Photoz bias $langle z_{obs} - z_{true} rangle$
        cluster_z["sigma0"] = 0.1  # Photoz scatter $sigma_z$
    else:
        raise ValueError(f"Invalid cluster redshift type: {cluster_redshift_type}")

    return cluster_z


def _set_mset_params(mset: Ncm.MSet, params: tuple[float, float, float]) -> None:
    """Set the parameters for the cosmology model ($Omega_c, w, sigma_8$).

    It calculates the primordial amplitude $ln(10^{10} A_s)$ required to match
    the target $sigma_8$ for the given $Omega_c$ and $w$. This ensures $sigma_8$ is
    always consistent with the input parameters.

    :param mset: The model set containing $text{NcHICosmo}$ and $text{NcHIPrim}$ models.
    :param params: The parameters for the cosmology model, $(Omega_c, w, sigma_8)$.
    """
    # Setup for internal $sigma_8$ calculation
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)
    Omegac, w, sigma8 = params

    # Retrieve submodels and use 'cast' for type hints (NumCosmo helper)
    cosmo: Nc.HICosmo = cast(Nc.HICosmo, mset["NcHICosmo"])
    prim: Nc.HIPrim = cast(Nc.HIPrim, mset["NcHIPrim"])

    # Set $Omega_c$ and $w$ directly
    mset["NcHICosmo"]["Omegac"] = Omegac
    mset["NcHICosmo"]["w"] = w

    # Calculate and set the new $A_s$ value to match the target $sigma_8$
    A_s_fid = np.exp(prim["ln10e10ASA"]) * 1.0e-10
    # Calculate $sigma_8$ for a normalized $A_s = 10^{-10}$
    sigma8_fid_norm = psml.sigma_tophat_R(cosmo, 1.0e-7, 0.0, 8.0 / cosmo.h())
    # Scaling factor: $(text{target } sigma_8 / text{fiducial } sigma_8)^2$
    fact = (sigma8 / sigma8_fid_norm) ** 2
    # Set the new $A_s$ value
    prim["ln10e10ASA"] = np.log(1.0e10 * A_s_fid * fact)


def generate_jpas_forecast_2024(
    area: float = 2959.1,
    z_min: float = 0.1,
    z_max: float = 0.8,
    znknots: int = 8,
    cluster_redshift_type: ClusterRedshiftType = ClusterRedshiftType.NODIST,
    lnM_obs_min: float | None = None,
    lnM_obs_max: float | None = None,
    lnMobsnknots: int = 2,
    cluster_mass_type: ClusterMassType = ClusterMassType.NODIST,
    use_fixed_cov: bool = False,
    fitting_model: tuple[float, float, float] = (0.2612, -1.0, 0.8159),
    resample_model: tuple[float, float, float] = (0.2612, -1.0, 0.8159),
    resample_seed: int = 1234,
    fitting_Sij_type: JpasSSCType = JpasSSCType.FULLSKY,
    resample_Sij_type: JpasSSCType = JpasSSCType.NO_SSC,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate J-Pas 2024 cluster abundance forecast experiment dictionary.

    This function assembles all model components, sets up the likelihood, generates
    mock data, and applies covariance matrices (including SSC) for a Fisher matrix
    forecast or likelihood analysis.

    :param area: The fixed survey area in square degrees to use. This is overridden if
        a masked SSC type (FULL/GUARANTEED) is used for resampling or fitting.
    :param z_min: Minimum redshift for the cluster bins.
    :param z_max: Maximum redshift for the cluster bins.
    :param znknots: Number of redshift bin boundaries.
    :param cluster_redshift_type: Model for cluster photometric redshift scatter.
    :param lnM_obs_min: Minimum log-observable for the bins. If None, uses the default
        from the cluster_mass_type model.
    :param lnM_obs_max: Maximum log-observable for the bins. If None, uses the default
        from the cluster_mass_type model.
    :param lnMobsnknots: Number of log-observable bin boundaries.
    :param cluster_mass_type: Model for mass-observable relation scatter/bias.
    :param use_fixed_cov: If True, the full covariance is computed once and fixed.
    :param fitting_model: $(Omega_c, w, sigma_8)$ for the model used to calculate the
        **theoretical covariance matrix**.
    :param resample_model: $(Omega_c, w, sigma_8)$ for the model used to **generate the
        mock data vector**.
    :param resample_seed: Seed for the random number generator used for mock data.
    :param fitting_Sij_type: The type of SSC matrix $S_{ij}$ to use for the fitting
        (theoretical) covariance.
    :param resample_Sij_type: The type of SSC matrix $S_{ij}$ to use for generating the
        mock data vector.
    :return: A tuple containing:
             - experiment (Ncm.ObjDictStr): Dictionary containing the likelihood and
               model set.
             - mfunc_oa (Ncm.ObjArray): Array of extra (derived) functions.
    """
    t_start = time.time()
    print("Starting J-PAS forecast generation...")

    # Adjust area if a masked SSC type is used, giving priority to resample model
    if resample_Sij_type in (JpasSSCType.FULL, JpasSSCType.GUARANTEED):
        area = survey_area(resample_Sij_type)
        print(
            f"sky cut is of type {resample_Sij_type} adjusting to "
            f"corresponding survey area = {area:.2f} sqd. "
        )
    elif fitting_Sij_type in (JpasSSCType.FULL, JpasSSCType.GUARANTEED):
        area = survey_area(fitting_Sij_type)
        print(
            f"sky cut is of type {fitting_Sij_type} adjusting to "
            f"corresponding survey area = {area:.2f} sqd. "
        )

    # --- Setup Core NumCosmo/PySSC Tools ---
    t_setup_start = time.time()
    dist = Nc.Distance.new(2.0)
    tf = Nc.TransferFuncEH()
    psml = Nc.PowspecMLTransfer.new(tf)
    psml.require_kmin(1.0e-6)
    psml.require_kmax(1.0e3)
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()  # Sets R=8 Mpc/h for $sigma_8$
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    mulf.set_Delta(200.0)  # Uses $Delta = 200c$ for mass definition
    hmf = Nc.HaloMassFunction.new(dist, psf, mulf)
    hbias_Tinker = Nc.HaloBiasTinker.new(hmf)
    cad = Nc.ClusterAbundance.new(hmf, hbias_Tinker)
    # Set area in steradians
    cad.set_area(area * (np.pi / 180) ** 2)

    # --- Instantiate Models ---
    cluster_m, default_lnM_obs_min, default_lnM_obs_max = create_cluster_mass(
        cluster_mass_type
    )
    cluster_z = create_cluster_redshift(cluster_redshift_type)
    cosmo = create_cosmo()

    # Use defaults if not specified
    if lnM_obs_min is None:
        lnM_obs_min = default_lnM_obs_min
    if lnM_obs_max is None:
        lnM_obs_max = default_lnM_obs_max

    mset = Ncm.MSet.new_array([cosmo, cluster_m, cluster_z])
    mset.prepare_fparam_map()  # Map fitting parameters for efficient access
    print(f"Model setup completed in {time.time() - t_setup_start:.2f}s")

    # --- Define Bins and Kernels ---
    kernel_z, kernels_T, z_bins_knots = create_zbins_kernels(
        z_min=z_min, z_max=z_max, nknots=znknots
    )
    lnM_bins_obs_knots = create_lnM_obs_bins(
        lnM_obs_min=lnM_obs_min, lnM_obs_max=lnM_obs_max, nknots=lnMobsnknots
    )
    z_bins_vec = Ncm.Vector.new_array(z_bins_knots)
    lnM_obs_bins_vec = Ncm.Vector.new_array(lnM_bins_obs_knots)

    # --- Likelihood: Cluster NCountsGauss ---
    ncounts_gauss = Nc.DataClusterNCountsGauss.new(cad)
    # Total number of bins = (N_z - 1) * (N_lnM - 1)
    n_bins = (z_bins_vec.len() - 1) * (lnM_obs_bins_vec.len() - 1)
    ncounts_gauss.set_size(n_bins)
    ncounts_gauss.set_init(True)
    ncounts_gauss.use_norma(True)  # Use normalized counts/volume
    ncounts_gauss.set_z_obs(z_bins_vec)
    ncounts_gauss.set_lnM_obs(lnM_obs_bins_vec)

    # Set fitting (theoretical) SSC matrix
    if fitting_Sij_type != JpasSSCType.NO_SSC:
        print(f"Computing fitting SSC matrix ({fitting_Sij_type})...", flush=True)
        t_fit_ssc_start = time.time()
        _set_mset_params(mset, fitting_model)
        fitting_S_ij = create_covariance_S(kernel_z, kernels_T, fitting_Sij_type, cosmo)
        ncounts_gauss.set_s_matrix(fitting_S_ij)
        print(f"Fitting SSC matrix computed in {time.time() - t_fit_ssc_start:.2f}s")

    # Set resampling (mock data) SSC matrix
    if resample_Sij_type != JpasSSCType.NO_SSC:
        print(f"Computing resampling SSC matrix ({resample_Sij_type})...", flush=True)
        t_resamp_ssc_start = time.time()
        _set_mset_params(mset, resample_model)
        resample_S_ij = create_covariance_S(
            kernel_z, kernels_T, resample_Sij_type, cosmo
        )
        ncounts_gauss.set_resample_s_matrix(resample_S_ij)
        print(
            f"Resampling SSC matrix computed in {time.time() - t_resamp_ssc_start:.2f}s"
        )

    dset = Ncm.Dataset.new_array([ncounts_gauss])
    likelihood = Ncm.Likelihood.new(dset)

    # --- Extra functions ---
    mfunc_oa = create_mfunc_array(psml)

    # --- Generate experiment mock data ---
    rng = Ncm.RNG.seeded_new(None, resample_seed)

    _set_mset_params(mset, resample_model)

    # Temporarily enable SSC for resampling if needed
    if resample_Sij_type != JpasSSCType.NO_SSC:
        ncounts_gauss.set_has_ssc(True)

    print("Generating mock data...", flush=True)
    t_resample_start = time.time()
    # Generate the mock data vector (this computes $mu_{resample}$ and resamples $N$)
    ncounts_gauss.resample(mset, rng)
    print(f"Mock data generated in {time.time() - t_resample_start:.2f}s")

    # Reset $text{has_ssc}$ based on the $text{fitting_Sij_type}$ for the actual
    # likelihood evaluation
    if fitting_Sij_type == JpasSSCType.NO_SSC:
        ncounts_gauss.set_has_ssc(False)
    else:
        ncounts_gauss.set_has_ssc(True)

    # Compute and fix the full covariance matrix (if requested, useful for Fisher
    # matrix)
    if use_fixed_cov:
        print("Computing fixed covariance...", flush=True)
        t_cov_start = time.time()
        _set_mset_params(mset, fitting_model)  # Use fitting model for fixed covariance
        cov, updated = ncounts_gauss.compute_cov(mset)

        assert updated
        ncounts_gauss.set_cov(cov)
        ncounts_gauss.set_fix_cov(True)
        print(f"Fixed covariance computed in {time.time() - t_cov_start:.2f}s")

    # Set the model back to the resample model (or fiducial) for initial fitting state
    _set_mset_params(mset, resample_model)

    # --- Save experiment ---
    experiment = Ncm.ObjDictStr()
    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    print(f"\nTotal experiment generation time: {time.time() - t_start:.2f}s")
    return experiment, mfunc_oa
