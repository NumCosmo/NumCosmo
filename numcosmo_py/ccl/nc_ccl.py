#
# comparison.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# comparison.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""NumCosmo and CCL comparison functions."""

import numpy as np
import pyccl

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import Cosmology

# CCL uses an older release of CODATA
# The following function is a workaround to use the latest CODATA values
pyccl.physical_constants.unfreeze()
pyccl.physical_constants.KBOLTZ = Ncm.C.kb()
pyccl.physical_constants.STBOLTZ = Ncm.C.stefan_boltzmann()
pyccl.physical_constants.GNEWT = Ncm.C.G()
pyccl.physical_constants.HPLANCK = Ncm.C.h()
pyccl.physical_constants.CLIGHT = Ncm.C.c()
pyccl.physical_constants.SOLAR_MASS = Ncm.C.mass_solar()
pyccl.physical_constants.MPC_TO_METER = Ncm.C.Mpc()
pyccl.physical_constants.RHO_CRITICAL = Ncm.C.crit_mass_density_h2_solar_mass_Mpc3()
pyccl.physical_constants.freeze()


def _get_neutrino_masses(ccl_cosmo: pyccl.Cosmology) -> tuple[int, list[float]]:
    """Get neutrino masses from CCL cosmology."""
    assert isinstance(ccl_cosmo["m_nu"], (list, np.ndarray))
    massnu_length = len(ccl_cosmo["m_nu"])
    m_nu = list(ccl_cosmo["m_nu"])

    return massnu_length, m_nu


# pylint:disable-next=too-many-arguments,too-many-locals
def create_nc_obj(
    ccl_cosmo: pyccl.Cosmology,
    prec: float = 1.0e-7,
    *,
    dist_z_max: float = 15.0,
    ps_nln_z_max: float = 10.0,
    k_min: float = 1.0e-6,
    k_max: float = 1.0e3,
) -> Cosmology:
    """Create a NumCosmo object from a CCL cosmology."""
    massnu_length, m_nu = _get_neutrino_masses(ccl_cosmo)
    cosmo = Nc.HICosmoDECpl(massnu_length=massnu_length)
    cosmo.props.CCL_comp = True
    cosmo.omega_x2omega_k()
    cosmo["H0"] = ccl_cosmo["h"] * 100.0
    cosmo["Omegak"] = ccl_cosmo["Omega_k"]
    cosmo["w0"] = ccl_cosmo["w0"]
    cosmo["w1"] = ccl_cosmo["wa"]
    cosmo["Omegab"] = ccl_cosmo["Omega_b"]
    cosmo["Omegac"] = ccl_cosmo["Omega_c"]
    cosmo["ENnu"] = ccl_cosmo["N_nu_rel"]
    cosmo["Tgamma0"] = ccl_cosmo["T_CMB"]
    for i, m in enumerate(m_nu):
        cosmo[f"massnu_{i}"] = m
        cosmo[f"Tnu_{i}"] = ccl_cosmo["T_ncdm"]

    # Creates the HI Primordial object
    hiprim = Nc.HIPrimPowerLaw.new()
    hiprim["n_SA"] = ccl_cosmo["n_s"]

    # Creates the HI Reionization object
    hireion = Nc.HIReionCamb.new()

    cosmo.add_submodel(hiprim)
    cosmo.add_submodel(hireion)

    dist = Nc.Distance.new(dist_z_max)
    dist.prepare(cosmo)

    # Creating the transfer/linear power spectrum
    tf = None
    ps_ml = None

    # pylint: disable=protected-access
    if ccl_cosmo._config_init_kwargs["transfer_function"] == "eisenstein_hu":
        tf = Nc.TransferFuncEH.new()
        tf.props.CCL_comp = True

        ps_ml = Nc.PowspecMLTransfer.new(tf)
    else:
        raise ValueError(
            "Transfer function type `"
            + ccl_cosmo._config_init_kwargs["transfer_function"]  # noqa: W503
            + "` not supported"  # noqa: W503
        )
    # pylint: enable=protected-access

    ps_ml.set_kmin(k_min)
    ps_ml.set_kmax(k_max)
    ps_ml.prepare(cosmo)

    if not np.isnan(ccl_cosmo["A_s"]):
        hiprim.param_set_by_name("ln10e10ASA", np.log(1.0e10 * ccl_cosmo["A_s"]))
    else:
        A_s = np.exp(hiprim.param_get_by_name("ln10e10ASA")) * 1.0e-10
        fact = (
            ccl_cosmo["sigma8"]
            / ps_ml.sigma_tophat_R(cosmo, prec, 0.0, 8.0 / cosmo.h())  # noqa: W503
        ) ** 2
        hiprim.param_set_by_name("ln10e10ASA", np.log(1.0e10 * A_s * fact))

    ps_mln = None
    # pylint: disable=protected-access
    if ccl_cosmo._config_init_kwargs["matter_power_spectrum"] == "halofit":
        ps_mln = Nc.PowspecMNLHaloFit.new(ps_ml, ps_nln_z_max, prec)
        ps_mln.set_kmin(k_min)
        ps_mln.set_kmax(k_max)
        ps_mln.prepare(cosmo)
    # pylint: enable=protected-access

    psf = None
    if ps_ml:
        psf = Ncm.PowspecFilter.new(ps_ml, Ncm.PowspecFilterType.TOPHAT)
        psf.set_best_lnr0()

    # pylint: enable=protected-access
    return Cosmology(cosmo=cosmo, dist=dist, ps_ml=ps_ml, ps_mnl=ps_mln, psf=psf)


class CCLParams:
    """CCL cosmology parameters."""

    DEFAULT_INTEGRATION_EPSREL: float = pyccl.gsl_params.INTEGRATION_EPSREL
    DEFAULT_INTEGRATION_DISTANCE_EPSREL: float = (
        pyccl.gsl_params.INTEGRATION_DISTANCE_EPSREL
    )
    DEFAULT_INTEGRATION_LIMBER_EPSREL: float = (
        pyccl.gsl_params.INTEGRATION_LIMBER_EPSREL
    )
    DEFAULT_EPS_SCALEFAC_GROWTH: float = pyccl.gsl_params.EPS_SCALEFAC_GROWTH
    DEFAULT_ODE_GROWTH_EPSREL: float = pyccl.gsl_params.ODE_GROWTH_EPSREL
    DEFAULT_N_ITERATION: int = pyccl.gsl_params.N_ITERATION
    DEFAULT_INTEGRATION_SIGMAR_EPSREL: float = (
        pyccl.gsl_params.INTEGRATION_SIGMAR_EPSREL
    )
    DEFAULT_A_SPLINE_NLOG: int = pyccl.spline_params.A_SPLINE_NLOG
    DEFAULT_A_SPLINE_NA: int = pyccl.spline_params.A_SPLINE_NA
    DEFAULT_A_SPLINE_NA_PK: int = pyccl.spline_params.A_SPLINE_NA_PK
    DEFAULT_A_SPLINE_NLOG_PK: int = pyccl.spline_params.A_SPLINE_NLOG_PK
    DEFAULT_N_K: int = pyccl.spline_params.N_K
    DEFAULT_K_MIN: float = pyccl.spline_params.K_MIN
    DEFAULT_K_MAX: float = pyccl.spline_params.K_MAX
    DEFAULT_A_SPLINE_NLOG_SM: int = pyccl.spline_params.A_SPLINE_NLOG_SM
    DEFAULT_A_SPLINE_NA_SM: int = pyccl.spline_params.A_SPLINE_NA_SM
    DEFAULT_LOGM_SPLINE_NM: int = pyccl.spline_params.LOGM_SPLINE_NM

    @staticmethod
    def set_default_params():
        """Set CCL parameters to default values."""
        pyccl.gsl_params.INTEGRATION_EPSREL = CCLParams.DEFAULT_INTEGRATION_EPSREL
        pyccl.gsl_params.INTEGRATION_DISTANCE_EPSREL = (
            CCLParams.DEFAULT_INTEGRATION_DISTANCE_EPSREL
        )
        pyccl.gsl_params.INTEGRATION_LIMBER_EPSREL = (
            CCLParams.DEFAULT_INTEGRATION_DISTANCE_EPSREL
        )
        pyccl.gsl_params.EPS_SCALEFAC_GROWTH = CCLParams.DEFAULT_EPS_SCALEFAC_GROWTH
        pyccl.gsl_params.ODE_GROWTH_EPSREL = CCLParams.DEFAULT_ODE_GROWTH_EPSREL
        pyccl.gsl_params.N_ITERATION = CCLParams.DEFAULT_N_ITERATION
        pyccl.gsl_params.INTEGRATION_SIGMAR_EPSREL = (
            CCLParams.DEFAULT_INTEGRATION_SIGMAR_EPSREL
        )
        pyccl.spline_params.A_SPLINE_NLOG = CCLParams.DEFAULT_A_SPLINE_NLOG
        pyccl.spline_params.A_SPLINE_NA = CCLParams.DEFAULT_A_SPLINE_NA
        pyccl.spline_params.A_SPLINE_NA_PK = CCLParams.DEFAULT_A_SPLINE_NA_PK
        pyccl.spline_params.A_SPLINE_NLOG_PK = CCLParams.DEFAULT_A_SPLINE_NLOG_PK
        pyccl.spline_params.N_K = CCLParams.DEFAULT_N_K
        pyccl.spline_params.K_MIN = CCLParams.DEFAULT_K_MIN
        pyccl.spline_params.K_MAX = CCLParams.DEFAULT_K_MAX
        pyccl.spline_params.A_SPLINE_NLOG_SM = CCLParams.DEFAULT_A_SPLINE_NLOG_SM
        pyccl.spline_params.A_SPLINE_NA_SM = CCLParams.DEFAULT_A_SPLINE_NA_SM
        pyccl.spline_params.LOGM_SPLINE_NM = CCLParams.DEFAULT_LOGM_SPLINE_NM

    @staticmethod
    def set_high_prec_params():
        """Set CCL parameters to high precision values."""
        pyccl.gsl_params.INTEGRATION_EPSREL = 1.0e-7
        pyccl.gsl_params.INTEGRATION_DISTANCE_EPSREL = 1.0e-7
        pyccl.gsl_params.INTEGRATION_LIMBER_EPSREL = 1.0e-6
        pyccl.gsl_params.EPS_SCALEFAC_GROWTH = 1.0e-30
        pyccl.gsl_params.ODE_GROWTH_EPSREL = 1.0e-8
        pyccl.gsl_params.N_ITERATION = 10000
        pyccl.gsl_params.INTEGRATION_SIGMAR_EPSREL = 1.0e-9
        pyccl.spline_params.A_SPLINE_NLOG = 1000
        pyccl.spline_params.A_SPLINE_NA = 8000
        pyccl.spline_params.A_SPLINE_NA_PK = 1000
        pyccl.spline_params.A_SPLINE_NLOG_PK = 1000
        pyccl.spline_params.N_K = 1000
        pyccl.spline_params.K_MIN = 1.0e-6
        pyccl.spline_params.K_MAX = 1.0e3
        pyccl.spline_params.A_SPLINE_NLOG_SM = 100
        pyccl.spline_params.A_SPLINE_NA_SM = 100
        pyccl.spline_params.LOGM_SPLINE_NM = 300


# Missing function in CCL
def dsigmaM_dlnM(cosmo, M, a):  # pylint: disable=invalid-name
    """Compute the logarithmic derivative of the mass variance.

    Compute the logarithmic derivative of the mass variance with respect to
    the natural logarithm of the mass.
    """
    cosmo.compute_sigma()

    logM = np.log10(np.atleast_1d(M))  # pylint: disable=invalid-name
    status = 0
    # pylint: disable-next=invalid-name
    dsigMdlnM, status = pyccl.lib.dlnsigM_dlogM_vec(
        cosmo.cosmo, a, logM, len(logM), status
    )
    if np.ndim(M) == 0:
        # pylint: disable-next=invalid-name
        dsigMdlnM = dsigMdlnM[0]
    return dsigMdlnM
