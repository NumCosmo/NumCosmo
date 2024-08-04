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

import math
import numpy as np
import pyccl

from numcosmo_py import Ncm, Nc

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


# pylint:disable-next=too-many-arguments,too-many-locals
def create_nc_obj(
    ccl_cosmo,
    prec=1.0e-7,
    dist_z_max=15.0,
    ps_nln_z_max=10.0,
    k_min=1.0e-6,
    k_max=1.0e3,
):
    """Create a NumCosmo object from a CCL cosmology."""
    cosmo = Nc.HICosmoDECpl(massnu_length=0)
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", ccl_cosmo["h"] * 100)
    cosmo.param_set_by_name("Omegak", ccl_cosmo["Omega_k"])
    cosmo.param_set_by_name("w0", ccl_cosmo["w0"])
    cosmo.param_set_by_name("w1", ccl_cosmo["wa"])
    cosmo.param_set_by_name("Omegab", ccl_cosmo["Omega_b"])
    cosmo.param_set_by_name("Omegac", ccl_cosmo["Omega_c"])
    cosmo.param_set_by_name("ENnu", ccl_cosmo["Neff"])
    cosmo.param_set_by_name("Tgamma0", ccl_cosmo["T_CMB"])

    # Creates the HI Primordial object
    hiprim = Nc.HIPrimPowerLaw.new()
    hiprim.param_set_by_name("n_SA", ccl_cosmo["n_s"])

    cosmo.add_submodel(hiprim)

    dist = Nc.Distance.new(dist_z_max)
    dist.prepare(cosmo)

    # Checking if neutrinos are compatible
    if isinstance(ccl_cosmo["m_nu"], (list, np.ndarray)):
        for m_nu_i in ccl_cosmo["m_nu"]:
            if m_nu_i != 0.0:
                raise ValueError("Massive neutrinos are not supported")
    else:
        if ccl_cosmo["m_nu"] != 0:
            raise ValueError("Massive neutrinos are not supported")

    # Creating the transfer/linear power spectrum
    tf = None  # pylint: disable=invalid-name
    ps_lin = None

    # pylint: disable=protected-access
    if ccl_cosmo._config_init_kwargs["transfer_function"] == "eisenstein_hu":
        tf = Nc.TransferFuncEH.new()  # pylint: disable=invalid-name
        tf.props.CCL_comp = True

        ps_lin = Nc.PowspecMLTransfer.new(tf)
    else:
        raise ValueError(
            "Transfer function type `"
            + ccl_cosmo._config_init_kwargs["transfer_function"]  # noqa: W503
            + "` not supported"  # noqa: W503
        )

    ps_lin.set_kmin(k_min)
    ps_lin.set_kmax(k_max)
    ps_lin.prepare(cosmo)

    if not math.isnan(ccl_cosmo["A_s"]):
        hiprim.param_set_by_name("ln10e10ASA", math.log(1.0e10 * ccl_cosmo["A_s"]))
    else:
        A_s = math.exp(hiprim.param_get_by_name("ln10e10ASA")) * 1.0e-10
        fact = (
            ccl_cosmo["sigma8"]
            / ps_lin.sigma_tophat_R(cosmo, prec, 0.0, 8.0 / cosmo.h())  # noqa: W503
        ) ** 2
        hiprim.param_set_by_name("ln10e10ASA", math.log(1.0e10 * A_s * fact))

    ps_nln = None
    if ccl_cosmo._config_init_kwargs["matter_power_spectrum"] == "halofit":
        ps_nln = Nc.PowspecMNLHaloFit.new(ps_lin, ps_nln_z_max, prec)
        ps_nln.set_kmin(k_min)
        ps_nln.set_kmax(k_max)
        ps_nln.prepare(cosmo)

    hmfunc = None

    if ps_lin:
        psf = Ncm.PowspecFilter.new(ps_lin, Ncm.PowspecFilterType.TOPHAT)
        psf.set_best_lnr0()

    # pylint: enable=protected-access
    return cosmo, dist, ps_lin, ps_nln, hmfunc


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
        pyccl.gsl_params.INTEGRATION_EPSREL = 1.0e-13
        pyccl.gsl_params.INTEGRATION_DISTANCE_EPSREL = 1.0e-13
        pyccl.gsl_params.INTEGRATION_LIMBER_EPSREL = 1.0e-9
        pyccl.gsl_params.EPS_SCALEFAC_GROWTH = 1.0e-30
        pyccl.gsl_params.ODE_GROWTH_EPSREL = 1.0e-13
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
