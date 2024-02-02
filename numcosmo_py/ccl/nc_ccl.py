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
        # pylint: disable-next=invalid-name
        A_s = math.exp(hiprim.param_get_by_name("ln10e10ASA")) * 1.0e-10
        fact = (
            ccl_cosmo["sigma8"]
            / ps_lin.sigma_tophat_R(cosmo, prec, 0.0, 8.0 / cosmo.h())  # noqa: W503
        ) ** 2
        hiprim.param_set_by_name("ln10e10ASA", math.log(1.0e10 * A_s * fact))

    ps_nln = None
    if ccl_cosmo._config_init_kwargs["matter_power_spectrum"] == "halofit":
        ps_nln = Nc.PowspecMNLHaloFit.new(ps_lin, ps_nln_z_max, prec)
        ps_nln.prepare(cosmo)

    hmfunc = None

    if ps_lin:
        psf = Ncm.PowspecFilter.new(ps_lin, Ncm.PowspecFilterType.TOPHAT)
        psf.set_best_lnr0()

        #if ccl_cosmo._config_init_kwargs["mass_function"] == "tinker10":
        #    # pylint: disable-next=invalid-name
        #    hmf_T10 = Nc.MultiplicityFuncTinkerMeanNormalized.new()
        #    hmfunc = Nc.HaloMassFunction.new(dist, psf, hmf_T10)
        #    hmfunc.prepare(cosmo)

    # pylint: enable=protected-access
    return cosmo, dist, ps_lin, ps_nln, hmfunc


def ccl_cosmo_set_high_prec():
    """Set CCL cosmology to high precision."""
    pyccl.gsl_params.INTEGRATION_EPSREL = 1.0e-13
    pyccl.gsl_params.ODE_GROWTH_EPSREL = 1.0e-13
    pyccl.gsl_params.N_ITERATION = 10000
    pyccl.gsl_params.INTEGRATION_SIGMAR_EPSREL = 1.0e-9
    pyccl.spline_params.A_SPLINE_NLOG = 1000
    pyccl.spline_params.A_SPLINE_NA = 1000
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
    """Derivative of the mass variance with respect to the logarithm of the
    mass."""

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
