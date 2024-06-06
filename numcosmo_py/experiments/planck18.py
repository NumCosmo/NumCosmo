#
# planck18.py
#
# Mon Feb 20 22:31:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# planck18.py
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

"""Factory functions to generate Planck18 likelihood and models."""

from enum import Enum

import numpy as np

from numcosmo_py import Ncm, Nc


class Planck18Types(str, Enum):
    """Planck 18 baseline data combinations."""

    TT = "TT"
    TTTEEE = "TTTEEE"


class Planck18HIPrimModel(str, Enum):
    """Planck 18 primordial model."""

    POWER_LAW = "power-law"
    TWO_FLUIDS = "two-fluids"


EXP_PARAMETERS: dict[tuple[str, str], dict[str, float]] = {
    (Planck18Types.TT, Planck18HIPrimModel.POWER_LAW): {
        "NcHICosmo:H0": 66.86,
        "NcHICosmo:omegac": 0.12068,
        "NcHICosmo:omegab": 0.022126,
        "NcHIPrim:ln10e10ASA": 3.0413,
        "NcHIPrim:n_SA": 0.9635,
        "NcHIReion:z_re": 7.54,
        "NcPlanckFI:A_cib_217": 48.5,
        "NcPlanckFI:xi_sz_cib": 0.32,
        "NcPlanckFI:A_sz": 7.03,
        "NcPlanckFI:ps_A_100_100": 254.9,
        "NcPlanckFI:ps_A_143_143": 49.8,
        "NcPlanckFI:ps_A_143_217": 47.3,
        "NcPlanckFI:ps_A_217_217": 119.9,
        "NcPlanckFI:ksz_norm": 0.00,
        "NcPlanckFI:gal545_A_100": 8.86,
        "NcPlanckFI:gal545_A_143": 10.80,
        "NcPlanckFI:gal545_A_143_217": 19.43,
        "NcPlanckFI:gal545_A_217": 94.8,
        "NcPlanckFI:calib_100T": 0.99965,
        "NcPlanckFI:calib_217T": 0.99825,
        "NcPlanckFI:A_planck": 1.00046,
    },
    (Planck18Types.TTTEEE, Planck18HIPrimModel.POWER_LAW): {
        "NcHICosmo:H0": 67.32,
        "NcHICosmo:omegac": 0.12010,
        "NcHICosmo:omegab": 0.022377,
        "NcHIPrim:ln10e10ASA": 3.0447,
        "NcHIPrim:n_SA": 0.96589,
        "NcHIReion:z_re": 7.68,
        "NcPlanckFI:A_cib_217": 47.2,
        "NcPlanckFI:xi_sz_cib": 0.42,
        "NcPlanckFI:A_sz": 7.23,
        "NcPlanckFI:ps_A_100_100": 250.5,
        "NcPlanckFI:ps_A_143_143": 47.4,
        "NcPlanckFI:ps_A_143_217": 47.3,
        "NcPlanckFI:ps_A_217_217": 119.8,
        "NcPlanckFI:ksz_norm": 0.01,
        "NcPlanckFI:gal545_A_100": 8.86,
        "NcPlanckFI:gal545_A_143": 11.10,
        "NcPlanckFI:gal545_A_143_217": 19.83,
        "NcPlanckFI:gal545_A_217": 95.1,
        "NcPlanckFI:calib_100T": 0.99969,
        "NcPlanckFI:calib_217T": 0.99816,
        "NcPlanckFI:A_planck": 1.00061,
        "NcPlanckFI:galf_TE_A_100": 0.1142,
        "NcPlanckFI:galf_TE_A_100_143": 0.1345,
        "NcPlanckFI:galf_TE_A_100_217": 0.482,
        "NcPlanckFI:galf_TE_A_143": 0.224,
        "NcPlanckFI:galf_TE_A_143_217": 0.664,
        "NcPlanckFI:galf_TE_A_217": 2.081,
    },
    (Planck18Types.TT, Planck18HIPrimModel.TWO_FLUIDS): {
        "NcHICosmo:H0": 66.86,
        "NcHICosmo:omegac": 0.12068,
        "NcHICosmo:omegab": 0.022126,
        "NcHIPrim:ln10e10ASA": 3.0413,
        "NcHIPrim:lnk0": -6.0,
        "NcHIPrim:lnw": np.log(1.0e-4),
        "NcHIReion:z_re": 7.54,
        "NcPlanckFI:A_cib_217": 48.5,
        "NcPlanckFI:xi_sz_cib": 0.32,
        "NcPlanckFI:A_sz": 7.03,
        "NcPlanckFI:ps_A_100_100": 254.9,
        "NcPlanckFI:ps_A_143_143": 49.8,
        "NcPlanckFI:ps_A_143_217": 47.3,
        "NcPlanckFI:ps_A_217_217": 119.9,
        "NcPlanckFI:ksz_norm": 0.00,
        "NcPlanckFI:gal545_A_100": 8.86,
        "NcPlanckFI:gal545_A_143": 10.80,
        "NcPlanckFI:gal545_A_143_217": 19.43,
        "NcPlanckFI:gal545_A_217": 94.8,
        "NcPlanckFI:calib_100T": 0.99965,
        "NcPlanckFI:calib_217T": 0.99825,
        "NcPlanckFI:A_planck": 1.00046,
    },
    (Planck18Types.TTTEEE, Planck18HIPrimModel.TWO_FLUIDS): {
        "NcHICosmo:H0": 67.32,
        "NcHICosmo:omegac": 0.12010,
        "NcHICosmo:omegab": 0.022377,
        "NcHIPrim:ln10e10ASA": 3.0447,
        "NcHIPrim:lnk0": -6.0,
        "NcHIPrim:lnw": np.log(1.0e-4),
        "NcHIReion:z_re": 7.68,
        "NcPlanckFI:A_cib_217": 47.2,
        "NcPlanckFI:xi_sz_cib": 0.42,
        "NcPlanckFI:A_sz": 7.23,
        "NcPlanckFI:ps_A_100_100": 250.5,
        "NcPlanckFI:ps_A_143_143": 47.4,
        "NcPlanckFI:ps_A_143_217": 47.3,
        "NcPlanckFI:ps_A_217_217": 119.8,
        "NcPlanckFI:ksz_norm": 0.01,
        "NcPlanckFI:gal545_A_100": 8.86,
        "NcPlanckFI:gal545_A_143": 11.10,
        "NcPlanckFI:gal545_A_143_217": 19.83,
        "NcPlanckFI:gal545_A_217": 95.1,
        "NcPlanckFI:calib_100T": 0.99969,
        "NcPlanckFI:calib_217T": 0.99816,
        "NcPlanckFI:A_planck": 1.00061,
        "NcPlanckFI:galf_TE_A_100": 0.1142,
        "NcPlanckFI:galf_TE_A_100_143": 0.1345,
        "NcPlanckFI:galf_TE_A_100_217": 0.482,
        "NcPlanckFI:galf_TE_A_143": 0.224,
        "NcPlanckFI:galf_TE_A_143_217": 0.664,
        "NcPlanckFI:galf_TE_A_217": 2.081,
    },
}


def set_mset_parameters(
    mset: Ncm.MSet, exp_type: Planck18Types, prim_model: Planck18HIPrimModel
):
    """Set the experiment parameters."""
    for param, value in EXP_PARAMETERS[(exp_type, prim_model)].items():
        pi = mset.fparam_get_pi_by_name(param)
        if pi is None:
            raise ValueError(f"Invalid parameter: {param}")

        mset.param_set(pi.mid, pi.pid, value)


def create_cosmo_for_cmb(
    massive_nu: bool = False,
    prim_model: Planck18HIPrimModel = Planck18HIPrimModel.POWER_LAW,
) -> Nc.HICosmo:
    """Create a cosmology for CMB experiments."""
    if massive_nu:
        cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    else:
        cosmo = Nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.cmb_params()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("omegab", 0.022)
    cosmo.param_set_by_name("omegac", 0.12)

    if massive_nu:
        cosmo.orig_param_set(Nc.HICosmoDESParams.ENNU, 2.0328)
        param_id = cosmo.vparam_index(Nc.HICosmoDEVParams.M, 0)
        cosmo.param_set_ftype(param_id, Ncm.ParamType.FREE)

    cosmo.set_property("H0_fit", True)
    cosmo.set_property("Omegac_fit", True)
    cosmo.set_property("Omegab_fit", True)
    cosmo.set_property("w_fit", False)
    cosmo.param_set_by_name("Omegak", 0.00)
    cosmo.set_property("Omegax_fit", False)

    if prim_model == Planck18HIPrimModel.POWER_LAW:
        prim = Nc.HIPrimPowerLaw.new()
        prim.set_property("ln10e10ASA_fit", True)
        prim.set_property("n_SA_fit", True)
    elif prim_model == Planck18HIPrimModel.TWO_FLUIDS:
        prim = Nc.HIPrimTwoFluids(use_default_calib=True)
        prim.set_property("ln10e10ASA_fit", True)
        prim.set_property("lnk0_fit", True)
        prim.set_property("lnw_fit", True)
    else:
        raise ValueError(f"Invalid primordial model: {prim_model}")

    reion = Nc.HIReionCamb.new()
    reion.set_property("z_re_fit", True)

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def create_mfunc_array_for_cmb(cbe: Nc.CBE) -> Ncm.ObjArray:
    """Create a list of extra functions for CMB experiments."""
    mfunc_oa = Ncm.ObjArray.new()

    psml = Nc.PowspecMLCBE.new_full(cbe=cbe)
    psml.set_kmin(1.0e-5)
    psml.set_kmax(1.0e1)
    psml.require_zi(0.0)
    psml.require_zf(1.0)

    psf_cbe = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf_cbe.set_best_lnr0()

    mfunc_sigma8 = Ncm.MSetFuncList.new("NcHICosmo:sigma8", psf_cbe)
    mfunc_oa.add(mfunc_sigma8)

    dist = Nc.Distance.new(10.0)

    mfunc_r_zd = Ncm.MSetFuncList.new("NcDistance:r_zd_Mpc", dist)
    mfunc_oa.add(mfunc_r_zd)

    mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
    mfunc_oa.add(mfunc_Omegam)

    mfunc_theta100 = Ncm.MSetFuncList.new("NcDistance:theta100CMB", dist)
    mfunc_oa.add(mfunc_theta100)

    return mfunc_oa


def generate_planck18_tt(
    massive_nu: bool = False,
    prim_model: Planck18HIPrimModel = Planck18HIPrimModel.POWER_LAW,
    use_lensing_likelihood: bool = False,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""
    # Likelihood
    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    if prim_model == Planck18HIPrimModel.TWO_FLUIDS:
        cbe = cbe_boltzmann.peek_cbe()
        cbe_prec = cbe.peek_precision()
        # Increase the precision for the two-fluids model
        cbe_prec.props.k_per_decade_primordial = 30.0

    b18_lowl_EE = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_EE, cbe_boltzmann
    )
    b18_lowl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe_boltzmann
    )
    b18_highl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TT, cbe_boltzmann
    )
    if not use_lensing_likelihood:
        dset = Ncm.Dataset.new_array([b18_lowl_EE, b18_lowl_TT, b18_highl_TT])
    else:
        b18_PP = Nc.DataPlanckLKL.full_new_id(
            Nc.DataPlanckLKLType.BASELINE_18_LENSING, cbe_boltzmann
        )
        dset = Ncm.Dataset.new_array([b18_lowl_EE, b18_lowl_TT, b18_highl_TT, b18_PP])
    likelihood = Ncm.Likelihood.new(dset)
    Nc.PlanckFICorTT.add_all_default18_priors(likelihood)

    # Models

    planck_model = Nc.PlanckFICorTT()
    planck_model.params_set_default_ftype()

    cosmo = create_cosmo_for_cmb(massive_nu=massive_nu, prim_model=prim_model)

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(cbe_boltzmann.peek_cbe())

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    return experiment, mfunc_oa


def generate_planck18_ttteee(
    massive_nu: bool = False,
    prim_model: Planck18HIPrimModel = Planck18HIPrimModel.POWER_LAW,
    use_lensing_likelihood: bool = False,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""
    # Likelihood
    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    if prim_model == Planck18HIPrimModel.TWO_FLUIDS:
        cbe = cbe_boltzmann.peek_cbe()
        cbe_prec = cbe.peek_precision()
        # Increase the precision for the two-fluids model
        cbe_prec.props.k_per_decade_primordial = 30.0

    b18_lowl_EE = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_EE, cbe_boltzmann
    )

    b18_lowl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe_boltzmann
    )

    b18_highl_TTTEEE = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TTTEEE, cbe_boltzmann
    )
    if not use_lensing_likelihood:
        dset = Ncm.Dataset.new_array([b18_lowl_EE, b18_lowl_TT, b18_highl_TTTEEE])
    else:
        b18_PP = Nc.DataPlanckLKL.full_new_id(
            Nc.DataPlanckLKLType.BASELINE_18_LENSING, cbe_boltzmann
        )
        dset = Ncm.Dataset.new_array(
            [b18_lowl_EE, b18_lowl_TT, b18_highl_TTTEEE, b18_PP]
        )

    likelihood = Ncm.Likelihood.new(dset)
    Nc.PlanckFICorTTTEEE.add_all_default18_priors(likelihood)

    # Models

    planck_model = Nc.PlanckFICorTTTEEE()
    planck_model.params_set_default_ftype()

    cosmo = create_cosmo_for_cmb(massive_nu=massive_nu, prim_model=prim_model)

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(cbe_boltzmann.peek_cbe())

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    return experiment, mfunc_oa
