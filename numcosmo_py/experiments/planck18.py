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

"""Factory functions to generate Planck18 likelihood and models.

NumCosmo provides independent reimplementations of the Planck likelihoods,
including support for alternative resampled representations. The underlying
likelihood methodology and original data products are due to the Planck
Collaboration. Analyses using these likelihoods should cite the corresponding
Planck publications in addition to NumCosmo; see ``planck_native_provenance.md``
for the source data and the reference list.
"""

from enum import StrEnum
from typing import Any, cast

import numpy as np

from numcosmo_py import Nc, Ncm
from numcosmo_py.cosmology import HIPrimModel, ParameterDesc, create_cosmo


class Planck18Types(StrEnum):
    """Planck 18 baseline data combinations."""

    TT = "TT"
    TTTEEE = "TTTEEE"


# The distributed Ncm.MSetPreset of each data combination: the Planck 2018 base LCDM
# best fits of plikHM_TT_lowl_lowE and plikHM_TTTEEE_lowl_lowE. Each has a cosmology
# preset and a "_nuisance" preset for the NcPlanckFI parameters.
PRESET_NAMES: dict[Planck18Types, str] = {
    Planck18Types.TT: "planck2018_tt_lowl_lowe",
    Planck18Types.TTTEEE: "planck2018_ttteee_lowl_lowe",
}

# Starting values of the two-fluids primordial parameters. They are NumCosmo choices,
# not published values.
TWO_FLUIDS_PARAMETERS: dict[str, ParameterDesc] = {
    "NcHIPrim:lnk0": ParameterDesc(value=-6.0),
    "NcHIPrim:lnw": ParameterDesc(value=float(np.log(1.0e-4))),
}


def mset_set_parameters(
    mset: Ncm.MSet, exp_type: Planck18Types, prim_model: HIPrimModel
) -> None:
    """Apply the Planck 2018 presets of `exp_type` to `mset`.

    `mset` holds the cosmology of create_cosmo() and the NcPlanckFI model of the data
    combination. The presets set the best fits, the marginalized standard deviations as
    scales and the fit flags (see Ncm.MSetPreset). The two-fluids primordial model has
    no n_SA: the cosmology preset is applied without it, and lnk0 and lnw take the
    values of TWO_FLUIDS_PARAMETERS. A failed preset raises GLib.Error.
    """
    name = PRESET_NAMES[exp_type]
    cosmo_preset = Ncm.MSetPreset.new_from_name(name)
    nuisance_preset = Ncm.MSetPreset.new_from_name(f"{name}_nuisance")
    assert cosmo_preset is not None and nuisance_preset is not None

    if prim_model == HIPrimModel.TWO_FLUIDS:
        assert cosmo_preset.remove("NcHIPrim:n_SA")

    cosmo_preset.apply(mset)
    nuisance_preset.apply(mset)

    if prim_model == HIPrimModel.TWO_FLUIDS:
        for param, desc in TWO_FLUIDS_PARAMETERS.items():
            model, pname = param.split(":")
            mset[model].param_set_desc(pname, cast(dict[str, Any], desc))


def create_mfunc_array_for_cmb(
    dist: Nc.Distance, psf: Ncm.PowspecFilter
) -> Ncm.ObjArray:
    """Create a list of extra functions for CMB experiments."""
    mfunc_oa = Ncm.ObjArray.new()

    mfunc_sigma8 = Ncm.MSetFuncList.new("NcHICosmo:sigma8", psf)
    mfunc_oa.add(mfunc_sigma8)

    mfunc_r_zd = Ncm.MSetFuncList.new("NcDistance:r_zd_Mpc", dist)
    mfunc_oa.add(mfunc_r_zd)

    mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
    mfunc_oa.add(mfunc_Omegam)

    mfunc_theta100 = Ncm.MSetFuncList.new("NcDistance:theta100CMB", dist)
    mfunc_oa.add(mfunc_theta100)

    return mfunc_oa


def generate_planck18_tt(
    massive_nu: bool = True,
    fit_nu_mass: bool = False,
    prim_model: HIPrimModel = HIPrimModel.POWER_LAW,
    use_lensing_likelihood: bool = False,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""
    # Likelihood
    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    if prim_model == HIPrimModel.TWO_FLUIDS:
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

    cosmo = create_cosmo(
        massive_nu=massive_nu, fit_nu_mass=fit_nu_mass, prim_model=prim_model
    )

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Distance
    dist = Nc.Distance.new(10.0)

    # Linear matter power spectrum
    cbe = cbe_boltzmann.peek_cbe()
    psml = Nc.PowspecMLCBE.new_full(cbe=cbe)
    psml.set_kmin(1.0e-5)
    psml.set_kmax(1.0e1)
    psml.require_zi(0.0)
    psml.require_zf(1.0)

    # Linear matter power spectrum filter
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(dist, psf)

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)
    experiment.set("distance", dist)
    experiment.set("ps-ml", psml)
    experiment.set("ps-ml-filter", psf)

    return experiment, mfunc_oa


def generate_planck18_ttteee(
    massive_nu: bool = True,
    fit_nu_mass: bool = False,
    prim_model: HIPrimModel = HIPrimModel.POWER_LAW,
    use_lensing_likelihood: bool = False,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""
    # Likelihood
    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    if prim_model == HIPrimModel.TWO_FLUIDS:
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

    cosmo = create_cosmo(
        massive_nu=massive_nu, fit_nu_mass=fit_nu_mass, prim_model=prim_model
    )

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Distance
    dist = Nc.Distance.new(10.0)

    # Linear matter power spectrum
    cbe = cbe_boltzmann.peek_cbe()
    psml = Nc.PowspecMLCBE.new_full(cbe=cbe)
    psml.set_kmin(1.0e-5)
    psml.set_kmax(1.0e1)
    psml.require_zi(0.0)
    psml.require_zf(1.0)

    # Linear matter power spectrum filter
    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(dist, psf)

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)
    experiment.set("distance", dist)
    experiment.set("ps-ml", psml)
    experiment.set("ps-ml-filter", psf)

    return experiment, mfunc_oa


def generate_planck18_native(
    data_type: Planck18Types,
    massive_nu: bool = True,
    fit_nu_mass: bool = False,
    prim_model: HIPrimModel = HIPrimModel.POWER_LAW,
    use_lensing_likelihood: bool = False,
    from_release: bool = False,
) -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate a Planck 2018 experiment from the native (clik-free) likelihoods.

    Assembles the low-l EE (simall), low-l TT (commander) and high-l TT/TTTEEE
    (SMICA) blocks -- plus optional lensing -- sharing one CBE. Unlike the legacy
    ``generate_planck18_*`` (which wraps the clik library through ``DataPlanckLKL``),
    the native objects embed their data and can resample; each self-configures the
    shared Boltzmann in prepare(), so the resulting file reloads and fits without
    the clik data files or the PLC library.

    With ``from_release=True`` the blocks are downloaded (and cached) from the
    NumCosmo native-Planck GitHub release instead of being built from a local
    ``plc_3.0`` clik tree -- no Planck data or PLC library needed at all.

    See ``numcosmo_py/experiments/planck_native_provenance.md`` for the data
    provenance and the required Planck Collaboration citations.
    """
    from numcosmo_py.experiments.planck_commander import (
        COMMANDER_RELPATH,
        build_commander,
    )
    from numcosmo_py.experiments.planck_lensing import (
        LENSING_FULL_RELPATH,
        build_lensing,
    )
    from numcosmo_py.experiments.planck_lite import find_baseline_file
    from numcosmo_py.experiments.planck_native_release import (
        PlanckReleaseId,
        load_planck_release,
    )
    from numcosmo_py.experiments.planck_simall import SIMALL_EE_RELPATH, build_simall
    from numcosmo_py.experiments.planck_smica import (
        PLIK_TT_RELPATH,
        PLIK_TTTEEE_RELPATH,
        build_smica_tt,
        build_smica_ttteee,
    )

    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    if prim_model == HIPrimModel.TWO_FLUIDS:
        cbe = cbe_boltzmann.peek_cbe()
        cbe.peek_precision().props.k_per_decade_primordial = 30.0

    def _block(release_id: PlanckReleaseId, relpath: str, builder) -> Ncm.Data:
        """Load a block from the release, or build it from the local clik data."""
        if from_release:
            return load_planck_release(release_id, cbe_boltzmann)
        path = find_baseline_file(relpath)
        if path is None:
            raise FileNotFoundError(
                f"Planck baseline data not found ({relpath}); build from the local "
                "clik data needs it, or pass from_release=True."
            )
        return builder(path, cbe_boltzmann)

    # Each native likelihood raises the shared Boltzmann's targets/lmax itself in
    # prepare() (nc_hipert_boltzmann_require), including after the experiment is
    # reloaded, so no manual CBE configuration is needed here.
    lowl_EE = _block(PlanckReleaseId.PR3_SIMALL_EE, SIMALL_EE_RELPATH, build_simall)
    lowl_TT = _block(PlanckReleaseId.PR3_COMMANDER, COMMANDER_RELPATH, build_commander)

    if data_type == Planck18Types.TT:
        highl = _block(PlanckReleaseId.PR3_PLIK_TT, PLIK_TT_RELPATH, build_smica_tt)
        planck_model: Nc.PlanckFI = Nc.PlanckFICorTT()
    elif data_type == Planck18Types.TTTEEE:
        highl = _block(
            PlanckReleaseId.PR3_PLIK_TTTEEE, PLIK_TTTEEE_RELPATH, build_smica_ttteee
        )
        planck_model = Nc.PlanckFICorTTTEEE()
    else:
        raise ValueError(f"Invalid data type: {data_type}")

    blocks = [lowl_EE, lowl_TT, highl]
    if use_lensing_likelihood:
        blocks.append(
            _block(PlanckReleaseId.PR3_LENSING, LENSING_FULL_RELPATH, build_lensing)
        )

    dset = Ncm.Dataset.new_array(blocks)
    likelihood = Ncm.Likelihood.new(dset)

    if data_type == Planck18Types.TT:
        Nc.PlanckFICorTT.add_all_default18_priors(likelihood)
    else:
        Nc.PlanckFICorTTTEEE.add_all_default18_priors(likelihood)

    planck_model.params_set_default_ftype()

    cosmo = create_cosmo(
        massive_nu=massive_nu, fit_nu_mass=fit_nu_mass, prim_model=prim_model
    )
    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    dist = Nc.Distance.new(10.0)

    cbe = cbe_boltzmann.peek_cbe()
    psml = Nc.PowspecMLCBE.new_full(cbe=cbe)
    psml.set_kmin(1.0e-5)
    psml.set_kmax(1.0e1)
    psml.require_zi(0.0)
    psml.require_zf(1.0)

    psf = Ncm.PowspecFilter.new(psml, Ncm.PowspecFilterType.TOPHAT)
    psf.set_best_lnr0()

    mfunc_oa = create_mfunc_array_for_cmb(dist, psf)

    experiment = Ncm.ObjDictStr()
    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)
    experiment.set("distance", dist)
    experiment.set("ps-ml", psml)
    experiment.set("ps-ml-filter", psf)

    return experiment, mfunc_oa
