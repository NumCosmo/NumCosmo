#
# generate.py
#
# Mon Feb 19 16:33:18 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# generate.py
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

"""NumCosmo APP subcommands generate experiment files."""

from typing import Annotated
from enum import Enum
import dataclasses
from pathlib import Path
import typer

from numcosmo_py import Ncm, Nc


class Planck18Types(str, Enum):
    """Planck 18 baseline data combinations."""

    TT = "TT"
    TTTEEE = "TTTEEE"


def create_cosmo_for_cmb() -> Nc.HICosmo:
    """Create a cosmology for CMB experiments."""
    cosmo = Nc.HICosmoDEXcdm(massnu_length=1)
    cosmo.params_set_default_ftype()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegac", 0.25)

    cosmo.orig_param_set(Nc.HICosmoDESParams.ENNU, 2.0328)
    param_id = cosmo.vparam_index(Nc.HICosmoDEVParams.M, 0)
    cosmo.param_set_ftype(param_id, Ncm.ParamType.FIXED)

    cosmo.set_property("H0_fit", True)
    cosmo.set_property("Omegac_fit", True)
    cosmo.set_property("Omegab_fit", True)
    cosmo.set_property("w_fit", True)
    cosmo.param_set_by_name("Omegak", 0.00)
    cosmo.set_property("Omegax_fit", False)

    prim = Nc.HIPrimPowerLaw.new()
    prim.set_property("ln10e10ASA_fit", True)
    prim.set_property("n_SA_fit", True)

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


def generate_planck18_tt() -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""

    # Likelihood

    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    b18_lowl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe_boltzmann
    )
    b18_highl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TT, cbe_boltzmann
    )

    dset = Ncm.Dataset.new_array([b18_lowl_TT, b18_highl_TT])
    likelihood = Ncm.Likelihood.new(dset)
    Nc.PlanckFICorTT.add_all_default18_priors(likelihood)

    # Models

    planck_model = Nc.PlanckFICorTT()
    planck_model.params_set_default_ftype()

    cosmo = create_cosmo_for_cmb()

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(cbe_boltzmann.peek_cbe())

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    return experiment, mfunc_oa


def generate_planck18_ttteee() -> tuple[Ncm.ObjDictStr, Ncm.ObjArray]:
    """Generate Planck 2018 TT baseline experiment dictionary."""

    # Likelihood

    cbe_boltzmann = Nc.HIPertBoltzmannCBE.new()
    b18_lowl_TT = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe_boltzmann
    )

    b18_lowl_EE = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_LOWL_EE, cbe_boltzmann
    )

    b18_highl_TTTEEE = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TTTEEE, cbe_boltzmann
    )

    dset = Ncm.Dataset.new_array([b18_lowl_TT, b18_lowl_EE, b18_highl_TTTEEE])
    likelihood = Ncm.Likelihood.new(dset)
    Nc.PlanckFICorTTTEEE.add_all_default18_priors(likelihood)

    # Models

    planck_model = Nc.PlanckFICorTTTEEE()
    planck_model.params_set_default_ftype()

    cosmo = create_cosmo_for_cmb()

    mset = Ncm.MSet.new_array([planck_model, cosmo])
    mset.prepare_fparam_map()

    # Extra functions

    mfunc_oa = create_mfunc_array_for_cmb(cbe_boltzmann.peek_cbe())

    # Save experiment

    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    return experiment, mfunc_oa


@dataclasses.dataclass
class GeneratePlanck:
    """Common block for commands that load an experiment. All commands that load an
    experiment should inherit from this class."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    data_type: Annotated[
        Planck18Types, typer.Option(help="Data type to use.", show_default=True)
    ] = Planck18Types.TT

    def __post_init__(self):
        """Generate Planck 2018 TT baseline experiment."""

        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        if self.data_type == Planck18Types.TT:
            exp, mfunc_array = generate_planck18_tt()
        elif self.data_type == Planck18Types.TTTEEE:
            exp, mfunc_array = generate_planck18_ttteee()
        else:
            raise ValueError(f"Invalid data type: {self.data_type}")

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())
        ser.array_to_yaml_file(
            mfunc_array,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )
