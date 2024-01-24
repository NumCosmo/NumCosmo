#
# cosmosis.py
#
# Tue Jan 23 08:22:00 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cosmosis.py
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


"""Cosmosis utilities."""

import os
import math
from typing import List, Dict, Tuple, Optional
from pathlib import Path
from enum import Enum

from cosmosis.runtime.config import Inifile
from cosmosis.runtime.parameter import Parameter
from cosmosis.runtime.prior import Prior, GaussianPrior
from cosmosis.runtime.pipeline import (
    PIPELINE_INI_SECTION,
    NO_LIKELIHOOD_NAMES,
    module,
    config_to_block,
)
from cosmosis.datablock.cosmosis_py.block import DataBlock
from firecrown.likelihood.likelihood import NamedParameters
from firecrown.connector.numcosmo.numcosmo import NumCosmoFactory, MappingNumCosmo

from numcosmo_py import Ncm, Nc, GObject, to_camel_case


class LinearMatterPowerSpectrum(str, Enum):
    """Possible linear matter power spectrum models."""

    NONE = "none"
    BBKS = "bbks"
    EISENSTEIN_HU = "eisenstein_hu"
    CLASS = "class"


class NonLinearMatterPowerSpectrum(str, Enum):
    """Possible non-linear matter power spectrum models."""

    NONE = "none"
    HALOFIT = "halofit"


def convert_parameter(p: Parameter) -> Ncm.SParam:
    """Converts a Cosmosis parameter to a NumCosmo parameter."""

    if p.limits[0] == p.limits[1]:
        if p.start == 0.0:
            p_lower = -1.0
            p_upper = 1.0
            p_scale = 0.01
        elif p.start > 0.0:
            p_lower = p.start * 0.8
            p_upper = p.start * 1.2
            p_scale = +p.start * 0.01
        else:
            p_lower = p.start * 1.2
            p_upper = p.start * 0.8
            p_scale = -p.start * 0.01
        ftype = Ncm.ParamType.FIXED
    else:
        p_lower = p.limits[0]
        p_upper = p.limits[1]
        p_scale = (p_upper - p_lower) * 0.01
        ftype = Ncm.ParamType.FREE
    return Ncm.SParam.new(
        name=p.name,
        symbol=p.name,
        lower_bound=p_lower,
        upper_bound=p_upper,
        scale=p_scale,
        abstol=0.0,
        default_val=p.start,
        ftype=ftype,
    )


def convert_single_model(
    sampling_parameters_section: str,
    parameters: List[Parameter],
) -> Tuple[str, Ncm.ModelBuilder, Ncm.Model]:
    """Converts a single sampling parameters section from a Cosmosis ini file
    to a NumCosmo model.

    Args:
        sampling_parameters_section (str): The name of the sampling parameters
            section to convert.
        parameters (List[Parameter]): The list of parameters in the ini file.

    Returns:
        Tuple[str, Ncm.ModelBuilder, Ncm.Model]: A tuple containing the name of
            the model, the model builder and the model.

    """

    model_name = to_camel_case(sampling_parameters_section)
    model_builder = Ncm.ModelBuilder.new(
        Ncm.Model,
        model_name,
        f"Model from cosmosis ini section {sampling_parameters_section}",
    )
    for p in parameters:
        if p.section == sampling_parameters_section:
            model_builder.add_sparam_obj(convert_parameter(p))

    FirecrownModel = model_builder.create()  # pylint: disable=invalid-name
    GObject.new(FirecrownModel)
    NcmFirecrownModel = FirecrownModel.pytype  # pylint: disable=invalid-name
    GObject.type_register(NcmFirecrownModel)
    return model_name, model_builder, NcmFirecrownModel()


def convert_models(
    sampling_parameters_sections: List[str],
    parameters: List[Parameter],
    model_builders: Ncm.ObjDictStr,
    mset: Ncm.MSet,
) -> List[str]:
    """Converts a list of sampling parameters sections from a Cosmosis ini file
    to NumCosmo models."""

    model_names_list = []
    for sampling_parameters_section in sampling_parameters_sections:
        model_name, model_builder, model = convert_single_model(
            sampling_parameters_section,
            parameters,
        )
        model_names_list.append(model_name)
        model_builders.add(sampling_parameters_section, model_builder)
        mset.set(model)

    return model_names_list


def convert_single_likelihood(
    config: DataBlock,
    likelihood_name: str,
    module_list: List[str],
    parameters: List[Parameter],
    priors: Dict[Tuple[str, str], Prior],
    mapping: MappingNumCosmo,
    model_builders: Ncm.ObjDictStr,
    mset: Ncm.MSet,
    likelihood: Ncm.Likelihood,
) -> None:
    """Converts a single likelihood from a Cosmosis ini file to a NumCosmo
    likelihood.

    Args:
        ini (Inifile): The Cosmosis ini file.
        likelihood_name (str): The name of the likelihood to convert.
        module_list (List[str]): The list of modules in the ini file.
        parameters (List[Parameter]): The list of parameters in the ini file.
        priors (Dict[Tuple[str, str], Prior]): The priors in the ini file.
        model_builders (Ncm.ObjDictStr): The model builders to add the models to.
        mset (Ncm.MSet): The model set to add the models to.
        likelihood (Ncm.Likelihood): The likelihood to add the data to.
    """
    # CosmoSIS adds '_likelihood' to the likelihood name to get the module name.
    likelihood_module = likelihood_name + "_likelihood"

    # Check if the likelihood is a firecrown likelihood.
    # Only firecrown likelihoods are supported at the moment and their
    # name must contain 'firecrown'.
    if "firecrown" not in likelihood_name:
        raise ValueError(f"Likelihood {likelihood_name} is not a firecrown likelihood.")

    # Check if the likelihood is in the module list.
    if likelihood_module not in module_list:
        raise ValueError(
            f"Likelihood {likelihood_name} is not in the module list {module_list}."
        )

    # Gets the sampling parameters sections.
    sampling_parameters_sections = config.get(
        likelihood_module, "sampling_parameters_sections"
    ).split()

    # Converts the sampling parameters sections to NumCosmo models.
    model_names_list = convert_models(
        sampling_parameters_sections,
        parameters,
        model_builders,
        mset,
    )

    # Gets the likelihood source.
    try:
        likelihood_source = config[likelihood_module, "likelihood_source"]
    except KeyError as e:
        raise ValueError(
            f"Likelihood {likelihood_name} does not have a likelihood_source."
        ) from e

    # Gets the likelihood parameters.
    build_parameters_dict = {
        name: config[likelihood_module, name]
        for _, name in config.keys(section=likelihood_module)
    }
    build_parameters = NamedParameters(build_parameters_dict)

    # Using Firecrown's NumCosmoFactory to create the NumCosmo likelihood.
    numcosmo_factory = NumCosmoFactory(
        likelihood_source=likelihood_source,
        build_parameters=build_parameters,
        mapping=mapping,
        model_list=model_names_list,
    )

    # Getting the NumCosmo likelihood and removing the data from the dataset.
    # This results in a clean serialization of the likelihood.
    firecrown_data = numcosmo_factory.get_data()
    if isinstance(firecrown_data, Ncm.DataGaussCov):
        firecrown_data.set_size(0)
    likelihood.peek_dataset().append_data(firecrown_data)

    for (section, name), prior in priors.items():
        # Adding only Gaussian priors, uniform and delta priors are included
        # in the limits of the parameters and by choosing FIXED or FREE.
        # If more non-trivial priors are added to Cosmosis, they should be
        # added here.
        if isinstance(prior, GaussianPrior):
            gp = Ncm.PriorGaussParam.new_name(
                f"{to_camel_case(section)}:{name}", prior.mu, prior.sigma
            )
            likelihood.priors_add(gp)


def create_numcosmo_mapping(
    matter_ps: LinearMatterPowerSpectrum = LinearMatterPowerSpectrum.NONE,
    nonlin_matter_ps: NonLinearMatterPowerSpectrum = NonLinearMatterPowerSpectrum.NONE,
    distance_max_z: float = 10.0,
) -> MappingNumCosmo:
    """Creates a NumCosmo mapping to be used in the likelihoods.
    converted from Cosmosis."""

    ps_ml = None
    ps_mnl = None
    dist = Nc.Distance.new(distance_max_z)

    if matter_ps == LinearMatterPowerSpectrum.BBKS:
        transfer_bbks = Nc.TransferFuncBBKS.new()
        ps_ml = Nc.PowspecMLTransfer.new(transfer_bbks)
    elif matter_ps == LinearMatterPowerSpectrum.EISENSTEIN_HU:
        transfer_eh = Nc.TransferFuncEH.new()
        ps_ml = Nc.PowspecMLTransfer.new(transfer_eh)
    elif matter_ps == LinearMatterPowerSpectrum.CLASS:
        ps_ml = Nc.PowspecMLCBE.new()

    if nonlin_matter_ps == NonLinearMatterPowerSpectrum.HALOFIT:
        if ps_ml is None:
            raise ValueError(
                "Non-linear matter power spectrum is HALOFIT but linear matter"
                " power spectrum is not set."
            )
        ps_mnl = Nc.PowspecMNLHaloFit.new(
            psml=ps_ml, zmaxnl=distance_max_z, reltol=1.0e-5
        )

    return MappingNumCosmo(p_ml=ps_ml, p_mnl=ps_mnl, dist=dist)


def convert_cosmology(
    params_datablock: DataBlock,
) -> Nc.HICosmo:
    """Converts the cosmology objects described in a Cosmosis ini file to a NumCosmo
    cosmology."""

    COSMOLOGY_SECTION = "cosmological_parameters"
    params_dict: Dict[str, float] = {}

    for section, name in params_datablock.keys():
        if section == COSMOLOGY_SECTION:
            params_dict[name] = params_datablock[section, name]

    massnu_length = 0
    if "mnu" in params_dict:
        massnu_length = 1

    cosmo = Nc.HICosmoDECpl(massnu_length=massnu_length)
    cosmo.omega_x2omega_k()

    if "h0" in params_dict:
        cosmo.param_set_by_name("H0", params_dict["h0"] * 100.0)
    if "omega_k" in params_dict:
        cosmo.param_set_by_name("Omegak", params_dict["omega_k"])
    if "omega_b" in params_dict:
        cosmo.param_set_by_name("Omegab", params_dict["omega_b"])
    if "omega_c" in params_dict:
        cosmo.param_set_by_name("Omegac", params_dict["omega_c"])

    if "mnu" in params_dict:
        cosmo.param_set_by_name("massnu_0", params_dict["mnu"])
    if "nnu" in params_dict:
        cosmo.param_set_by_name("ENnu", params_dict["nnu"])
    if "tcmb" in params_dict:
        cosmo.param_set_by_name("Tgamma0", params_dict["tcmb"])

    if "w" in params_dict:
        cosmo.param_set_by_name("w0", params_dict["w"])
    if "wa" in params_dict:
        cosmo.param_set_by_name("w1", params_dict["wa"])

    prim = Nc.HIPrimPowerLaw.new()
    if "a_s" in params_dict:
        prim.param_set_by_name("ln10e10ASA", math.log(1.0e10 * params_dict["a_s"]))

    if "n_s" in params_dict:
        prim.param_set_by_name("n_SA", params_dict["n_s"])

    reion = Nc.HIReionCamb.new()
    reion.set_z_from_tau(cosmo, 0.0561)

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def convert_likelihoods(
    inifile: Path,
    mapping: Optional[MappingNumCosmo] = None,
) -> Tuple[Ncm.ObjDictStr, Ncm.MSet, Ncm.Likelihood]:
    """Converts the likelihoods from a Cosmosis ini file to NumCosmo
    likelihoods.

    Args:
        inifile (Path): Path to the Cosmosis ini file.
    """

    # Loads the ini file.
    ini = Inifile(inifile)

    # Converts the ini file to a datablock.
    relevant_sections = ini.sections()
    global_sections = ini.get("runtime", "global", fallback=" ")
    for global_section in global_sections.split():
        relevant_sections.append(global_section)

    config_block = config_to_block(relevant_sections, ini)

    # Gets the list of modules.
    module_list = ini.get(PIPELINE_INI_SECTION, "modules", fallback="").split()

    # We need to load the consistency module to convert the cosmological
    # parameters to NumCosmo.
    consistency = module.Module.from_options(
        "consistency",
        ini,
        root_directory=ini.get("runtime", "root", fallback=os.getcwd()),
    )

    # Gets the likelihood names.
    likelihood_names: List[str] = ini.get(
        PIPELINE_INI_SECTION, "likelihoods", fallback=NO_LIKELIHOOD_NAMES
    ).split()

    # Gets the values and priors files.
    values_file = ini.get(PIPELINE_INI_SECTION, "values", fallback="")
    priors_files = ini.get(PIPELINE_INI_SECTION, "priors", fallback="").split()

    # Loads the parameters and priors.
    parameters: List[Parameter] = Parameter.load_parameters(values_file, priors_files)
    priors: Dict[Tuple[str, str], Prior] = Prior.load_priors(priors_files)

    # Putting the cosmological parameters in a datablock.
    params_datablock = DataBlock()
    for param in parameters:
        params_datablock[param.section, param.name] = param.start

    # Running the consistency module.
    config_block[PIPELINE_INI_SECTION, "current_module"] = consistency.name
    consistency.setup(config_block)
    consistency.execute(params_datablock)

    model_builders = Ncm.ObjDictStr.new()
    mset = Ncm.MSet.empty_new()
    dataset = Ncm.Dataset.new()
    likelihood = Ncm.Likelihood.new(dataset)
    if mapping is None:
        mapping = create_numcosmo_mapping()

    # Converts the cosmology.
    cosmo = convert_cosmology(params_datablock)
    mset.set(cosmo)

    for likelihood_name in likelihood_names:
        convert_single_likelihood(
            config_block,
            likelihood_name,
            module_list,
            parameters,
            priors,
            mapping,
            model_builders,
            mset,
            likelihood,
        )

    return model_builders, mset, likelihood
