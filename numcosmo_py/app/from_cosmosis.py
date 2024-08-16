#
# from_cosmosis.py
#
# Wed Feb 14 18:59:34 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# from_cosmosis.py
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

"""NumCosmo APP subcommand to convert CosmoSIS likelihoods to NumCosmo."""

from typing import Optional, Annotated
from pathlib import Path
import typer

from numcosmo_py import Ncm

try:
    from numcosmo_py.external.cosmosis import (
        convert_likelihoods,
        create_numcosmo_mapping,
        LinearMatterPowerSpectrum,
        NonLinearMatterPowerSpectrum,
    )
except ImportError:
    COSMOSIS = False
else:
    COSMOSIS = True


if COSMOSIS:
    # Only add the cosmosis commands if cosmosis is available

    def numcosmo_from_cosmosis(
        inifile: Annotated[Path, typer.Argument(help="Path to the Cosmosis ini file.")],
        *,
        outfile: Annotated[
            Optional[Path],
            typer.Option(
                help="Path to the output file, if not given,"
                " the input file name is used with the extension .yaml."
            ),
        ] = None,
        matter_ps: Annotated[
            LinearMatterPowerSpectrum,
            typer.Option(help="Matter power spectrum to use."),
        ] = LinearMatterPowerSpectrum.NONE,
        nonlin_matter_ps: Annotated[
            NonLinearMatterPowerSpectrum,
            typer.Option(help="Non-linear matter power spectrum to use."),
        ] = NonLinearMatterPowerSpectrum.NONE,
        distance_max_z: Annotated[
            float,
            typer.Option(
                help="Max distance to optimize distance computations", min=0.0
            ),
        ] = 10.0,
        reltol: Annotated[
            float,
            typer.Option(
                help="Relative tolerance for the distance computation", min=1.0e-14
            ),
        ] = 1.0e-4,
        mute_cosmosis: Annotated[
            bool,
            typer.Option(
                help="Mute Cosmosis output.",
            ),
        ] = False,
    ):
        """Convert a Cosmosis ini file to a NumCosmo yaml file.

        The NumCosmo yaml will contain the model builders and experiment matching the
        Cosmosis ini file. The NumCosmo yaml file can be used to run the same
        likelihoods in NumCosmo.

        Due to the differences between the two frameworks, some likelihoods may not be
        converted correctly, usually due to the different parameter names or the
        different parameterizations of the models. In this case, the user should
        manually adjust the model builders.
        """
        Ncm.cfg_init()

        if outfile is None:
            outfile = Path(inifile.stem + ".yaml")

        mapping = create_numcosmo_mapping(
            matter_ps=matter_ps,
            nonlin_matter_ps=nonlin_matter_ps,
            distance_max_z=distance_max_z,
            reltol=reltol,
        )

        model_builders, mset, likelihood = convert_likelihoods(
            inifile, mapping=mapping, mute_cosmosis=mute_cosmosis
        )

        builders_file = outfile.with_suffix(".builders.yaml")

        experiment = Ncm.ObjDictStr.new()
        experiment.add("likelihood", likelihood)
        experiment.add("model-set", mset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.dict_str_to_yaml_file(model_builders, builders_file.absolute().as_posix())
        ser.dict_str_to_yaml_file(experiment, outfile.absolute().as_posix())
