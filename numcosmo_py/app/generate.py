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

from typing import Annotated, Optional
import dataclasses
from pathlib import Path
import typer

from numcosmo_py import Ncm
from numcosmo_py.experiments.planck18 import (
    Planck18Types,
    generate_planck18_tt,
    generate_planck18_ttteee,
    set_mset_parameters,
)
from numcosmo_py.experiments.jpas_forecast24 import (
    JpasSSCType,
    JpasModelTypes,
    generate_jpas_forecast_2024,
)


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

    massive_nu: Annotated[
        bool, typer.Option(help="Use massive neutrinos.", show_default=True)
    ] = False

    def __post_init__(self):
        """Generate Planck 2018 TT baseline experiment."""

        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        if self.data_type == Planck18Types.TT:
            exp, mfunc_array = generate_planck18_tt(massive_nu=self.massive_nu)
        elif self.data_type == Planck18Types.TTTEEE:
            exp, mfunc_array = generate_planck18_ttteee(massive_nu=self.massive_nu)
        else:
            raise ValueError(f"Invalid data type: {self.data_type}")

        set_mset_parameters(exp.peek("model-set"), self.data_type)

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


@dataclasses.dataclass
class GenerateJpasForecast:
    """Common block for commands that load an experiment. All commands that load an
    experiment should inherit from this class."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    fitting_sky_cut: Annotated[
        Optional[JpasSSCType],
        typer.Option(
            help="Super Sample Covariance method for fitting.", show_default=True
        ),
    ] = JpasSSCType.FULLSKY

    resample_sky_cut: Annotated[
        JpasSSCType,
        typer.Option(
            help="Super Sample Covariance method for resampling.",
            show_default=True,
        ),
    ] = JpasSSCType.NO_SSC

    model: Annotated[
        JpasModelTypes,
        typer.Option(
            help="J-Pas model parameters to use during resample.", show_default=True
        ),
    ] = JpasModelTypes.FIDUCIAL

    model_for_cov: Annotated[
        Optional[JpasModelTypes],
        typer.Option(
            help=(
                "J-Pas model parameters to use during covariance matrix generation. "
                "This option is only used if use_fixed_cov is True."
            ),
            show_default=True,
        ),
    ] = None

    use_fixed_cov: Annotated[
        bool, typer.Option(help="Use fixed covariance matrix.", show_default=True)
    ] = False

    def __post_init__(self):
        """Generate JPAS 2024 forecast experiment."""

        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        exp, mfunc_array = generate_jpas_forecast_2024(
            use_fixed_cov=self.use_fixed_cov,
            fitting_Sij_type=self.fitting_sky_cut,
            resample_Sij_type=self.resample_sky_cut,
            model_for_resample=self.model,
            model_for_cov=self.model_for_cov,
        )

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
