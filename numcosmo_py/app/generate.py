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

import numpy as np
import typer

from numcosmo_py import Ncm
from numcosmo_py.experiments.planck18 import (
    Planck18Types,
    Planck18HIPrimModel,
    generate_planck18_tt,
    generate_planck18_ttteee,
    set_mset_parameters,
)
from numcosmo_py.experiments.jpas_forecast24 import (
    ClusterRedshiftType,
    ClusterMassType,
    JpasSSCType,
    generate_jpas_forecast_2024,
)


@dataclasses.dataclass(kw_only=True)
class GeneratePlanck:
    """Common block for commands that load an experiment.

    All commands that load an experiment should inherit from this class.
    """

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    data_type: Annotated[
        Planck18Types, typer.Option(help="Data type to use.", show_default=True)
    ] = Planck18Types.TT

    prim_model: Annotated[
        Planck18HIPrimModel,
        typer.Option(help="Primordial model to use.", show_default=True),
    ] = Planck18HIPrimModel.POWER_LAW

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
            exp, mfunc_array = generate_planck18_tt(
                massive_nu=self.massive_nu, prim_model=self.prim_model
            )
        elif self.data_type == Planck18Types.TTTEEE:
            exp, mfunc_array = generate_planck18_ttteee(
                massive_nu=self.massive_nu, prim_model=self.prim_model
            )
        else:
            raise ValueError(f"Invalid data type: {self.data_type}")

        set_mset_parameters(exp.peek("model-set"), self.data_type, self.prim_model)

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


@dataclasses.dataclass(kw_only=True)
class GenerateJpasForecast:
    """Common block for commands that load an experiment.

    All commands that load an experiment should inherit from this class.
    """

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

    resample_model: Annotated[
        tuple[float, float, float],
        typer.Option(
            help=(
                "Model for "
                "resample.(NcHICosmo:Omegac,NcHICosmo:w,NcHIPrim:ln10e10ASA)"
            ),
            show_default=True,
        ),
    ] = (0.2612, -1.0, 3.027)

    fitting_model: Annotated[
        tuple[float, float, float],
        typer.Option(
            help=(
                "Model for fitting. "
                "(NcHICosmo:Omegac,NcHICosmo:w,NcHIPrim:ln10e10ASA)"
            ),
            show_default=True,
        ),
    ] = (0.2612, -1.0, 3.027)

    use_fixed_cov: Annotated[
        bool, typer.Option(help="Use fixed covariance matrix.", show_default=True)
    ] = False

    z_min: Annotated[
        float,
        typer.Option(help="Jpas minimum redshift.", show_default=True, min=0),
    ] = 0.1

    z_max: Annotated[
        float,
        typer.Option(help="Jpas maximum redshift.", show_default=True, max=2.0),
    ] = 0.8

    znknots: Annotated[
        int,
        typer.Option(help="Jpas number of redshift bins.", show_default=True, min=2),
    ] = 8

    cluster_redshift_type: Annotated[
        Optional[ClusterRedshiftType],
        typer.Option(help="Cluster photoz relation.", show_default=True),
    ] = ClusterRedshiftType.NODIST

    lnM_min: Annotated[
        float,
        typer.Option(
            help="Jpas minimum mass.",
            show_default=True,
        ),
    ] = (
        np.log(10.0) * 14.0
    )

    lnM_max: Annotated[
        float,
        typer.Option(
            help="Jpas maximum mass.",
            show_default=True,
        ),
    ] = (
        np.log(10.0) * 15.0
    )

    lnMnknots: Annotated[
        int,
        typer.Option(help="Jpas number of mass bins.", show_default=True, min=2),
    ] = 2

    cluster_mass_type: Annotated[
        Optional[ClusterMassType],
        typer.Option(help="Cluster mass-observable relation.", show_default=True),
    ] = ClusterMassType.NODIST

    survey_area: Annotated[
        float,
        typer.Option(
            help=(
                "Jpas survey area. This option is unvailable "
                "for the partial sky cases."
            ),
            show_default=True,
            min=0,
        ),
    ] = 2959.1

    def __post_init__(self):
        """Generate JPAS 2024 forecast experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        exp, mfunc_array = generate_jpas_forecast_2024(
            area=self.survey_area,
            z_min=self.z_min,
            z_max=self.z_max,
            znknots=self.znknots,
            cluster_redshift_type=self.cluster_redshift_type,
            lnM_min=self.lnM_min,
            lnM_max=self.lnM_max,
            lnMnknots=self.lnMnknots,
            cluster_mass_type=self.cluster_mass_type,
            use_fixed_cov=self.use_fixed_cov,
            fitting_Sij_type=self.fitting_sky_cut,
            resample_Sij_type=self.resample_sky_cut,
            resample_model=self.resample_model,
            fitting_model=self.fitting_model,
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
