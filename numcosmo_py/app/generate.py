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

from numcosmo_py import Ncm, Nc
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
from numcosmo_py.experiments.cluster_wl import (
    generate_lsst_cluster_wl,
    GalaxySDShapeDist,
    GalaxyZDist,
)
from numcosmo_py.experiments.xcdm_no_perturbations import SNIaID, add_snia_likelihood


@dataclasses.dataclass(kw_only=True)
class GeneratePlanck:
    """Generate Planck 2018 experiment."""

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

    include_lens_lkl: Annotated[
        bool, typer.Option(help="Include lensing likelihood.", show_default=True)
    ] = False

    include_snia: Annotated[
        Optional[SNIaID], typer.Option(help="Include SNIa data.", show_default=True)
    ] = None

    include_des_y3_S8_prior: Annotated[
        bool, typer.Option(help="Include DES Year 3 S8 prior.", show_default=True)
    ] = False

    def __post_init__(self) -> None:
        """Generate Planck 2018 TT baseline experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        if self.data_type == Planck18Types.TT:
            exp, mfunc_array = generate_planck18_tt(
                massive_nu=self.massive_nu,
                prim_model=self.prim_model,
                use_lensing_likelihood=self.include_lens_lkl,
            )
        elif self.data_type == Planck18Types.TTTEEE:
            exp, mfunc_array = generate_planck18_ttteee(
                massive_nu=self.massive_nu,
                prim_model=self.prim_model,
                use_lensing_likelihood=self.include_lens_lkl,
            )
        else:
            raise ValueError(f"Invalid data type: {self.data_type}")

        mset = exp.peek("model-set")
        assert isinstance(mset, Ncm.MSet)

        likelihood = exp.peek("likelihood")
        assert isinstance(likelihood, Ncm.Likelihood)

        dataset = likelihood.peek_dataset()
        assert isinstance(dataset, Ncm.Dataset)

        set_mset_parameters(mset, self.data_type, self.prim_model)

        if self.include_snia is not None:
            dist = exp.get("distance")
            assert isinstance(dist, Nc.Distance)
            add_snia_likelihood(dataset, mset, dist, self.include_snia)
            cosmo = mset.peek(Nc.HICosmo.id())
            cosmo.set_property("w_fit", True)

        if self.include_des_y3_S8_prior:
            psf = exp.get("ps-ml-filter")
            assert isinstance(psf, Ncm.PowspecFilter)
            func = Ncm.MSetFuncList.new("NcHICosmo:S8", psf)
            prior = Ncm.PriorGaussFunc.new(func, 0.775, 0.025, 0.0)
            likelihood.priors_add(prior)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        ser.to_binfile(
            dataset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())
        ser.array_to_yaml_file(
            mfunc_array,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )


@dataclasses.dataclass(kw_only=True)
class GenerateJpasForecast:
    """Generate JPAS 2024 forecast experiment."""

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
                "Model for resample. (NcHICosmo:Omegac,NcHICosmo:w,NcHICosmo:sigma8)"
            ),
            show_default=True,
        ),
    ] = (0.2612, -1.0, 0.8159)

    fitting_model: Annotated[
        tuple[float, float, float],
        typer.Option(
            help="Model for fitting. (NcHICosmo:Omegac,NcHICosmo:w,NcHICosmo:sigma8)",
            show_default=True,
        ),
    ] = (0.2612, -1.0, 0.8159)

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
    
    resample_seed: Annotated[
        int,
        typer.Option(
            help=(
                "Seed used to generate experiment."
            ),
            show_default=True,
            min=0,
        ),
    ] = 1234

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
            resample_seed=self.resample_seed,
            fitting_model=self.fitting_model,
        )

        mset = exp.peek("model-set")
        assert isinstance(mset, Ncm.MSet)

        likelihood = exp.peek("likelihood")
        assert isinstance(likelihood, Ncm.Likelihood)

        dataset = likelihood.peek_dataset()
        assert isinstance(dataset, Ncm.Dataset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        ser.to_binfile(
            dataset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())
        ser.array_to_yaml_file(
            mfunc_array,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )


@dataclasses.dataclass(kw_only=True)
class GenerateClusterWL:
    """Generate LSST cluster weak lensing experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    cluster_ra: Annotated[
        float, typer.Option(help="Cluster right ascension.", show_default=True)
    ] = 12.34

    cluster_dec: Annotated[
        float, typer.Option(help="Cluster declination.", show_default=True)
    ] = -55.123

    cluster_z: Annotated[
        float, typer.Option(help="Cluster redshift.", show_default=True)
    ] = 0.2

    cluster_mass: Annotated[
        float, typer.Option(help="Cluster mass.", show_default=True)
    ] = 1.0e14

    cluster_mass_min: Annotated[
        float, typer.Option(help="Minimum cluster mass.", show_default=True)
    ] = 1.0e13

    cluster_mass_max: Annotated[
        float, typer.Option(help="Maximum cluster mass.", show_default=True)
    ] = 1.0e15

    cluster_c: Annotated[
        float, typer.Option(help="Cluster concentration.", show_default=True)
    ] = 4.0

    ra_min: Annotated[
        float, typer.Option(help="Minimum right ascension.", show_default=True)
    ] = 12.14

    ra_max: Annotated[
        float, typer.Option(help="Maximum right ascension.", show_default=True)
    ] = 12.54

    dec_min: Annotated[
        float, typer.Option(help="Minimum declination.", show_default=True)
    ] = -55.323

    dec_max: Annotated[
        float, typer.Option(help="Maximum declination.", show_default=True)
    ] = -54.923

    z_min: Annotated[
        float, typer.Option(help="Minimum redshift.", show_default=True)
    ] = 0.01

    z_max: Annotated[
        float, typer.Option(help="Maximum redshift.", show_default=True)
    ] = 1.6

    z_dist: Annotated[
        GalaxyZDist,
        typer.Option(help="Galaxy redshift distribution.", show_default=True),
    ] = GalaxyZDist.GAUSS

    sigma_z: Annotated[
        float, typer.Option(help="Galaxy redshift dispersion.", show_default=True)
    ] = 0.03

    shape_dist: Annotated[
        GalaxySDShapeDist,
        typer.Option(help="Galaxy shape distribution.", show_default=True),
    ] = GalaxySDShapeDist.GAUSS

    galaxy_shape_e_rms: Annotated[
        float, typer.Option(help="Galaxy shape rms.", show_default=True)
    ] = 1.5e-1

    galaxy_shape_e_sigma: Annotated[
        float, typer.Option(help="Galaxy shape sigma.", show_default=True)
    ] = 1.0e-2

    parameter_list: Annotated[
        list[str],
        typer.Option(
            help="Parameter to fit.",
            show_default=True,
            default_factory=lambda: [
                "NcHaloMassSummary:log10MDelta",
                "NcHaloPosition:ra",
                "NcHaloPosition:dec",
            ],
        ),
    ]

    seed: Annotated[
        Optional[int], typer.Option(help="Random seed.", show_default=True)
    ] = None

    summary: Annotated[
        bool, typer.Option(help="Print experiment summary.", show_default=True)
    ] = False

    def __post_init__(self):
        """Generate LSST cluster weak lensing experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        exp = generate_lsst_cluster_wl(
            cluster_ra=self.cluster_ra,
            cluster_dec=self.cluster_dec,
            cluster_z=self.cluster_z,
            cluster_mass=self.cluster_mass,
            cluster_c=self.cluster_c,
            cluster_mass_min=self.cluster_mass_min,
            cluster_mass_max=self.cluster_mass_max,
            ra_min=self.ra_min,
            ra_max=self.ra_max,
            dec_min=self.dec_min,
            dec_max=self.dec_max,
            z_min=self.z_min,
            z_max=self.z_max,
            z_dist=self.z_dist,
            sigma_z=self.sigma_z,
            shape_dist=self.shape_dist,
            galaxy_shape_e_rms=self.galaxy_shape_e_rms,
            galaxy_shape_e_sigma=self.galaxy_shape_e_sigma,
            seed=self.seed,
            summary=self.summary,
        )

        mset = exp.peek("model-set")
        assert isinstance(mset, Ncm.MSet)

        likelihood = exp.peek("likelihood")
        assert isinstance(likelihood, Ncm.Likelihood)

        dataset = likelihood.peek_dataset()
        assert isinstance(dataset, Ncm.Dataset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

        for param in self.parameter_list:
            pindex = mset.param_get_by_full_name(param)
            if pindex is None:
                raise ValueError(f"Invalid parameter: {param}")
            mset.param_set_ftype(pindex.mid, pindex.pid, Ncm.ParamType.FREE)

        mset.prepare_fparam_map()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        ser.to_binfile(
            dataset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())
