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
import shlex

import numpy as np
import typer

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck18 import (
    Planck18Types,
    HIPrimModel,
    generate_planck18_tt,
    generate_planck18_ttteee,
    mset_set_parameters,
)
from numcosmo_py.experiments.jpas_forecast24 import (
    ClusterRedshiftType,
    ClusterMassType,
    JpasSSCType,
    generate_jpas_forecast_2024,
)
from numcosmo_py.experiments.cluster_wl import (
    generate_lsst_cluster_wl,
    GalaxyShapeGen,
    GalaxyZGen,
)
from numcosmo_py.datasets.hicosmo import (
    SNIaID,
    BAOID,
    HID,
    add_bao_likelihood,
    add_h_likelihood,
    add_snia_likelihood,
)


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
        HIPrimModel,
        typer.Option(help="Primordial model to use.", show_default=True),
    ] = HIPrimModel.POWER_LAW

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
        """Generate Planck 2018 TT baseline experiment.

        Raises:
            ValueError: Invalid experiment file suffix.
            ValueError: Invalid data type.
        """
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
        mset_set_parameters(mset, self.data_type, self.prim_model)

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
            help=("Model for fitting. (NcHICosmo:Omegac,NcHICosmo:w,NcHICosmo:sigma8)"),
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

    def __post_init__(self):
        """Generate JPAS 2024 forecast experiment.

        Raises:
            ValueError: Invalid experiment file suffix.
        """
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

    r_min: Annotated[float, typer.Option(help="Minimum radius.", show_default=True)] = (
        0.3 / 0.7
    )

    r_max: Annotated[float, typer.Option(help="Maximum radius.", show_default=True)] = (
        3.0 / 0.7
    )

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

    profile_type: Annotated[
        str, typer.Option(help="Cluster profile to use.", show_default=True)
    ] = "nfw"

    z_dist: Annotated[
        str,
        typer.Option(
            help=GalaxyZGen.get_help_text(),
            show_default=True,
            metavar=GalaxyZGen.get_help_metavar(),
            rich_help_panel="Galaxy redshift source distribution",
        ),
    ] = "gauss zp_min=0.0 zp_max=5.0 sigma0=0.03"

    shape_dist: Annotated[
        str,
        typer.Option(
            help=GalaxyShapeGen.get_help_text(),
            show_default=True,
            metavar=GalaxyShapeGen.get_help_metavar(),
            rich_help_panel="Galaxy shape source distribution",
        ),
    ] = "gauss ellip_conv=trace-det ellip_coord=celestial sigma=0.3 std_noise=0.1"

    galaxy_density: Annotated[
        float, typer.Option(help="Galaxy density.", show_default=True)
    ] = 18.0

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
        """Generate LSST cluster weak lensing experiment.

        Raises:
            ValueError: Invalid experiment file suffix.
        """
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        _z_dist, z_gen = self.parse_z_dist()
        _shape_dist, shape_gen = self.parse_shape_dist()

        exp = generate_lsst_cluster_wl(
            cluster_ra=self.cluster_ra,
            cluster_dec=self.cluster_dec,
            cluster_z=self.cluster_z,
            cluster_mass=self.cluster_mass,
            cluster_c=self.cluster_c,
            r_min=self.r_min,
            r_max=self.r_max,
            cluster_mass_min=self.cluster_mass_min,
            cluster_mass_max=self.cluster_mass_max,
            ra_min=self.ra_min,
            ra_max=self.ra_max,
            dec_min=self.dec_min,
            dec_max=self.dec_max,
            profile_type=self.profile_type,
            z_gen=z_gen,
            shape_gen=shape_gen,
            density=self.galaxy_density,
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

        ser.to_binfile(
            dataset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())

    def parse_z_dist(self):
        """Parse the z_dist string."""
        z_dist_list = shlex.split(self.z_dist)
        z_dist_type = z_dist_list.pop(0)
        try:
            z_dist = GalaxyZGen(z_dist_type)
        except ValueError as e:
            raise typer.BadParameter(e)

        try:
            z_dist_args = z_dist.model_cls.from_args(z_dist_list)
        except ValueError as e:
            raise typer.BadParameter(e)
        return z_dist, z_dist_args

    def parse_shape_dist(self):
        """Parse the shape_dist string."""
        shape_dist_list = shlex.split(self.shape_dist)
        shape_dist_type = shape_dist_list.pop(0)
        try:
            shape_dist = GalaxyShapeGen(shape_dist_type)
        except ValueError as e:
            raise typer.BadParameter(e)

        try:
            shape_dist_args = shape_dist.model_cls.from_args(shape_dist_list)
        except ValueError as e:
            raise typer.BadParameter(e)
        return shape_dist, shape_dist_args


@dataclasses.dataclass(kw_only=True)
class GenerateQSpline:
    """Generate QSpline experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    n_knots: Annotated[
        int, typer.Option(help="Number of knots.", show_default=True, min=6)
    ] = 6

    z_max: Annotated[
        float, typer.Option(help="Maximum redshift.", show_default=True)
    ] = 2.1

    include_snia: Annotated[
        Optional[SNIaID], typer.Option(help="Include SNIa data.", show_default=True)
    ] = None

    include_bao: Annotated[
        Optional[BAOID], typer.Option(help="Include BAO data.", show_default=True)
    ] = None

    include_hubble: Annotated[
        Optional[HID], typer.Option(help="Include Hubble data.", show_default=True)
    ] = None

    def __post_init__(self):
        """Generate QSpline experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        cosmo = Nc.HICosmoQSpline.new(
            Ncm.SplineCubicNotaknot(), self.n_knots, self.z_max
        )
        for i in range(self.n_knots):
            cosmo.param_set_desc(f"qparam_{i}", {"fit": True})
        mset = Ncm.MSet.new_array([cosmo])
        dset = Ncm.Dataset.new()

        dist = Nc.Distance.new(10.0)

        if self.include_snia is not None:
            add_snia_likelihood(dset, mset, dist, self.include_snia)

        if self.include_bao is not None:
            add_bao_likelihood(dset, mset, dist, self.include_bao)
            cosmo.param_set_desc("asdrag", {"fit": True})

        if self.include_hubble is not None:
            add_h_likelihood(dset, mset, self.include_hubble)
            cosmo.param_set_desc("H0", {"fit": True})

        if dset.get_length() == 0:
            raise ValueError("No data included in the experiment.")

        mset.prepare_fparam_map()
        likelihood = Ncm.Likelihood.new(dset)
        # Save experiment
        experiment = Ncm.ObjDictStr()

        experiment.set("distance", dist)
        experiment.set("likelihood", likelihood)
        experiment.set("model-set", mset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.to_binfile(
            dset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(experiment, self.experiment.absolute().as_posix())


@dataclasses.dataclass(kw_only=True)
class GenerateXCDM:
    """Generate XCDM experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    curvature: Annotated[
        bool, typer.Option(help="Include curvature.", show_default=False)
    ] = False

    include_snia: Annotated[
        Optional[SNIaID], typer.Option(help="Include SNIa data.", show_default=True)
    ] = None

    include_bao: Annotated[
        Optional[BAOID], typer.Option(help="Include BAO data.", show_default=True)
    ] = None

    include_hubble: Annotated[
        Optional[HID], typer.Option(help="Include Hubble data.", show_default=True)
    ] = None

    def __post_init__(self):
        """Generate XCDM experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        cosmo = Nc.HICosmoDEXcdm.new()
        if self.curvature is False:
            cosmo.omega_x2omega_k()
            cosmo["Omegak"] = 0.0

        cosmo.param_set_desc("w", {"fit": True})
        cosmo.param_set_desc("Omegac", {"fit": True})
        # cosmo.param_set_desc(f"H0", {"fit": True})

        mset = Ncm.MSet.new_array([cosmo])
        dset = Ncm.Dataset.new()

        dist = Nc.Distance.new(10.0)

        if self.include_snia is not None:
            add_snia_likelihood(dset, mset, dist, self.include_snia)

        if self.include_bao is not None:
            add_bao_likelihood(dset, mset, dist, self.include_bao)

        if self.include_hubble is not None:
            add_h_likelihood(dset, mset, self.include_hubble)
            cosmo.param_set_desc("H0", {"fit": True})

        if dset.get_length() == 0:
            raise ValueError("No data included in the experiment.")

        mset.prepare_fparam_map()
        likelihood = Ncm.Likelihood.new(dset)
        # Save experiment
        experiment = Ncm.ObjDictStr()

        experiment.set("distance", dist)
        experiment.set("likelihood", likelihood)
        experiment.set("model-set", mset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.to_binfile(
            dset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(experiment, self.experiment.absolute().as_posix())

        mfunc_oa = Ncm.ObjArray.new()
        mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
        mfunc_oa.add(mfunc_Omegam)

        ser.array_to_yaml_file(
            mfunc_oa,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )


@dataclasses.dataclass(kw_only=True)
class GenerateDEWSpline:
    """Generate DE WSpline experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    n_knots: Annotated[
        int, typer.Option(help="Number of knots.", show_default=True, min=5)
    ] = 5

    z_max: Annotated[
        float, typer.Option(help="Maximum redshift.", show_default=True)
    ] = 2.33

    include_snia: Annotated[
        Optional[SNIaID], typer.Option(help="Include SNIa data.", show_default=True)
    ] = None

    include_bao: Annotated[
        Optional[BAOID], typer.Option(help="Include BAO data.", show_default=True)
    ] = None

    include_hubble: Annotated[
        Optional[HID], typer.Option(help="Include Hubble data.", show_default=True)
    ] = None

    def __post_init__(self):
        """Generate DE WSpline experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        cosmo = Nc.HICosmoDEWSpline.new(self.n_knots, self.z_max)
        for i in range(self.n_knots):
            cosmo.param_set_desc(f"w_{i}", {"fit": True})

        cosmo.param_set_desc("Omegac", {"fit": True})
        mset = Ncm.MSet.new_array([cosmo])
        dset = Ncm.Dataset.new()

        dist = Nc.Distance.new(10.0)

        if self.include_snia is not None:
            add_snia_likelihood(dset, mset, dist, self.include_snia)

        if self.include_bao is not None:
            add_bao_likelihood(dset, mset, dist, self.include_bao)
            # cosmo.param_set_desc("asdrag", {"fit": True})

        if self.include_hubble is not None:
            add_h_likelihood(dset, mset, self.include_hubble)
            cosmo.param_set_desc("H0", {"fit": True})

        if dset.get_length() == 0:
            raise ValueError("No data included in the experiment.")

        mset.prepare_fparam_map()
        likelihood = Ncm.Likelihood.new(dset)
        # Save experiment
        experiment = Ncm.ObjDictStr()

        experiment.set("distance", dist)
        experiment.set("likelihood", likelihood)
        experiment.set("model-set", mset)

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.to_binfile(
            dset, self.experiment.with_suffix(".dataset.gvar").absolute().as_posix()
        )
        ser.dict_str_to_yaml_file(experiment, self.experiment.absolute().as_posix())

        # It can be useful to add Omega_m as a function when fitting Omegac
        mfunc_oa = Ncm.ObjArray.new()
        mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
        mfunc_oa.add(mfunc_Omegam)

        ser.array_to_yaml_file(
            mfunc_oa,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )
