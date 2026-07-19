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
from enum import StrEnum, auto
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
from numcosmo_py.experiments.cluster_richness_count import (
    generate_cluster_richness_count,
    load_cluster_richness_count,
)
from numcosmo_py.datasets.hicosmo import (
    SNIaID,
    BAOID,
    HID,
    add_bao_likelihood,
    add_h_likelihood,
    add_snia_likelihood,
)
from numcosmo_py.experiments.curvature_weight import (
    wspline_curvature_weight,
    qspline_curvature_weight,
)


class KnotPlacement(StrEnum):
    """Knot distribution for a spline reconstruction (q(z) or w(z)).

    DEFAULT keeps the model's own default (Chebyshev for the w-spline, uniform
    for the q-spline); UNIFORM / CHEBYSHEV force the placement explicitly.
    """

    DEFAULT = auto()
    UNIFORM = auto()
    CHEBYSHEV = auto()


def _spline_knots(placement: "KnotPlacement"):
    """Map a KnotPlacement to the Nc enum, or None to keep the model default."""
    if placement is KnotPlacement.UNIFORM:
        return Nc.HICosmoSplineKnots.UNIFORM
    if placement is KnotPlacement.CHEBYSHEV:
        return Nc.HICosmoSplineKnots.CHEBYSHEV
    return None


class CurvaturePriorType(StrEnum):
    """Curvature prior functional for a spline reconstruction (q(z) or w(z)).

    NONE disables the prior; MEAN_KAPPA is the geometric L2 curvature (default);
    LP_KAPPA / LP_D2 are the Lp norms of the geometric curvature and of the second
    derivative respectively, with order p given by ``--curvature-p`` (large p
    approaches the maximum curvature). LOCAL_KAPPA / LOCAL_D2 are the data-driven
    *local* counterparts: the curvature norm is weighted by W(x), built from the
    data Fisher information so the prior relaxes where the data constrains the
    reconstruction and stays strong where it is blind (see
    ``numcosmo_py.experiments.curvature_weight``).
    """

    NONE = auto()
    MEAN_KAPPA = auto()
    LP_KAPPA = auto()
    LP_D2 = auto()
    LOCAL_KAPPA = auto()
    LOCAL_D2 = auto()


def _add_curvature_prior(
    likelihood: Ncm.Likelihood,
    *,
    namespace: str,
    d2_name: str,
    prior_type: CurvaturePriorType,
    sigma: float,
    p: float,
    weight: "Ncm.Spline | None" = None,
) -> None:
    """Add the selected curvature Gaussian prior to ``likelihood``.

    The norm order p is fed through the PriorGaussFunc variable slot, which the
    nvar=1 ``lp_*``/``wlp_*`` functions read as ``x[0]``. ``namespace`` selects
    the model (e.g. ``NcHICosmoQSpline``); ``d2_name`` is its second-derivative
    function (``lp_q2`` or ``lp_w2``). For the LOCAL_* variants the weighted
    ``wlp_*`` function carries ``weight`` (a precomputed W(x) spline) as its
    associated object. The prior func is added to the likelihood only, not to the
    derived-function array (those are evaluated with nvar=0).
    """
    if prior_type is CurvaturePriorType.NONE:
        return

    obj: "Ncm.Spline | None" = None
    if prior_type is CurvaturePriorType.MEAN_KAPPA:
        func_name, var = f"{namespace}:mean_kappa", 0.0
    elif prior_type is CurvaturePriorType.LP_KAPPA:
        func_name, var = f"{namespace}:lp_kappa", p
    elif prior_type is CurvaturePriorType.LP_D2:
        func_name, var = f"{namespace}:{d2_name}", p
    elif prior_type is CurvaturePriorType.LOCAL_KAPPA:
        func_name, var, obj = f"{namespace}:wlp_kappa", p, weight
    else:  # LOCAL_D2
        func_name, var, obj = f"{namespace}:w{d2_name}", p, weight

    if prior_type in (CurvaturePriorType.LOCAL_KAPPA, CurvaturePriorType.LOCAL_D2):
        if weight is None:
            raise ValueError(f"{prior_type} requires a precomputed weight spline.")

    func = Ncm.MSetFuncList.new(func_name, obj)
    likelihood.priors_add(Ncm.PriorGaussFunc.new(func, 0.0, sigma, var))


def _add_function_grid(
    mfunc_oa: Ncm.ObjArray,
    func_name: str,
    z_nodes: "np.ndarray",
) -> None:
    """Append nvar=0 functions evaluating ``func_name`` at each redshift node.

    Each node binds its redshift through ``set_eval_x`` (which serializes with the
    function and is read by ``eval0``), so an MC/MCMC run records the function at
    every node as a catalog column. Sampling a whole grid yields the
    reconstruction band and its full node-to-node covariance for free.
    """
    for z in z_nodes:
        func = Ncm.MSetFuncList.new(func_name, None)
        func.set_eval_x([float(z)])
        mfunc_oa.add(func)


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

    cluster_redshift_sigma0: Annotated[
        float,
        typer.Option(
            help="Photo-z scatter sigma0 in sigma_z(z)=sigma0(1+z) (GAUSS only).",
            show_default=True,
            min=0.0,
        ),
    ] = 0.1

    lnM_obs_min: Annotated[
        float,
        typer.Option(
            help="Jpas minimum observed mass.",
            show_default=True,
        ),
    ] = (
        np.log(10.0) * 14.0
    )

    lnM_obs_max: Annotated[
        float,
        typer.Option(
            help="Jpas maximum observed mass.",
            show_default=True,
        ),
    ] = (
        np.log(10.0) * 15.0
    )

    lnMobsnknots: Annotated[
        int,
        typer.Option(
            help="Jpas number of observed mass bins.", show_default=True, min=2
        ),
    ] = 2

    cluster_mass_type: Annotated[
        Optional[ClusterMassType],
        typer.Option(help="Cluster mass-observable relation.", show_default=True),
    ] = ClusterMassType.NODIST

    survey_area: Annotated[
        float,
        typer.Option(
            help=(
                "Jpas survey area. This option is unavailable "
                "for the partial sky cases."
            ),
            show_default=True,
            min=0,
        ),
    ] = 2959.1

    resample_seed: Annotated[
        int,
        typer.Option(
            help=("Seed used to generate experiment."),
            show_default=True,
            min=0,
        ),
    ] = 1234

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
            cluster_redshift_sigma0=self.cluster_redshift_sigma0,
            lnM_obs_min=self.lnM_obs_min,
            lnM_obs_max=self.lnM_obs_max,
            lnMobsnknots=self.lnMobsnknots,
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
    ] = "composed variant=y1-source zp_min=0.0 zp_max=5.0 sigma0=0.03"

    shape_dist: Annotated[
        str,
        typer.Option(
            help=GalaxyShapeGen.get_help_text(),
            show_default=True,
            metavar=GalaxyShapeGen.get_help_metavar(),
            rich_help_panel="Galaxy shape source distribution",
        ),
    ] = "factor scheme=var_add ellip_conv=trace-det ellip_coord=celestial sigma=0.3 std_noise=0.1 c1_sigma=0.05 c2_sigma=0.05 m_sigma=0.05"

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

    use_lnint: Annotated[
        bool, typer.Option(help="Use log-integration.", show_default=True)
    ] = False

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

        if self.galaxy_density <= 0:
            raise typer.BadParameter(
                f"galaxy_density must be positive, got {self.galaxy_density}"
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
            use_lnint=self.use_lnint,
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
class GenerateClusterRichnessCount:
    """Generate a cluster mass-richness-count experiment.

    By default generates a mock: true log-masses and redshifts are drawn
    uniformly, and the observed richness (galaxy count) of each cluster is a
    Poisson-Lognormal realization implied by a NcClusterMassAscaso model,
    keeping only clusters above the selection cut.

    If --data-file is given, real data is loaded instead (e.g. true halo
    mass, redshift and galaxy count from a simulation catalog such as DC2,
    with no projection/background contamination applied); the mock
    generation options (n-clusters, mass-min/max, z-min/max, seed) are then
    ignored.
    """

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to write.")
    ]

    data_file: Annotated[
        Optional[Path],
        typer.Option(
            help=(
                "Path to a FITS file with real cluster data (mass, redshift "
                "and galaxy count columns). If given, real data is loaded "
                "instead of generating a mock."
            ),
            exists=True,
            readable=True,
        ),
    ] = None

    mass_column: Annotated[
        str,
        typer.Option(help="Name of the halo mass column (--data-file only)."),
    ] = "halo_mass"

    redshift_column: Annotated[
        str,
        typer.Option(help="Name of the redshift column (--data-file only)."),
    ] = "redshift"

    n_gal_column: Annotated[
        str,
        typer.Option(help="Name of the true galaxy count column (--data-file only)."),
    ] = "richness"

    mass_in_log10: Annotated[
        bool,
        typer.Option(
            help=("Mass column is log10(M) instead of linear M " "(--data-file only).")
        ),
    ] = False

    hdu: Annotated[int, typer.Option(help="HDU number to read from the FITS file.")] = 1

    n_clusters: Annotated[
        int,
        typer.Option(help="Number of mock clusters (before selection).", min=1),
    ] = 500

    mass_min: Annotated[
        float, typer.Option(help="Minimum halo mass.", show_default=True)
    ] = 1.0e13

    mass_max: Annotated[
        float, typer.Option(help="Maximum halo mass.", show_default=True)
    ] = 1.0e15

    z_min: Annotated[
        float, typer.Option(help="Minimum redshift.", show_default=True)
    ] = 0.1

    z_max: Annotated[
        float, typer.Option(help="Maximum redshift.", show_default=True)
    ] = 0.9

    mup0: Annotated[
        float, typer.Option(help="Mean ln-richness: constant term.", show_default=True)
    ] = 4.0

    mup1: Annotated[
        float,
        typer.Option(help="Mean ln-richness: ln-mass coefficient.", show_default=True),
    ] = 1.0

    mup2: Annotated[
        float,
        typer.Option(help="Mean ln-richness: ln(1+z) coefficient.", show_default=True),
    ] = 0.2

    sigmap0: Annotated[
        float,
        typer.Option(help="Ln-richness scatter: constant term.", show_default=True),
    ] = 0.5

    sigmap1: Annotated[
        float,
        typer.Option(
            help="Ln-richness scatter: ln-mass coefficient.", show_default=True
        ),
    ] = 0.03

    sigmap2: Annotated[
        float,
        typer.Option(
            help="Ln-richness scatter: ln(1+z) coefficient.", show_default=True
        ),
    ] = 0.15

    cut: Annotated[
        float,
        typer.Option(
            help=(
                "Richness selection cut in log units; clusters with "
                "N < round(exp(cut)) are discarded."
            ),
            show_default=True,
        ),
    ] = 0.0

    parameter_list: Annotated[
        list[str],
        typer.Option(
            help="Parameter to fit.",
            show_default=True,
            default_factory=lambda: [
                "NcClusterMass:mup0",
                "NcClusterMass:mup1",
                "NcClusterMass:mup2",
                "NcClusterMass:sigmap0",
                "NcClusterMass:sigmap1",
                "NcClusterMass:sigmap2",
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
        """Generate or load a cluster mass-richness-count experiment.

        Raises:
            ValueError: Invalid experiment file suffix, data column or fit
                parameter.
        """
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        if self.data_file is not None:
            exp = self._load_from_file()
        else:
            exp = generate_cluster_richness_count(
                n_clusters=self.n_clusters,
                lnM_min=float(np.log(self.mass_min)),
                lnM_max=float(np.log(self.mass_max)),
                z_min=self.z_min,
                z_max=self.z_max,
                mup0=self.mup0,
                mup1=self.mup1,
                mup2=self.mup2,
                sigmap0=self.sigmap0,
                sigmap1=self.sigmap1,
                sigmap2=self.sigmap2,
                cut=self.cut,
                seed=self.seed,
                summary=self.summary,
            )

        mset = exp.peek("model-set")
        assert isinstance(mset, Ncm.MSet)

        for param in self.parameter_list:
            pindex = mset.param_get_by_full_name(param)
            if pindex is None:
                raise ValueError(f"Invalid parameter: {param}")
            mset.param_set_ftype(pindex.mid, pindex.pid, Ncm.ParamType.FREE)

        mset.prepare_fparam_map()

        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        ser.dict_str_to_yaml_file(exp, self.experiment.absolute().as_posix())

    def _load_from_file(self) -> Ncm.ObjDictStr:
        """Load real cluster data from self.data_file into an experiment.

        Raises:
            ValueError: A requested column is not present in the table.
        """
        from astropy.table import Table  # pylint: disable=import-outside-toplevel

        assert self.data_file is not None
        table = Table.read(self.data_file, hdu=self.hdu)

        try:
            mass = np.array(table[self.mass_column], dtype=float)
            z = np.array(table[self.redshift_column], dtype=float)
            n_gal = np.array(table[self.n_gal_column], dtype=float)
        except KeyError as e:
            available_cols = ", ".join(table.colnames)
            raise ValueError(
                f"Column {e} not found. Available columns: {available_cols}"
            ) from e

        lnM = mass * np.log(10.0) if self.mass_in_log10 else np.log(mass)

        return load_cluster_richness_count(
            lnM=lnM,
            z=z,
            N=n_gal,
            mup0=self.mup0,
            mup1=self.mup1,
            mup2=self.mup2,
            sigmap0=self.sigmap0,
            sigmap1=self.sigmap1,
            sigmap2=self.sigmap2,
            cut=self.cut,
            summary=self.summary,
        )


@dataclasses.dataclass(kw_only=True)
class GenerateQSpline:
    """Generate QSpline experiment."""

    experiment: Annotated[
        Path, typer.Argument(help="Path to the experiment file to fit.")
    ]

    n_knots: Annotated[
        int, typer.Option(help="Number of knots.", show_default=True, min=6)
    ] = 6

    knots: Annotated[
        KnotPlacement,
        typer.Option(help="Knot placement in z.", show_default=True),
    ] = KnotPlacement.DEFAULT

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

    curvature_prior: Annotated[
        CurvaturePriorType,
        typer.Option(help="Curvature prior functional.", show_default=True),
    ] = CurvaturePriorType.MEAN_KAPPA

    curvature_sigma: Annotated[
        float,
        typer.Option(help="Curvature prior standard deviation.", show_default=True),
    ] = 3.0

    curvature_p: Annotated[
        float,
        typer.Option(
            help="Lp norm order p for the lp_* curvature priors "
            "(large p approaches the maximum curvature).",
            show_default=True,
        ),
    ] = 2.0

    curvature_ref_factor: Annotated[
        float,
        typer.Option(
            help="Crossover-scale knob for the LOCAL_* data-driven curvature "
            "weight (smaller relaxes the prior more eagerly where data informs).",
            show_default=True,
        ),
    ] = 1.0

    band_nodes: Annotated[
        int,
        typer.Option(
            help="Number of redshift nodes for the q(z) reconstruction band "
            "recorded in MC/MCMC catalogs (0 disables).",
            show_default=True,
            min=0,
        ),
    ] = 20

    def __post_init__(self):
        """Generate QSpline experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        knots = _spline_knots(self.knots)
        if knots is None:
            cosmo = Nc.HICosmoQSpline.new(
                Ncm.SplineCubicNotaknot(), self.n_knots, self.z_max
            )
        else:
            cosmo = Nc.HICosmoQSpline(
                spline=Ncm.SplineCubicNotaknot(),
                qparam_length=self.n_knots,
                zf=self.z_max,
                knots=knots,
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

        # Expose mean_kappa and the transition redshift as derived (nvar=0)
        # functions for post-processing.
        mfunc_oa = Ncm.ObjArray.new()
        mfunc_oa.add(Ncm.MSetFuncList.new("NcHICosmoQSpline:mean_kappa", None))
        mfunc_oa.add(Ncm.MSetFuncList.new("NcHICosmoQSpline:q_transition", None))

        # Sample q(z) on a redshift grid -> reconstruction band + covariance.
        if self.band_nodes > 0:
            _add_function_grid(
                mfunc_oa, "NcHICosmo:q", np.linspace(0.0, self.z_max, self.band_nodes)
            )

        weight = None
        if self.curvature_prior in (
            CurvaturePriorType.LOCAL_KAPPA,
            CurvaturePriorType.LOCAL_D2,
        ):
            weight = qspline_curvature_weight(
                dset,
                mset,
                cosmo,
                z_max=self.z_max,
                ref_factor=self.curvature_ref_factor,
            )

        _add_curvature_prior(
            likelihood,
            namespace="NcHICosmoQSpline",
            d2_name="lp_q2",
            prior_type=self.curvature_prior,
            sigma=self.curvature_sigma,
            p=self.curvature_p,
            weight=weight,
        )

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

        ser.array_to_yaml_file(
            mfunc_oa,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )


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

    knots: Annotated[
        KnotPlacement,
        typer.Option(help="Knot placement in alpha=ln(1+z).", show_default=True),
    ] = KnotPlacement.DEFAULT

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

    curvature_prior: Annotated[
        CurvaturePriorType,
        typer.Option(help="Curvature prior functional.", show_default=True),
    ] = CurvaturePriorType.MEAN_KAPPA

    curvature_sigma: Annotated[
        float,
        typer.Option(help="Curvature prior standard deviation.", show_default=True),
    ] = 3.0

    curvature_p: Annotated[
        float,
        typer.Option(
            help="Lp norm order p for the lp_* curvature priors "
            "(large p approaches the maximum curvature).",
            show_default=True,
        ),
    ] = 2.0

    curvature_ref_factor: Annotated[
        float,
        typer.Option(
            help="Crossover-scale knob for the LOCAL_* data-driven curvature "
            "weight (smaller relaxes the prior more eagerly where data informs).",
            show_default=True,
        ),
    ] = 1.0

    band_nodes: Annotated[
        int,
        typer.Option(
            help="Number of redshift nodes for the w(z) reconstruction band "
            "recorded in MC/MCMC catalogs (0 disables).",
            show_default=True,
            min=0,
        ),
    ] = 20

    fit_h0: Annotated[
        bool,
        typer.Option(
            help="Fit H0 even without Hubble data (e.g. with SH0ES-calibrated "
            "SNIa, which constrains H0).",
            show_default=True,
        ),
    ] = False

    def __post_init__(self):
        """Generate DE WSpline experiment."""
        Ncm.cfg_init()

        if self.experiment.suffix != ".yaml":
            raise ValueError(
                f"Invalid experiment file suffix: {self.experiment.suffix}"
            )

        knots = _spline_knots(self.knots)
        if knots is None:
            cosmo = Nc.HICosmoDEWSpline.new(self.n_knots, self.z_max)
        else:
            cosmo = Nc.HICosmoDEWSpline(
                zf=self.z_max, w_length=self.n_knots, knots=knots
            )
        cosmo.omega_x2omega_k()
        cosmo["Omegak"] = 0.0

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

        if self.include_hubble is not None:
            add_h_likelihood(dset, mset, self.include_hubble)
            cosmo.param_set_desc("H0", {"fit": True})

        if self.fit_h0:
            cosmo.param_set_desc("H0", {"fit": True})

        if dset.get_length() == 0:
            raise ValueError("No data included in the experiment.")

        mset.prepare_fparam_map()
        likelihood = Ncm.Likelihood.new(dset)
        # It can be useful to add Omega_m as a function when fitting Omegac
        mfunc_oa = Ncm.ObjArray.new()
        mfunc_Omegam = Ncm.MSetFuncList.new("NcHICosmo:Omega_m0", None)
        mfunc_oa.add(mfunc_Omegam)
        # Always expose mean_kappa as a derived (nvar=0) function for post-processing.
        mfunc_mean_kappa = Ncm.MSetFuncList.new("NcHICosmoDEWSpline:mean_kappa", None)
        mfunc_oa.add(mfunc_mean_kappa)

        # Sample w(z) on a redshift grid -> reconstruction band + covariance.
        if self.band_nodes > 0:
            _add_function_grid(
                mfunc_oa,
                "NcHICosmoDE:wDE_z",
                np.linspace(0.0, self.z_max, self.band_nodes),
            )

        weight = None
        if self.curvature_prior in (
            CurvaturePriorType.LOCAL_KAPPA,
            CurvaturePriorType.LOCAL_D2,
        ):
            weight = wspline_curvature_weight(
                dset, mset, cosmo, ref_factor=self.curvature_ref_factor
            )

        _add_curvature_prior(
            likelihood,
            namespace="NcHICosmoDEWSpline",
            d2_name="lp_w2",
            prior_type=self.curvature_prior,
            sigma=self.curvature_sigma,
            p=self.curvature_p,
            weight=weight,
        )

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

        ser.array_to_yaml_file(
            mfunc_oa,
            self.experiment.with_suffix(".functions.yaml").absolute().as_posix(),
        )
