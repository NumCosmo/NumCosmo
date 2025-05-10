#
# cluster_wl.py
#
# Tue Oct 22 12:50:10 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# cluster_wl.py
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

"""Factory functions for LSST Weak Lensing cluster forecast.

This module provides factory functions to create the experiment dictionary
for the LSS Weak Lensing cluster forecast.
"""

from typing import Annotated, Any
from enum import StrEnum, auto
from dataclasses import dataclass
import numpy as np
from pydantic import BaseModel, Field, ConfigDict, PrivateAttr
from pydantic_core import core_schema
from rich.console import Console
from rich.table import Table
from tabulate import tabulate

from numcosmo_py import Ncm, Nc, parse_options_strict, GEnum

DEFAULT_TABLE_FMT = "rounded_grid"

DEFAULT_SPEC_Z_MIN = 0.0
DEFAULT_SPEC_Z_MAX = 1.0


class EllipConv(GEnum):
    """Ellipticity convention."""

    # pylint: disable=no-member
    TRACE = Nc.GalaxyWLObsEllipConv.TRACE
    TRACE_DET = Nc.GalaxyWLObsEllipConv.TRACE_DET

    @classmethod
    def __get_pydantic_core_schema__(
        cls, _source_type: Any, _handler: Any
    ) -> core_schema.CoreSchema:
        """Get the Pydantic core schema for the TypeSource class."""
        return core_schema.no_info_before_validator_function(
            lambda v: cls(v) if isinstance(v, str) else v,
            core_schema.enum_schema(cls, list(cls), sub_type="str"),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda v: str(v.value)
            ),
        )


class GalaxyZGenSpec(BaseModel):
    """Galaxy redshift source distribution parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    z_min: Annotated[float, Field(ge=0.0)] = DEFAULT_SPEC_Z_MIN
    z_max: Annotated[float, Field(gt=0.0)] = DEFAULT_SPEC_Z_MAX

    _true: Nc.GalaxySDTrueRedshiftLSSTSRD = PrivateAttr()
    _spec: Nc.GalaxySDObsRedshiftSpec = PrivateAttr()

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyZGenSpec":
        """Create a GalaxyZDistSpec from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the galaxy redshift source distribution."""
        return [
            "GalaxyZDistSpec",
            f"z_min={DEFAULT_SPEC_Z_MIN}, z_max={DEFAULT_SPEC_Z_MAX}",
        ]

    def model_post_init(self, _: Any, /) -> None:
        """Check that z_min is less than z_max."""
        if self.z_min > self.z_max:
            raise ValueError(f"z_min {self.z_min} must be less than z_max {self.z_max}")

        self._true = Nc.GalaxySDTrueRedshiftLSSTSRD.new()
        self._spec = Nc.GalaxySDObsRedshiftSpec.new(self._true, self.z_min, self.z_max)

    def gen_z(
        self, mset: Ncm.MSet, z_data: Nc.GalaxySDObsRedshiftData, rng: Ncm.RNG
    ) -> bool:
        """Generate the galaxy redshift source distribution data observations."""
        return self._spec.gen1(mset, z_data, rng)

    def get_z_dist(self) -> Nc.GalaxySDObsRedshift:
        """Return the galaxy redshift source distribution data observations."""
        return self._spec


DEFAULT_GAUSS_ZP_MIN = 0.0
DEFAULT_GAUSS_ZP_MAX = 1.0
DEFAULT_GAUSS_SIGMA0 = 0.03


class GalaxyZGenGauss(BaseModel):
    """Galaxy redshift source distribution parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    zp_min: Annotated[float, Field(ge=0.0)] = DEFAULT_GAUSS_ZP_MIN
    zp_max: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_ZP_MAX
    sigma0: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_SIGMA0

    _true: Nc.GalaxySDTrueRedshiftLSSTSRD = PrivateAttr()
    _gauss: Nc.GalaxySDObsRedshiftGauss = PrivateAttr()

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyZGenGauss":
        """Create a GalaxyZDistGauss from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the galaxy redshift source distribution."""
        return [
            "GalaxyZDistGauss",
            (
                f"zp_min={DEFAULT_GAUSS_ZP_MIN}, "
                f"zp_max={DEFAULT_GAUSS_ZP_MAX}, "
                f"sigma0={DEFAULT_GAUSS_SIGMA0}"
            ),
        ]

    def model_post_init(self, _: Any, /) -> None:
        """Check that zp_min is less than zp_max."""
        if self.zp_min > self.zp_max:
            raise ValueError(
                f"zp_min {self.zp_min} must be less than zp_max {self.zp_max}"
            )

        self._true = Nc.GalaxySDTrueRedshiftLSSTSRD.new()
        self._gauss = Nc.GalaxySDObsRedshiftGauss.new(
            self._true, self.zp_min, self.zp_max
        )

    def gen_z(
        self, mset: Ncm.MSet, z_data: Nc.GalaxySDObsRedshiftData, rng: Ncm.RNG
    ) -> bool:
        """Generate the galaxy redshift source distribution data observations."""
        return self._gauss.gen1(mset, z_data, self.sigma0, rng)

    def get_z_dist(self) -> Nc.GalaxySDObsRedshift:
        """Return the galaxy redshift source distribution data observations."""
        return self._gauss


GalaxyZGenTypes = GalaxyZGenSpec | GalaxyZGenGauss


class GalaxyZGen(StrEnum):
    """Galaxy redshift source distribution types."""

    SPEC = (auto(), GalaxyZGenSpec)
    GAUSS = (auto(), GalaxyZGenGauss)

    def __new__(cls, value: str, _model_cls: type[GalaxyZGenTypes]) -> "GalaxyZGen":
        """Create a new instance of the enum.

        Initialize the enum including help text and model class.
        """
        obj = str.__new__(cls, value)
        obj._value_ = value
        return obj

    def __init__(self, _value: str, model_cls: type[GalaxyZGenTypes]) -> None:
        """Initialize the enum."""
        self.model_cls = model_cls

    @classmethod
    def get_help_text(cls) -> str:
        """Return the help text for the galaxy redshift source distribution."""
        headers = ["Enum", "Model", "Params"]
        rows = []

        for value in cls:
            rows.append([value] + value.model_cls.help_text())

        return "\b" + tabulate(rows, headers=headers, tablefmt=DEFAULT_TABLE_FMT)

    @classmethod
    def get_help_metavar(cls) -> str:
        """Return the help metavar for the galaxy redshift source distribution."""
        metavar = "[" + "|".join([value.value for value in cls]) + "]"
        return metavar


DEFAULT_SHAPE_GAUSS_SIGMA = 0.15
DEFAULT_SHAPE_GAUSS_STD_NOISE = 0.01


class GalaxyShapeGenGauss(BaseModel):
    """Galaxy shape Gaussian parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    ellip_conv: Annotated[EllipConv, Field()] = EllipConv.TRACE_DET
    sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_GAUSS_SIGMA
    std_noise: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_GAUSS_STD_NOISE

    _gauss: Nc.GalaxySDShapeGauss = PrivateAttr()

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the galaxy shape distribution."""
        return [
            "GalaxyShapeGauss",
            (
                f"ellip_conv={EllipConv.TRACE_DET.value}, "
                f"sigma={DEFAULT_SHAPE_GAUSS_SIGMA}, "
                f"std_noise={DEFAULT_SHAPE_GAUSS_STD_NOISE}"
            ),
        ]

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyShapeGenGauss":
        """Create a GalaxyShapeGauss from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def model_post_init(self, _: Any, /) -> None:
        """Check that sigma is less than std_noise."""
        self._gauss = Nc.GalaxySDShapeGauss.new(self.ellip_conv.genum)
        self._gauss["sigma"] = self.sigma

    def gen_shape(
        self, mset: Ncm.MSet, shape_data: Nc.GalaxySDShapeData, rng: Ncm.RNG
    ) -> None:
        """Generate the galaxy shape source distribution data observations."""
        coord = Nc.GalaxyWLObsCoord.EUCLIDEAN
        self._gauss.gen(mset, shape_data, self.std_noise, coord, rng)

    def get_shape_dist(self) -> Nc.GalaxySDShapeGauss:
        """Return the galaxy shape source distribution data observations."""
        return self._gauss

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


DEFAULT_GAUSS_HSC_STD_SHAPE = 0.15
DEFAULT_GAUSS_HSC_STD_NOISE = 0.01
DEFAULT_GAUSS_HSC_STD_SHAPE_SIGMA = 0.01
DEFAULT_GAUSS_HSC_C1_SIGMA = 0.01
DEFAULT_GAUSS_HSC_C2_SIGMA = 0.01
DEFAULT_GAUSS_HSC_M_SIGMA = 0.08


class GalaxyShapeGenGaussHSC(BaseModel):
    """Galaxy shape Gaussian parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    ellip_conv: Annotated[EllipConv, Field()] = EllipConv.TRACE_DET
    std_shape: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_STD_SHAPE
    std_noise: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_STD_NOISE
    sigma_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_STD_SHAPE_SIGMA
    c1_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_C1_SIGMA
    c2_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_C2_SIGMA
    m_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_GAUSS_HSC_M_SIGMA

    _gauss_hsc: Nc.GalaxySDShapeGaussHSC = PrivateAttr()
    _k: float = PrivateAttr()
    _theta: float = PrivateAttr()

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the galaxy shape distribution."""
        return [
            "GalaxyShapeGaussHSC",
            (
                f"ellip_conv={EllipConv.TRACE_DET.value}, "
                f"sigma={DEFAULT_GAUSS_HSC_STD_SHAPE}, "
                f"std_noise={DEFAULT_GAUSS_HSC_STD_NOISE}, "
                f"sigma_sigma={DEFAULT_GAUSS_HSC_STD_SHAPE_SIGMA}, "
                f"c1_sigma={DEFAULT_GAUSS_HSC_C1_SIGMA}, "
                f"c2_sigma={DEFAULT_GAUSS_HSC_C2_SIGMA}, "
                f"m_sigma={DEFAULT_GAUSS_HSC_M_SIGMA}"
            ),
        ]

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyShapeGenGaussHSC":
        """Create a GalaxyShapeGaussHSC from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def model_post_init(self, _, /):
        """Check that sigma is less than std_noise."""
        self._gauss_hsc = Nc.GalaxySDShapeGaussHSC.new()
        self._k = ((self.std_shape**2) / self.sigma_sigma) ** 2
        self._theta = self.sigma_sigma**2 / self.std_shape**2

    def gen_shape(self, mset: Ncm.MSet, shape_data: Nc.GalaxySDShapeData, rng: Ncm.RNG):
        """Generate the galaxy shape source distribution data observations."""
        std_noise1 = np.sqrt(rng.gamma_gen(self._k, self._theta))
        std_noise2 = np.sqrt(rng.gamma_gen(self._k, self._theta))
        c1 = rng.gaussian_gen(0.0, self.c1_sigma)
        c2 = rng.gaussian_gen(0.0, self.c2_sigma)
        m = np.exp(rng.gaussian_gen(0.0, self.m_sigma))
        coord = Nc.GalaxyWLObsCoord.EUCLIDEAN
        self._gauss_hsc.gen(
            mset, shape_data, std_noise1, std_noise2, c1, c2, m, coord, rng
        )

    def get_shape_dist(self) -> Nc.GalaxySDShapeGaussHSC:
        """Return the galaxy shape source distribution data observations."""
        return self._gauss_hsc

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


GalaxyShapeDistGenTypes = GalaxyShapeGenGauss | GalaxyShapeGenGaussHSC


class GalaxyShapeGen(StrEnum):
    """Galaxy shape source distribution types."""

    GAUSS = (auto(), GalaxyShapeGenGauss)
    GAUSS_HSC = (auto(), GalaxyShapeGenGaussHSC)

    def __new__(cls, value: str, _model_cls: type[GalaxyShapeDistGenTypes]):
        """Create a new instance of the enum.

        Initialize the enum including help text and model class.
        """
        obj = str.__new__(cls, value)
        obj._value_ = value
        return obj

    def __init__(self, _value: str, model_cls: type[GalaxyShapeDistGenTypes]):
        """Initialize the enum."""
        self.model_cls = model_cls

    @classmethod
    def get_help_text(cls):
        """Return the help text for the galaxy shape source distribution."""
        headers = ["Enum", "Model", "Params"]
        rows = []

        for value in cls:
            rows.append([value] + value.model_cls.help_text())

        return "\b" + tabulate(rows, headers=headers, tablefmt=DEFAULT_TABLE_FMT)

    @classmethod
    def get_help_metavar(cls):
        """Return the help metavar for the galaxy shape source distribution."""
        metavar = "[" + "|".join([value.value for value in cls]) + "]"
        return metavar


class HaloProfileType(StrEnum):
    """Halo density profile types."""

    NFW = auto()
    EINASTO = auto()
    HERNQUIST = auto()


@dataclass(frozen=True, kw_only=True)
class HaloPositionData:
    """Cluster position parameters."""

    ra: float
    dec: float
    z: float


@dataclass(frozen=True, kw_only=True)
class GalaxyDistributionData:
    """Galaxy distribution parameters.

    The galaxy distribution is defined by the minimum and maximum values of right
    ascension and declination, and the galaxy density. The galaxy density is given in
    deg^-2.
    """

    ra_min: float
    ra_max: float
    dec_min: float
    dec_max: float
    density: float = 18 * 60 * 60  # 18 arcmin^-2 in deg^-2


class ClusterModel:
    """Cluster model parameters."""

    def __init__(
        self,
        *,
        mass_def: Nc.HaloMassSummaryMassDef = Nc.HaloMassSummaryMassDef.CRITICAL,
        delta: float = 200.0,
        profile_type: HaloProfileType = HaloProfileType.NFW,
        position: HaloPositionData = HaloPositionData(ra=12.34, dec=-55.123, z=0.2),
        cluster_c: float = 4.0,
        r_min: float = 0.3,
        r_max: float = 3.0,
        cluster_mass: float = 1.0e14,
        dist: None | Nc.Distance = None,
    ) -> None:
        """Initialize the cluster model."""
        if dist is None:
            dist = Nc.Distance.new(5.0)

        self.mass_def = mass_def
        self.mass_delta = delta
        self.profile_type = profile_type
        self.r_min = r_min
        self.r_max = r_max

        halo_mass_summary = Nc.HaloCMParam.new(mass_def, delta)
        match profile_type:
            case HaloProfileType.NFW:
                self.density_profile = Nc.HaloDensityProfileNFW.new(halo_mass_summary)
            case HaloProfileType.EINASTO:
                self.density_profile = Nc.HaloDensityProfileEinasto.new(
                    halo_mass_summary
                )
            case HaloProfileType.HERNQUIST:
                self.density_profile = Nc.HaloDensityProfileHernquist.new(
                    halo_mass_summary
                )
            case _:
                raise ValueError(f"Invalid halo profile type: {profile_type}")

        self.halo_mass_summary = halo_mass_summary
        self.surface_mass_density = Nc.WLSurfaceMassDensity.new(dist)
        self.halo_position = Nc.HaloPosition.new(dist)

        self.halo_position["ra"] = position.ra
        self.halo_position["dec"] = position.dec
        self.halo_position["z"] = position.z

        self.halo_mass_summary["cDelta"] = cluster_c
        self.halo_mass_summary["log10MDelta"] = np.log10(cluster_mass)

    @property
    def position_data(self) -> HaloPositionData:
        """Return the cluster position."""
        return HaloPositionData(
            ra=self.halo_position["ra"],
            dec=self.halo_position["dec"],
            z=self.halo_position["z"],
        )

    @property
    def mass(self) -> float:
        """Return the cluster mass."""
        return 10.0 ** self.halo_mass_summary["log10MDelta"]

    @property
    def concentration(self) -> float:
        """Return the cluster concentration."""
        return self.halo_mass_summary["cDelta"]

    def prepare(self, cosmo: Nc.HICosmo) -> None:
        """Prepare the cluster model."""
        self.surface_mass_density.prepare(cosmo)
        self.halo_position.prepare(cosmo)


class GalaxyDistributionModel:
    """Galaxy distribution model parameters."""

    def __init__(
        self,
        galaxies: GalaxyDistributionData,
        z_gen: GalaxyZGenTypes = GalaxyZGenGauss(),
        shape_gen: GalaxyShapeDistGenTypes = GalaxyShapeGenGauss(),
    ) -> None:
        """Initialize the galaxy distribution model."""
        self.galaxy_position = Nc.GalaxySDPositionFlat.new(
            galaxies.ra_min, galaxies.ra_max, galaxies.dec_min, galaxies.dec_max
        )
        alpha_max = np.deg2rad(galaxies.ra_max)
        alpha_min = np.deg2rad(galaxies.ra_min)
        delta_max = np.deg2rad(galaxies.dec_max)
        delta_min = np.deg2rad(galaxies.dec_min)

        self.sky_area = (
            (alpha_max - alpha_min) * (np.sin(delta_max) - np.sin(delta_min))
        ) * (180.0 / np.pi) ** 2
        self.n_galaxies = int(galaxies.density * self.sky_area)
        self.z_gen = z_gen
        self.shape_gen = shape_gen

    def generate_data(
        self, cosmo: Nc.HICosmo, cluster: ClusterModel, rng: Ncm.RNG
    ) -> tuple[Nc.DataClusterWL, Ncm.MSet]:
        """Generate the galaxy data."""
        cluster.prepare(cosmo)
        cluster_data = Nc.DataClusterWL.new()
        galaxy_redshift = self.z_gen.get_z_dist()
        galaxy_shape = self.shape_gen.get_shape_dist()

        cluster_data.props.r_min = cluster.r_min
        cluster_data.props.r_max = cluster.r_max

        mset = Ncm.MSet.new_array(
            [
                cosmo,
                cluster.density_profile,
                cluster.surface_mass_density,
                cluster.halo_position,
                galaxy_redshift,
                self.galaxy_position,
                galaxy_shape,
            ]
        )

        cluster.halo_position.prepare_if_needed(cosmo)
        s_data_array = []
        for i in range(self.n_galaxies):
            z_data = Nc.GalaxySDObsRedshiftData.new(galaxy_redshift)
            p_data = Nc.GalaxySDPositionData.new(self.galaxy_position, z_data)
            s_data = Nc.GalaxySDShapeData.new(galaxy_shape, p_data)

            valid_z = self.z_gen.gen_z(mset, z_data, rng)
            self.galaxy_position.gen(mset, p_data, rng)
            theta, _ = cluster.halo_position.polar_angles(p_data.ra, p_data.dec)
            radius = cluster.halo_position.projected_radius(cosmo, theta)
            self.shape_gen.gen_shape(mset, s_data, rng)

            if (cluster.r_min <= radius <= cluster.r_max) and valid_z:
                s_data_array.append(s_data)

        self.n_galaxies = len(s_data_array)
        obs = Nc.GalaxyWLObs.new(
            galaxy_shape.get_ellip_conv(),
            Nc.GalaxyWLObsCoord.EUCLIDEAN,
            self.n_galaxies,
            list(s_data.required_columns()),
        )
        for i, s_data in enumerate(s_data_array):
            s_data.write_row(obs, i)

        cluster_data.set_obs(obs)

        return cluster_data, mset


def create_cosmo() -> Nc.HICosmo:
    """Create a cosmology for the cluster model."""
    cosmo = Nc.HICosmoDEXcdm()

    cosmo.params_set_default_ftype()
    cosmo.omega_x2omega_k()
    cosmo["H0"] = 67.81
    cosmo["Omegab"] = 0.0486
    cosmo["Omegac"] = 0.2612
    cosmo["w"] = -1.0
    cosmo["Omegak"] = 0.00

    prim = Nc.HIPrimPowerLaw.new()
    prim["ln10e10ASA"] = 3.02745
    prim["n_SA"] = 0.9660

    reion = Nc.HIReionCamb.new()

    cosmo.param_set_desc("H0", {"fit": False})
    cosmo.param_set_desc("Omegac", {"fit": False})
    cosmo.param_set_desc("Omegab", {"fit": False})
    cosmo.param_set_desc("w", {"fit": False})
    cosmo.param_set_desc("Omegak", {"fit": False})
    prim.param_set_desc("ln10e10ASA", {"fit": False})
    prim.param_set_desc("n_SA", {"fit": False})
    reion.param_set_desc("z_re", {"fit": False})

    cosmo.add_submodel(prim)
    cosmo.add_submodel(reion)

    return cosmo


def generate_lsst_cluster_wl(
    *,
    cluster_ra: float,
    cluster_dec: float,
    cluster_z: float,
    cluster_mass: float,
    cluster_mass_min: float,
    cluster_mass_max: float,
    cluster_c: float,
    r_min: float,
    r_max: float,
    ra_min: float,
    ra_max: float,
    dec_min: float,
    dec_max: float,
    z_gen: GalaxyZGenTypes,
    shape_gen: GalaxyShapeDistGenTypes,
    density: float,
    seed: None | int,
    summary: bool,
) -> Ncm.ObjDictStr:
    """Generate J-Pas forecast 2024 experiment dictionary."""
    if (
        (cluster_ra < ra_min)
        or (cluster_ra > ra_max)
        or (cluster_dec < dec_min)
        or (cluster_dec > dec_max)
    ):
        raise ValueError(f"Cluster position out of bounds: {cluster_ra} {cluster_dec}")

    halo_position = HaloPositionData(ra=cluster_ra, dec=cluster_dec, z=cluster_z)
    cluster = ClusterModel(
        position=halo_position,
        cluster_mass=cluster_mass,
        cluster_c=cluster_c,
        r_min=r_min,
        r_max=r_max,
    )
    galaxy_distribution = GalaxyDistributionData(
        ra_min=ra_min,
        ra_max=ra_max,
        dec_min=dec_min,
        dec_max=dec_max,
        density=density * 60 * 60,
    )
    galaxy_model = GalaxyDistributionModel(
        galaxies=galaxy_distribution,
        z_gen=z_gen,
        shape_gen=shape_gen,
    )
    if seed is not None:
        rng = Ncm.RNG.seeded_new("mt19937", seed)
    else:
        rng = Ncm.RNG.new("mt19937")
        rng.set_random_seed(True)

    cluster.halo_position.param_set_desc(
        "ra", {"lower-bound": ra_min, "upper-bound": ra_max}
    )
    cluster.halo_position.param_set_desc(
        "dec", {"lower-bound": dec_min, "upper-bound": dec_max}
    )
    cluster.halo_mass_summary.param_set_desc(
        "log10MDelta",
        {
            "lower-bound": float(np.log10(cluster_mass_min)),
            "upper-bound": float(np.log10(cluster_mass_max)),
        },
    )

    cluster_data, mset = galaxy_model.generate_data(create_cosmo(), cluster, rng)

    dset = Ncm.Dataset.new_array([cluster_data])
    likelihood = Ncm.Likelihood.new(dset)

    # Save experiment
    experiment = Ncm.ObjDictStr()

    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    if summary:
        # Here we use typer and rich to print the summary
        # We print the Halo parameters
        # Then the Galaxy parameters
        # The number of galaxies
        console = Console()

        table = Table(title="Cluster Parameters")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("RA", f"{halo_position.ra}")
        table.add_row("DEC", f"{halo_position.dec}")
        table.add_row("z", f"{halo_position.z}")
        table.add_row("Mass Profile", f"{cluster.profile_type}")
        table.add_row("Mass Definition", f"{cluster.mass_def.value_name}")
        table.add_row("Mass Delta", f"{cluster.mass_delta}")
        table.add_row("Mass", f"{cluster.mass:.2e}")
        table.add_row("Concentration", f"{cluster.concentration}")
        table.add_row("Minimum Radius", f"{cluster.r_min}")
        table.add_row("Maximum Radius", f"{cluster.r_max}")
        console.print(table)

        table = Table(title="Galaxy Sample Parameters")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("RA min", f"{galaxy_distribution.ra_min}")
        table.add_row("RA max", f"{galaxy_distribution.ra_max}")
        table.add_row("DEC min", f"{galaxy_distribution.dec_min}")
        table.add_row("DEC max", f"{galaxy_distribution.dec_max}")
        table.add_row("Density", f"{galaxy_distribution.density}")
        table.add_row("Redshift Generator", repr(z_gen))
        table.add_row("Shape Generator", repr(shape_gen))

        console.print(table)

        console.print(f"Number of galaxies: {galaxy_model.n_galaxies}")

    return experiment
