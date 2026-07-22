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

"""Factory functions for cluster Weak Lensing experiments.

This module builds ``NcDataClusterWLFactor`` experiments, on the
``NcGalaxy*Factor`` pipeline, in two ways:

- ``generate_lsst_cluster_wl``: a mock, following the LSST-SRD redshift
  distribution, with galaxy shapes and positions resampled from the chosen
  ``z_dist``/``shape_dist`` schemes.
- ``load_cluster_wl``: from a real ``NcGalaxyWLObs`` catalog (e.g. the
  curated Subaru HSC-SSP PDR1 catalogs shipped in the NumCosmo data file
  release, see ``NcGalaxyWLObsCatalogId``), pairing the observed galaxies
  (including their per-galaxy photometric redshift ``NcmSpline``) with a
  ``NcGalaxyRedshiftFactorSpline`` and the chosen ``shape_dist`` scheme.
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

# Keys used in a real NcGalaxyWLObs catalog's inherited NcmCatalog metadata
# (see ncm_catalog_peek_meta()) to carry the lensing cluster's own position,
# redshift and concentration -- information the galaxy catalog itself cannot
# provide, since it only describes the observed background galaxies.
CATALOG_META_CLUSTER_RA = "cluster_ra"
CATALOG_META_CLUSTER_DEC = "cluster_dec"
CATALOG_META_CLUSTER_Z = "cluster_z"
CATALOG_META_CLUSTER_C = "cluster_c"


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


class EllipCoord(GEnum):
    """Ellipticity coordinate system."""

    # pylint: disable=no-member
    CARTESIAN = Nc.WLEllipticityFrame.CARTESIAN
    CELESTIAL = Nc.WLEllipticityFrame.CELESTIAL

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


class MassDef(GEnum):
    """Mass definition for the halo mass summary."""

    # pylint: disable=no-member
    CRITICAL = Nc.HaloMassSummaryMassDef.CRITICAL
    MEAN = Nc.HaloMassSummaryMassDef.MEAN
    VIRIAL = Nc.HaloMassSummaryMassDef.VIRIAL

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


class LSSTVariant(GEnum):
    """LSST-SRD redshift-distribution variant."""

    # pylint: disable=no-member
    Y1_SOURCE = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_SOURCE
    Y1_LENS = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_LENS
    Y10_SOURCE = Nc.GalaxyRedshiftPopLSSTSRDType.Y10_SOURCE
    Y10_LENS = Nc.GalaxyRedshiftPopLSSTSRDType.Y10_LENS

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


class HWLCatalogID(GEnum):
    """Curated Subaru HSC-SSP PDR1 weak lensing catalog (NumCosmo data release)."""

    # pylint: disable=no-member
    HWL16A_002 = Nc.GalaxyWLObsCatalogId(0)
    HWL16A_007 = Nc.GalaxyWLObsCatalogId(1)
    HWL16A_060 = Nc.GalaxyWLObsCatalogId(2)
    HWL16A_064 = Nc.GalaxyWLObsCatalogId(3)
    HWL16A_094 = Nc.GalaxyWLObsCatalogId(4)

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


DEFAULT_COMPOSED_ZP_MIN = 0.0
DEFAULT_COMPOSED_ZP_MAX = 5.0
DEFAULT_COMPOSED_SIGMA0 = 0.03


class GalaxyZGenComposed(BaseModel):
    """Galaxy redshift source distribution parameters.

    Builds a ``NcGalaxyRedshiftFactorComposed`` scheme: an LSST-SRD true-z
    population (``NcGalaxyRedshiftPopLSSTSRD``) observed through a Gaussian
    photo-z error kernel (``NcGalaxyRedshiftObsGauss``), parameterised by a
    per-galaxy ``sigma0`` (the photo-z error scales as ``sigma0*(1+z)``).
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    variant: Annotated[LSSTVariant, Field()] = LSSTVariant.Y1_SOURCE
    zp_min: Annotated[float, Field(ge=0.0)] = DEFAULT_COMPOSED_ZP_MIN
    zp_max: Annotated[float, Field(gt=0.0)] = DEFAULT_COMPOSED_ZP_MAX
    sigma0: Annotated[float, Field(gt=0.0)] = DEFAULT_COMPOSED_SIGMA0

    _pop: Nc.GalaxyRedshiftPopLSSTSRD = PrivateAttr()
    _obs_sel: Nc.GalaxyRedshiftObsGauss = PrivateAttr()
    _factor: Nc.GalaxyRedshiftFactorComposed = PrivateAttr()

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyZGenComposed":
        """Create a GalaxyZGenComposed from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the galaxy redshift source distribution."""
        return [
            "GalaxyZGenComposed",
            (
                f"variant={LSSTVariant.Y1_SOURCE.value}, "
                f"zp_min={DEFAULT_COMPOSED_ZP_MIN}, "
                f"zp_max={DEFAULT_COMPOSED_ZP_MAX}, "
                f"sigma0={DEFAULT_COMPOSED_SIGMA0}"
            ),
        ]

    def model_post_init(self, _: Any, /) -> None:
        """Check that zp_min is less than zp_max and build the models."""
        if self.zp_min > self.zp_max:
            raise ValueError(
                f"zp_min {self.zp_min} must be less than zp_max {self.zp_max}"
            )

        self._pop = Nc.GalaxyRedshiftPopLSSTSRD.new_from_type(self.variant.genum)
        self._obs_sel = Nc.GalaxyRedshiftObsGauss.new()
        self._factor = Nc.GalaxyRedshiftFactorComposed.new(self.zp_min, self.zp_max)

    def register_models(self, mset: Ncm.MSet) -> None:
        """Register this scheme's population/observable models into mset."""
        mset.set(self._pop)
        mset.set(self._obs_sel)

    def get_redshift_factor(self) -> Nc.GalaxyRedshiftFactor:
        """Return the redshift factor for the Data orchestrator."""
        return self._factor

    def write_calib(self, obs: Nc.GalaxyWLObs, i: int) -> None:
        """Write this galaxy's fixed calibration input (sigma0) into obs."""
        obs.set(Nc.GALAXY_REDSHIFT_OBS_GAUSS_COL_SIGMA0, i, self.sigma0)

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


GalaxyZGenTypes = GalaxyZGenComposed


class GalaxyZGen(StrEnum):
    """Galaxy redshift source distribution types."""

    COMPOSED = (auto(), GalaxyZGenComposed)

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


DEFAULT_SHAPE_STD_NOISE = 0.1
DEFAULT_SHAPE_STD_SIGMA = 0.01
DEFAULT_SHAPE_C1_SIGMA = 0.01
DEFAULT_SHAPE_C2_SIGMA = 0.01
DEFAULT_SHAPE_M_SIGMA = 0.08


DEFAULT_SHAPE_SERIES_TRUNC_ORDER = 4


class GalaxyShapeFactorGenBase(BaseModel):
    """Shared measurement-noise/calibration parameters for every
    ``NcGalaxyShapeFactor*`` marginalization scheme.

    Per-galaxy measurement noise and additive/multiplicative calibration
    biases are drawn once (Gamma-distributed noise around a target mean,
    Gaussian calibration biases), then written as fixed per-galaxy inputs
    ahead of ``NcDataClusterWLFactor.resample()``. The intrinsic-ellipticity
    population is a separate, orthogonal choice -- see ``--pop-dist``/
    ``GalaxyPopGen``, and ``check_shape_pop_compat()`` for the one
    scheme/population constraint that crosses the two.

    Each concrete scheme is its own subclass (own ``--shape-factor`` tag, own
    ``help_text()``) so a scheme-specific option (e.g. SeriesLensed's
    ``series_trunc_order``) is only ever shown/accepted for that scheme --
    not a dead option under every other one.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    ellip_conv: Annotated[EllipConv, Field()] = EllipConv.TRACE_DET
    ellip_coord: Annotated[EllipCoord, Field()] = EllipCoord.CELESTIAL
    std_noise: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_STD_NOISE
    std_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_STD_SIGMA
    c1_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_C1_SIGMA
    c2_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_C2_SIGMA
    m_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_SHAPE_M_SIGMA

    _factor: Nc.GalaxyShapeFactor = PrivateAttr()
    _k_noise: float = PrivateAttr()
    _theta_noise: float = PrivateAttr()

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyShapeFactorGenBase":
        """Create this shape factor scheme from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        """Build the concrete NcGalaxyShapeFactor*. Overridden by each scheme."""
        raise NotImplementedError  # pragma: no cover

    def requires_sigma(self) -> bool:
        """Whether this scheme requires nc_galaxy_shape_pop_get_sigma().

        True (the default) for schemes that linearize around a Gaussian --
        overridden to False by schemes that don't (see
        ``check_shape_pop_compat()``).
        """
        return True

    def model_post_init(self, _: Any, /) -> None:
        """Build the shape-marginalization factor."""
        self._factor = self._build_shape_factor()
        self._k_noise = ((self.std_noise**2) / self.std_sigma) ** 2
        self._theta_noise = self.std_sigma**2 / self.std_noise**2

    def get_shape_factor(self) -> Nc.GalaxyShapeFactor:
        """Return the shape factor for the Data orchestrator."""
        return self._factor

    def write_calib(self, obs: Nc.GalaxyWLObs, i: int, rng: Ncm.RNG) -> None:
        """Draw and write this galaxy's fixed measurement-calibration inputs into obs."""
        std_noise = min(np.sqrt(rng.gamma_gen(self._k_noise, self._theta_noise)), 0.5)
        c1 = rng.gaussian_gen(0.0, self.c1_sigma)
        c2 = rng.gaussian_gen(0.0, self.c2_sigma)
        m = np.exp(rng.gaussian_gen(0.0, self.m_sigma))

        obs.set("std_noise", i, std_noise)
        obs.set("c1", i, c1)
        obs.set("c2", i, c2)
        obs.set("m", i, m)

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


_SHARED_SHAPE_FACTOR_HELP = (
    f"ellip_conv={EllipConv.TRACE_DET.value}, "
    f"ellip_coord={EllipCoord.CELESTIAL.value}, "
    f"std_noise={DEFAULT_SHAPE_STD_NOISE}, "
    f"std_sigma={DEFAULT_SHAPE_STD_SIGMA}, \n"
    f"c1_sigma={DEFAULT_SHAPE_C1_SIGMA}, "
    f"c2_sigma={DEFAULT_SHAPE_C2_SIGMA}, "
    f"m_sigma={DEFAULT_SHAPE_M_SIGMA}"
)


class GalaxyShapeFactorGenVarAdd(GalaxyShapeFactorGenBase):
    """Variance-add linearization (``NcGalaxyShapeFactorVarAdd``)."""

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for this scheme."""
        return ["GalaxyShapeFactorGenVarAdd", _SHARED_SHAPE_FACTOR_HELP]

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        return Nc.GalaxyShapeFactorVarAdd.new(self.ellip_conv.genum)


class GalaxyShapeFactorGenSeriesLensed(GalaxyShapeFactorGenBase):
    """g-Taylor series around the lens (``NcGalaxyShapeFactorSeriesLensed``)."""

    series_trunc_order: Annotated[int, Field(gt=0)] = DEFAULT_SHAPE_SERIES_TRUNC_ORDER

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for this scheme."""
        return [
            "GalaxyShapeFactorGenSeriesLensed",
            f"{_SHARED_SHAPE_FACTOR_HELP}, \n"
            f"series_trunc_order={DEFAULT_SHAPE_SERIES_TRUNC_ORDER}",
        ]

    def requires_sigma(self) -> bool:
        """SeriesLensed works with any population (no Gaussian linearization)."""
        return False

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        return Nc.GalaxyShapeFactorSeriesLensed.new(
            self.ellip_conv.genum, self.series_trunc_order
        )


class GalaxyShapeFactorGenQuad(GalaxyShapeFactorGenBase):
    """Direct source-domain quadrature (``NcGalaxyShapeFactorQuad``)."""

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for this scheme."""
        return ["GalaxyShapeFactorGenQuad", _SHARED_SHAPE_FACTOR_HELP]

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        return Nc.GalaxyShapeFactorQuad.new(self.ellip_conv.genum)


class GalaxyShapeFactorGenFixedQuad(GalaxyShapeFactorGenBase):
    """Fixed lens-domain quadrature (``NcGalaxyShapeFactorFixedQuad``)."""

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for this scheme."""
        return ["GalaxyShapeFactorGenFixedQuad", _SHARED_SHAPE_FACTOR_HELP]

    def requires_sigma(self) -> bool:
        """FixedQuad works with any population (no Gaussian linearization)."""
        return False

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        return Nc.GalaxyShapeFactorFixedQuad.new(self.ellip_conv.genum)


class GalaxyShapeFactorGenLaplace(GalaxyShapeFactorGenBase):
    """Laplace-approximation marginal (``NcGalaxyShapeFactorLaplace``)."""

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for this scheme."""
        return ["GalaxyShapeFactorGenLaplace", _SHARED_SHAPE_FACTOR_HELP]

    def _build_shape_factor(self) -> Nc.GalaxyShapeFactor:
        return Nc.GalaxyShapeFactorLaplace.new(self.ellip_conv.genum)


GalaxyShapeFactorGenTypes = (
    GalaxyShapeFactorGenVarAdd
    | GalaxyShapeFactorGenSeriesLensed
    | GalaxyShapeFactorGenQuad
    | GalaxyShapeFactorGenFixedQuad
    | GalaxyShapeFactorGenLaplace
)


class ShapeFactorGen(StrEnum):
    """Galaxy shape-marginalization factor types (``NcGalaxyShapeFactor*``)."""

    VAR_ADD = (auto(), GalaxyShapeFactorGenVarAdd)
    SERIES_LENSED = (auto(), GalaxyShapeFactorGenSeriesLensed)
    QUAD = (auto(), GalaxyShapeFactorGenQuad)
    FIXED_QUAD = (auto(), GalaxyShapeFactorGenFixedQuad)
    LAPLACE = (auto(), GalaxyShapeFactorGenLaplace)

    def __new__(
        cls, value: str, _model_cls: type[GalaxyShapeFactorGenTypes]
    ) -> "ShapeFactorGen":
        """Create a new instance of the enum.

        Initialize the enum including help text and model class.
        """
        obj = str.__new__(cls, value)
        obj._value_ = value
        return obj

    def __init__(self, _value: str, model_cls: type[GalaxyShapeFactorGenTypes]) -> None:
        """Initialize the enum."""
        self.model_cls = model_cls

    @classmethod
    def get_help_text(cls) -> str:
        """Return the help text for the galaxy shape-marginalization factor."""
        headers = ["Enum", "Model", "Params"]
        rows = []

        for value in cls:
            rows.append([value] + value.model_cls.help_text())

        return "\b" + tabulate(rows, headers=headers, tablefmt=DEFAULT_TABLE_FMT)

    @classmethod
    def get_help_metavar(cls) -> str:
        """Return the help metavar for the galaxy shape-marginalization factor."""
        metavar = "[" + "|".join([value.value for value in cls]) + "]"
        return metavar


DEFAULT_POP_SIGMA = 0.3
DEFAULT_POP_STD_SIGMA = 0.01
DEFAULT_POP_BETA_ALPHA = Nc.GALAXY_SHAPE_POP_BETA_DEFAULT_ALPHA
DEFAULT_POP_BETA_BETA = Nc.GALAXY_SHAPE_POP_BETA_DEFAULT_BETA


class GalaxyPopGenGauss(BaseModel):
    """Global Gaussian intrinsic-ellipticity population (``NcGalaxyShapePopGauss``)."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_POP_SIGMA

    _pop: Nc.GalaxyShapePopGauss = PrivateAttr()

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the Gaussian population."""
        return ["GalaxyPopGenGauss", f"sigma={DEFAULT_POP_SIGMA}"]

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyPopGenGauss":
        """Create a GalaxyPopGenGauss from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def model_post_init(self, _: Any, /) -> None:
        """Build the population model."""
        self._pop = Nc.GalaxyShapePopGauss.new()
        self._pop["sigma"] = self.sigma

    def register_models(self, mset: Ncm.MSet) -> None:
        """Register this population model into mset."""
        mset.set(self._pop)

    def supports_sigma(self) -> bool:
        """Whether this population supports nc_galaxy_shape_pop_get_sigma()."""
        return True

    def get_mfuncs(self) -> list[Ncm.MSetFuncList]:
        """Derived NcmMSetFuncList entries worth recording in an MC/MCMC catalog.

        None here: sigma is already a directly fitted model parameter, so
        there is nothing else meaningful to derive from it.
        """
        return []

    def write_calib(self, obs: Nc.GalaxyWLObs, i: int, rng: Ncm.RNG) -> None:
        """No per-galaxy calibration input: sigma is a single global parameter."""

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


class GalaxyPopGenGaussLocal(BaseModel):
    """Per-galaxy Gaussian population (``NcGalaxyShapePopGaussLocal``).

    Reads a per-galaxy ``e_rms`` catalog column instead of a global model
    parameter. For mock generation, ``e_rms`` is drawn the same way
    ``GalaxyShapeFactorGenBase``'s own ``std_noise`` is (Gamma-distributed
    around ``sigma`` with spread ``std_sigma``, both here); for real
    catalogs the catalog itself must carry an ``e_rms`` column, or loading
    fails with a clear error.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_POP_SIGMA
    std_sigma: Annotated[float, Field(gt=0.0)] = DEFAULT_POP_STD_SIGMA

    _pop: Nc.GalaxyShapePopGaussLocal = PrivateAttr()
    _k_e_rms: float = PrivateAttr()
    _theta_e_rms: float = PrivateAttr()

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the per-galaxy Gaussian population."""
        return [
            "GalaxyPopGenGaussLocal",
            f"sigma={DEFAULT_POP_SIGMA}, std_sigma={DEFAULT_POP_STD_SIGMA}",
        ]

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyPopGenGaussLocal":
        """Create a GalaxyPopGenGaussLocal from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def model_post_init(self, _: Any, /) -> None:
        """Build the population model and the mock e_rms sampling distribution."""
        self._pop = Nc.GalaxyShapePopGaussLocal.new()
        # Gamma-on-variance reparametrization, same as GalaxyShapeFactorGenBase's
        # own std_noise: mean sigma**2, spread controlled by std_sigma.
        self._k_e_rms = ((self.sigma**2) / self.std_sigma) ** 2
        self._theta_e_rms = self.std_sigma**2 / self.sigma**2

    def register_models(self, mset: Ncm.MSet) -> None:
        """Register this population model into mset."""
        mset.set(self._pop)

    def supports_sigma(self) -> bool:
        """Whether this population supports nc_galaxy_shape_pop_get_sigma()."""
        return True

    def get_mfuncs(self) -> list[Ncm.MSetFuncList]:
        """Derived NcmMSetFuncList entries worth recording in an MC/MCMC catalog.

        None here: this population has no global model parameters at all
        (e_rms is a fixed per-galaxy catalog input), so there is nothing to
        derive.
        """
        return []

    def write_calib(self, obs: Nc.GalaxyWLObs, i: int, rng: Ncm.RNG) -> None:
        """Draw and write this galaxy's fixed e_rms input into obs."""
        # e_rms = sqrt(<x>/2) with x restricted to [0, 1] saturates at 0.5 as
        # sigma -> infinity (nc_galaxy_shape_pop_gauss_local.c's own
        # _e_rms_of_sigma()/_sigma_from_e_rms()); NOT 1.0 -- an e_rms draw at
        # or above 0.5 has no finite sigma to invert to and aborts the whole
        # process via that function's own g_assert_cmpfloat(). Clip well
        # inside the true (0.5 - 1e-9) ceiling, not at it, since sigma
        # diverges as e_rms -> 0.5 (bisection there is numerically fragile,
        # not just technically in-bounds).
        e_rms = min(np.sqrt(rng.gamma_gen(self._k_e_rms, self._theta_e_rms)), 0.49)
        obs.set(Nc.GALAXY_SHAPE_POP_GAUSS_LOCAL_COL_E_RMS, i, e_rms)

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


class GalaxyPopGenBeta(BaseModel):
    """Beta intrinsic-ellipticity population (``NcGalaxyShapePopBeta``).

    Global model over $x=|\\chi_I|^2$, parameterized directly by the Beta
    distribution's own shape parameters ``alpha``/``beta``. ``beta`` is
    bounded ``>=1`` by the C model itself (the density would otherwise
    diverge at $x=1$); ``alpha``'s own floor is looser (``>=0.5001``) since
    an ``alpha<1`` divergence at $x=0$ is only a real problem for shape
    schemes that Taylor-expand in the shear (SeriesLensed), not for
    FixedQuad's direct quadrature -- see ``NcGalaxyShapePopBeta``'s own
    class documentation for the full rationale (in particular: don't let a
    hard floor sit exactly where a fit's own posterior wants its mass).
    Incompatible with shape schemes that linearize around a Gaussian (see
    ``check_shape_pop_compat()``).

    ``mean``/``concentration`` (``alpha/(alpha+beta)``/``alpha+beta``) are
    reporting quantities, not fields here -- read them off a built instance
    via ``get_shape_pop().get_mean()``/``.get_concentration()``, or from an
    ``NcmMSetFuncList`` entry (``NcGalaxyShapePopBeta:mean``/
    ``:concentration``) once the population is registered into an ``NcmMSet``.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    # Lower bounds only -- loose sanity checks matching this class's own
    # rationale (see the docstring above), not a mirror of the C model's
    # advisory fit-bounds (which a diagnostic run may legitimately want to
    # exceed, e.g. testing a concentrated population past its own declared
    # range -- see test_galaxy_shape_pop.py's own overflow-robustness
    # tests). Must track nc_galaxy_shape_pop_beta.c's own alpha/beta floors.
    alpha: Annotated[float, Field(ge=0.5001)] = DEFAULT_POP_BETA_ALPHA
    beta: Annotated[float, Field(ge=1.0)] = DEFAULT_POP_BETA_BETA

    _pop: Nc.GalaxyShapePopBeta = PrivateAttr()

    @staticmethod
    def help_text() -> list[str]:
        """Return the help text for the Beta population."""
        return [
            "GalaxyPopGenBeta",
            f"alpha={DEFAULT_POP_BETA_ALPHA}, beta={DEFAULT_POP_BETA_BETA}",
        ]

    @classmethod
    def from_args(cls, args: list[str]) -> "GalaxyPopGenBeta":
        """Create a GalaxyPopGenBeta from command line arguments."""
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    def model_post_init(self, _: Any, /) -> None:
        """Build the population model."""
        self._pop = Nc.GalaxyShapePopBeta.new()
        self._pop["alpha"] = self.alpha
        self._pop["beta"] = self.beta

    def get_shape_pop(self) -> Nc.GalaxyShapePopBeta:
        """Return the built NcGalaxyShapePopBeta (e.g. for get_mean()/get_concentration())."""
        return self._pop

    def register_models(self, mset: Ncm.MSet) -> None:
        """Register this population model into mset."""
        mset.set(self._pop)

    def supports_sigma(self) -> bool:
        """Whether this population supports nc_galaxy_shape_pop_get_sigma()."""
        return False

    def get_mfuncs(self) -> list[Ncm.MSetFuncList]:
        """Derived NcmMSetFuncList entries worth recording in an MC/MCMC catalog.

        alpha/beta are the fitted parameters, but mean/std of x = |chi_I|^2
        are the physically meaningful reported quantities -- recorded as
        catalog columns via NcmFitESMCMC.new_funcs_array()/an MC run's own
        equivalent, the same convention every other generate command in this
        CLI uses for its own derived quantities (see e.g. GenerateQSpline's
        mean_kappa/q_transition).
        """
        return [
            Ncm.MSetFuncList.new("NcGalaxyShapePopBeta:mean", None),
            Ncm.MSetFuncList.new("NcGalaxyShapePopBeta:std", None),
        ]

    def write_calib(self, obs: Nc.GalaxyWLObs, i: int, rng: Ncm.RNG) -> None:
        """No per-galaxy calibration input: alpha/beta are single global parameters."""

    def __repr__(self):
        """Return a string representation of the model."""
        args = ", ".join(f"{k}={v!r}" for k, v in self.model_dump().items())
        return f"{self.__class__.__name__}({args})"


GalaxyPopGenTypes = GalaxyPopGenGauss | GalaxyPopGenGaussLocal | GalaxyPopGenBeta


def check_shape_pop_compat(
    shape_gen: GalaxyShapeFactorGenTypes, pop_gen: GalaxyPopGenTypes
) -> None:
    """Check that a shape-marginalization factor and a population are compatible.

    Raises:
        ValueError: shape_gen requires nc_galaxy_shape_pop_get_sigma() and
            pop_gen doesn't support it (e.g. --shape-factor var_add with
            --pop-dist beta) -- letting this reach the C layer instead
            aborts the whole process, not just this check.
    """
    if shape_gen.requires_sigma() and not pop_gen.supports_sigma():
        raise ValueError(
            f"{type(shape_gen).__name__} requires a population supporting "
            f"nc_galaxy_shape_pop_get_sigma(); got {type(pop_gen).__name__}, "
            "which doesn't (pop_dist=gauss or gauss_local do)."
        )


class GalaxyPopGen(StrEnum):
    """Galaxy intrinsic-ellipticity population distribution types."""

    GAUSS = (auto(), GalaxyPopGenGauss)
    GAUSS_LOCAL = (auto(), GalaxyPopGenGaussLocal)
    BETA = (auto(), GalaxyPopGenBeta)

    def __new__(cls, value: str, _model_cls: type[GalaxyPopGenTypes]) -> "GalaxyPopGen":
        """Create a new instance of the enum.

        Initialize the enum including help text and model class.
        """
        obj = str.__new__(cls, value)
        obj._value_ = value
        return obj

    def __init__(self, _value: str, model_cls: type[GalaxyPopGenTypes]) -> None:
        """Initialize the enum."""
        self.model_cls = model_cls

    @classmethod
    def get_help_text(cls) -> str:
        """Return the help text for the galaxy population distribution."""
        headers = ["Enum", "Model", "Params"]
        rows = []

        for value in cls:
            rows.append([value] + value.model_cls.help_text())

        return "\b" + tabulate(rows, headers=headers, tablefmt=DEFAULT_TABLE_FMT)

    @classmethod
    def get_help_metavar(cls) -> str:
        """Return the help metavar for the galaxy population distribution."""
        metavar = "[" + "|".join([value.value for value in cls]) + "]"
        return metavar


class HaloProfileType(StrEnum):
    """Halo density profile types."""

    NFW = auto()
    EINASTO = auto()
    HERNQUIST = auto()


class IntegMethod(GEnum):
    """Cluster WL redshift-integral method (``NcDataClusterWLIntegMethod``)."""

    # pylint: disable=no-member
    LNINT = Nc.DataClusterWLIntegMethod.LNINT
    FIXED_NODES = Nc.DataClusterWLIntegMethod.FIXED_NODES
    CUBATURE = Nc.DataClusterWLIntegMethod.CUBATURE


DEFAULT_INTEG_N_NODES = 10
DEFAULT_INTEG_RULE_N = 5
DEFAULT_INTEG_NODE_RELTOL = 1.0e-4
DEFAULT_INTEG_MAX_TOTAL_NODES = 2000


@dataclass(frozen=True, kw_only=True)
class IntegMethodOptions:
    """Redshift-integral method and its FIXED_NODES-only tuning knobs.

    ``n_nodes``/``rule_n``/``auto_nodes``/``node_reltol``/``max_total_nodes``
    are plain ``NcDataClusterWLFactor`` properties that only affect
    FIXED_NODES; setting them under LNINT/CUBATURE is harmless (ignored).
    """

    integ_method: IntegMethod = IntegMethod.LNINT
    auto_nodes: bool = False
    n_nodes: int = DEFAULT_INTEG_N_NODES
    rule_n: int = DEFAULT_INTEG_RULE_N
    node_reltol: float = DEFAULT_INTEG_NODE_RELTOL
    max_total_nodes: int = DEFAULT_INTEG_MAX_TOTAL_NODES

    def apply(self, cluster_data: Nc.DataClusterWLFactor) -> None:
        """Apply these settings to a NcDataClusterWLFactor instance."""
        cluster_data.set_integ_method(self.integ_method.genum)
        cluster_data.set_auto_nodes(self.auto_nodes)
        cluster_data.set_n_nodes(self.n_nodes)
        cluster_data.set_rule_n(self.rule_n)
        cluster_data.set_node_reltol(self.node_reltol)
        cluster_data.set_max_total_nodes(self.max_total_nodes)


class HaloPositionData(BaseModel):
    """Cluster position parameters."""

    model_config = ConfigDict(extra="forbid", frozen=True)

    ra: Annotated[float, Field(ge=-180.0, le=180.0)]
    dec: Annotated[float, Field(ge=-90.0, le=90.0)]
    z: Annotated[float, Field(ge=0.0, le=5.0)] = 0.2


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


class ClusterModel(BaseModel):
    """Cluster model parameters."""

    model_config = ConfigDict(extra="forbid", arbitrary_types_allowed=True)

    mass_def: Annotated[MassDef, Field()] = MassDef.CRITICAL
    delta: Annotated[float, Field(gt=0.0)] = 200.0
    profile_type: Annotated[HaloProfileType, Field()] = HaloProfileType.NFW
    position: Annotated[HaloPositionData, Field()] = HaloPositionData(
        ra=12.34, dec=-55.123, z=0.2
    )
    cluster_c: Annotated[float, Field(ge=0.1, le=30.0)] = 4.0
    r_min: Annotated[float, Field(gt=0.0)] = 0.3
    r_max: Annotated[float, Field(gt=0.0)] = 3.0
    cluster_mass: Annotated[float, Field(ge=1.0e10, le=1.0e17)] = 1.0e14
    dist: Annotated[None | Nc.Distance, Field()] = Nc.Distance.new(5.0)

    _halo_mass_summary: Nc.HaloCMParam = PrivateAttr()
    _density_profile: Nc.HaloDensityProfile = PrivateAttr()
    _surface_mass_density: Nc.WLSurfaceMassDensity = PrivateAttr()
    _halo_position: Nc.HaloPosition = PrivateAttr()

    def model_post_init(self, _, /):
        """Initialize the cluster model."""
        self._halo_mass_summary = Nc.HaloCMParam.new(self.mass_def.genum, self.delta)
        match self.profile_type:
            case HaloProfileType.NFW:
                self._density_profile = Nc.HaloDensityProfileNFW.new(
                    self._halo_mass_summary
                )
            case HaloProfileType.EINASTO:
                self._density_profile = Nc.HaloDensityProfileEinasto.new(
                    self._halo_mass_summary
                )
            case HaloProfileType.HERNQUIST:
                self._density_profile = Nc.HaloDensityProfileHernquist.new(
                    self._halo_mass_summary
                )
            case _:  # pragma: no cover
                raise ValueError(f"Invalid halo profile type: {self.profile_type}")
        self._surface_mass_density = Nc.WLSurfaceMassDensity.new(self.dist)
        self._halo_position = Nc.HaloPosition.new(self.dist)
        self._halo_position["ra"] = self.position.ra
        self._halo_position["dec"] = self.position.dec
        self._halo_position["z"] = self.position.z
        self._halo_mass_summary["cDelta"] = self.cluster_c
        self._halo_mass_summary["log10MDelta"] = np.log10(self.cluster_mass)

    @property
    def position_data(self) -> HaloPositionData:
        """Return the cluster position."""
        return HaloPositionData(
            ra=self._halo_position["ra"],
            dec=self._halo_position["dec"],
            z=self._halo_position["z"],
        )

    @property
    def mass(self) -> float:
        """Return the cluster mass."""
        return 10.0 ** self._halo_mass_summary["log10MDelta"]

    @property
    def concentration(self) -> float:
        """Return the cluster concentration."""
        return self._halo_mass_summary["cDelta"]

    @property
    def density_profile(self) -> Nc.HaloDensityProfile:
        """Return the density profile."""
        return self._density_profile

    @property
    def surface_mass_density(self) -> Nc.WLSurfaceMassDensity:
        """Return the surface mass density."""
        return self._surface_mass_density

    @property
    def halo_position(self) -> Nc.HaloPosition:
        """Return the halo position."""
        return self._halo_position

    @property
    def halo_mass_summary(self) -> Nc.HaloCMParam:
        """Return the halo mass summary."""
        return self._halo_mass_summary

    def prepare(self, cosmo: Nc.HICosmo) -> None:
        """Prepare the cluster model."""
        self._surface_mass_density.prepare(cosmo)
        self._halo_position.prepare(cosmo)


class GalaxyDistributionModel:
    """Galaxy distribution model parameters."""

    def __init__(
        self,
        galaxies: GalaxyDistributionData,
        z_gen: GalaxyZGenTypes = GalaxyZGenComposed(),
        shape_gen: GalaxyShapeFactorGenTypes = GalaxyShapeFactorGenVarAdd(),
        pop_gen: GalaxyPopGenTypes = GalaxyPopGenGauss(),
    ) -> None:
        """Initialize the galaxy distribution model."""
        check_shape_pop_compat(shape_gen, pop_gen)

        self.position_factor = Nc.GalaxyPositionFactorFlat.new(
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
        self.pop_gen = pop_gen

    def generate_data(
        self,
        cosmo: Nc.HICosmo,
        cluster: ClusterModel,
        rng: Ncm.RNG,
        integ_options: IntegMethodOptions = IntegMethodOptions(),
    ) -> tuple[Nc.DataClusterWLFactor, Ncm.MSet]:
        """Generate the galaxy data.

        ``resample()`` retries each row's position draw until it lands in
        [r_min, r_max], so every requested row is kept -- ``self.n_galaxies``
        is the exact final count, not an upper bound.
        """
        cluster.prepare(cosmo)
        redshift_factor = self.z_gen.get_redshift_factor()
        shape_factor = self.shape_gen.get_shape_factor()

        mset = Ncm.MSet.new_array(
            [
                cosmo,
                cluster.density_profile,
                cluster.surface_mass_density,
                cluster.halo_position,
            ]
        )
        self.z_gen.register_models(mset)
        self.pop_gen.register_models(mset)
        mset.prepare_fparam_map()

        cluster.halo_position.prepare_if_needed(cosmo)

        pos_data = Nc.GalaxyPositionFactorData.new(self.position_factor, mset)
        z_data = Nc.GalaxyRedshiftFactorData.new(redshift_factor, mset)
        s_data = Nc.GalaxyShapeFactorData.new(shape_factor, mset, pos_data, z_data)
        cols = Nc.GalaxyShapeFactorData.required_columns(s_data)

        obs = Nc.GalaxyWLObs.new(
            self.shape_gen.ellip_conv.genum,
            self.shape_gen.ellip_coord.genum,
            self.n_galaxies,
            cols,
        )
        for i in range(self.n_galaxies):
            for c in cols:
                obs.set(c, i, 0.0)
            self.z_gen.write_calib(obs, i)
            self.shape_gen.write_calib(obs, i, rng)
            self.pop_gen.write_calib(obs, i, rng)

        cluster_data = Nc.DataClusterWLFactor.new(
            self.position_factor, redshift_factor, shape_factor
        )
        cluster_data.set_cut(cluster.r_min, cluster.r_max)
        integ_options.apply(cluster_data)
        cluster_data.set_obs(obs)
        cluster_data.resample(mset, rng)

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
    profile_type: HaloProfileType = HaloProfileType.NFW,
    z_gen: GalaxyZGenTypes,
    shape_gen: GalaxyShapeFactorGenTypes,
    pop_gen: GalaxyPopGenTypes,
    density: float,
    seed: None | int,
    integ_options: IntegMethodOptions,
    summary: bool,
) -> Ncm.ObjDictStr:
    """Generate a cluster weak lensing experiment.

    This function generates a cluster with galaxies following the LSST SRD redshift
    distribution.
    """
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
        profile_type=profile_type,
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
        pop_gen=pop_gen,
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

    cluster_data, mset = galaxy_model.generate_data(
        create_cosmo(), cluster, rng, integ_options=integ_options
    )

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
        table.add_row("Mass Definition", f"{cluster.mass_def}")
        table.add_row("Mass Delta", f"{cluster.delta}")
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


def _resolve_cluster_value(
    obs: Nc.GalaxyWLObs, key: str, override: float | None, human_name: str
) -> float:
    """Resolve a cluster property from an explicit override or catalog metadata.

    Raises:
        ValueError: Neither @override nor matching catalog metadata (see
            ncm_catalog_peek_meta()) is available.
    """
    if override is not None:
        return override

    meta = obs.peek_meta()

    if meta is not None:
        ok, val = meta.get_double(key)
        if ok:
            return val

    raise ValueError(
        f"{human_name} was not given and the catalog carries no '{key}' "
        "metadata; pass it explicitly."
    )


def load_cluster_wl(
    *,
    obs: Nc.GalaxyWLObs,
    cluster_ra: float | None = None,
    cluster_dec: float | None = None,
    cluster_z: float | None = None,
    cluster_mass: float,
    cluster_mass_min: float,
    cluster_mass_max: float,
    cluster_c: float | None = None,
    r_min: float,
    r_max: float,
    profile_type: HaloProfileType = HaloProfileType.NFW,
    shape_gen: GalaxyShapeFactorGenTypes,
    pop_gen: GalaxyPopGenTypes,
    integ_options: IntegMethodOptions,
    summary: bool,
) -> Ncm.ObjDictStr:
    """Build a cluster weak lensing experiment from a real NcGalaxyWLObs catalog.

    Pairs the observed galaxies -- including their per-galaxy photometric
    redshift ``NcmSpline`` -- with a ``NcGalaxyRedshiftFactorSpline`` and the
    chosen ``shape_gen`` scheme, around a cluster at (cluster_ra, cluster_dec,
    cluster_z, cluster_c). The sky footprint for the position factor is the
    bounding box of the observed galaxies' RA/Dec.

    cluster_ra/cluster_dec/cluster_z/cluster_c fall back to the catalog's own
    metadata (see ncm_catalog_peek_meta(), keys ``cluster_ra``/``cluster_dec``/
    ``cluster_z``/``cluster_c``) when not given -- the galaxy catalog itself
    only describes the observed background galaxies, not the lensing cluster.
    """
    if not (obs.has_column("ra") and obs.has_column("dec")):
        raise ValueError("obs must have 'ra' and 'dec' columns")

    obs_ellip_conv = obs.get_ellip_conv()
    obs_coord = obs.get_coord()

    if (shape_gen.ellip_conv.genum != obs_ellip_conv) or (
        shape_gen.ellip_coord.genum != obs_coord
    ):
        raise ValueError(
            "shape_dist ellip_conv/ellip_coord must match the catalog's own "
            f"convention: ellip_conv={obs_ellip_conv.value_nick}, "
            f"ellip_coord={obs_coord.value_nick}."
        )

    check_shape_pop_compat(shape_gen, pop_gen)

    cluster_ra = _resolve_cluster_value(
        obs, CATALOG_META_CLUSTER_RA, cluster_ra, "cluster_ra"
    )
    cluster_dec = _resolve_cluster_value(
        obs, CATALOG_META_CLUSTER_DEC, cluster_dec, "cluster_dec"
    )
    cluster_z = _resolve_cluster_value(
        obs, CATALOG_META_CLUSTER_Z, cluster_z, "cluster_z"
    )
    cluster_c = _resolve_cluster_value(
        obs, CATALOG_META_CLUSTER_C, cluster_c, "cluster_c"
    )

    columns = list(obs.peek_columns())
    data = obs.peek_data().to_numpy()
    # NcHaloPosition's "ra" parameter (and HaloPositionData) are defined over
    # the signed [-180, 180] convention, not the standard [0, 360) RA used by
    # real survey catalogs -- wrap both the galaxies and the cluster into it.
    ra = ((data[:, columns.index("ra")] + 180.0) % 360.0) - 180.0
    dec = data[:, columns.index("dec")]
    cluster_ra = ((cluster_ra + 180.0) % 360.0) - 180.0

    gal_ra_min = float(ra.min())
    gal_ra_max = float(ra.max())
    gal_dec_min = float(dec.min())
    gal_dec_max = float(dec.max())

    # The catalog only carries galaxy positions, not the lensing cluster's own
    # RA/Dec -- cluster_ra/cluster_dec is an external estimate. Outside the
    # galaxy window there is no lensing signal to constrain it, so the fit
    # bounds are exactly that window; reject an estimate that falls outside
    # it up front, rather than let the fit silently explore uninformative
    # space (or leave an unreachable initial point for ESMCMC). Report RA and
    # DEC together so a bad guess on both axes doesn't take two round-trips.
    ra_bad = not (gal_ra_min <= cluster_ra <= gal_ra_max)
    dec_bad = not (gal_dec_min <= cluster_dec <= gal_dec_max)

    if ra_bad or dec_bad:
        raise ValueError(
            f"Cluster position ({cluster_ra}, {cluster_dec}) outside the "
            f"catalog's RA/DEC window ([{gal_ra_min}, {gal_ra_max}], "
            f"[{gal_dec_min}, {gal_dec_max}])"
        )

    halo_position_data = HaloPositionData(ra=cluster_ra, dec=cluster_dec, z=cluster_z)
    cluster = ClusterModel(
        position=halo_position_data,
        cluster_mass=cluster_mass,
        cluster_c=cluster_c,
        r_min=r_min,
        r_max=r_max,
        profile_type=profile_type,
    )
    cosmo = create_cosmo()
    cluster.prepare(cosmo)

    cluster.halo_position.param_set_desc(
        "ra", {"lower-bound": gal_ra_min, "upper-bound": gal_ra_max}
    )
    cluster.halo_position.param_set_desc(
        "dec", {"lower-bound": gal_dec_min, "upper-bound": gal_dec_max}
    )
    cluster.halo_mass_summary.param_set_desc(
        "log10MDelta",
        {
            "lower-bound": float(np.log10(cluster_mass_min)),
            "upper-bound": float(np.log10(cluster_mass_max)),
        },
    )

    mset = Ncm.MSet.new_array(
        [
            cosmo,
            cluster.density_profile,
            cluster.surface_mass_density,
            cluster.halo_position,
        ]
    )
    pop_gen.register_models(mset)
    mset.prepare_fparam_map()

    cluster.halo_position.prepare_if_needed(cosmo)

    position_factor = Nc.GalaxyPositionFactorFlat.new(
        float(ra.min()), float(ra.max()), float(dec.min()), float(dec.max())
    )
    redshift_factor = Nc.GalaxyRedshiftFactorSpline.new()
    shape_factor = shape_gen.get_shape_factor()

    cluster_data = Nc.DataClusterWLFactor.new(
        position_factor, redshift_factor, shape_factor
    )
    cluster_data.set_cut(cluster.r_min, cluster.r_max)
    integ_options.apply(cluster_data)
    cluster_data.set_obs(obs)

    dset = Ncm.Dataset.new_array([cluster_data])
    likelihood = Ncm.Likelihood.new(dset)

    experiment = Ncm.ObjDictStr()
    experiment.set("likelihood", likelihood)
    experiment.set("model-set", mset)

    if summary:
        console = Console()

        table = Table(title="Cluster Parameters")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("RA", f"{halo_position_data.ra}")
        table.add_row("DEC", f"{halo_position_data.dec}")
        table.add_row("z", f"{halo_position_data.z}")
        table.add_row("Mass Profile", f"{cluster.profile_type}")
        table.add_row("Mass Definition", f"{cluster.mass_def}")
        table.add_row("Mass Delta", f"{cluster.delta}")
        table.add_row("Mass", f"{cluster.mass:.2e}")
        table.add_row("Concentration", f"{cluster.concentration}")
        table.add_row("Minimum Radius", f"{cluster.r_min}")
        table.add_row("Maximum Radius", f"{cluster.r_max}")
        console.print(table)

        table = Table(title="Galaxy Sample")
        table.add_column("Parameter")
        table.add_column("Value")
        table.add_row("Number of galaxies", f"{obs.len()}")
        table.add_row("RA range", f"[{ra.min():.4f}, {ra.max():.4f}]")
        table.add_row("DEC range", f"[{dec.min():.4f}, {dec.max():.4f}]")
        table.add_row("Shape Generator", repr(shape_gen))
        console.print(table)

    return experiment
