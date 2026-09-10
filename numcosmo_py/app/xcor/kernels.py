#
# kernels.py
#
# Wed Mar 12 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# kernels.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Configuration classes for cross-correlation kernels.

This module provides Pydantic-based configuration classes for various
cross-correlation kernel types used in cosmological analyses.

Each kernel type has a dedicated configuration class that:
- Validates CLI parameters
- Provides default values
- Implements parsing from command-line arguments
- Offers human-readable help text

The module also provides utilities for parsing kernel specifications from
command-line strings.
"""

import shlex
from typing import Annotated, Any, Union, Type, cast

from pydantic import BaseModel, BeforeValidator, Field, ConfigDict
from pydantic_core import core_schema
from tabulate import tabulate

from numcosmo_py import Nc, parse_options_strict, GEnum


class LSSTBinType(GEnum):
    """LSST survey year and redshift bin type.

    :cvar Y1_LENS: Year 1 lens sample (5 bins) - CLI value: 'y1-lens'
    :cvar Y1_SOURCE: Year 1 source sample (5 bins) - CLI value: 'y1-source'
    :cvar Y10_LENS: Year 10 lens sample - CLI value: 'y10-lens'
    :cvar Y10_SOURCE: Year 10 source sample - CLI value: 'y10-source'
    """

    # pylint: disable=no-member
    Y1_LENS = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_LENS
    Y1_SOURCE = Nc.GalaxyRedshiftPopLSSTSRDType.Y1_SOURCE
    Y10_LENS = Nc.GalaxyRedshiftPopLSSTSRDType.Y10_LENS
    Y10_SOURCE = Nc.GalaxyRedshiftPopLSSTSRDType.Y10_SOURCE

    @classmethod
    def __get_pydantic_core_schema__(
        cls, _source_type: Any, _handler: Any
    ) -> core_schema.CoreSchema:
        """Get the Pydantic core schema for LSSTBinType."""
        return core_schema.no_info_before_validator_function(
            lambda v: cls(v) if isinstance(v, str) else v,
            core_schema.enum_schema(cls, list(cls), sub_type="str"),
            serialization=core_schema.plain_serializer_function_ser_schema(
                lambda v: str(v.value)
            ),
        )


class KernelCMBLensingConfig(BaseModel):
    """CMB lensing kernel configuration.

    This kernel represents the CMB lensing convergence field, which traces
    the integrated matter distribution along the line of sight.

    :ivar lmax: Maximum multipole for noise power spectrum.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    lmax: Annotated[int, Field(gt=0)] = 3000

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelCMBLensingConfig":
        """Create a KernelCMBLensingConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for CMB lensing kernel.

        :return: List containing [model name, parameter description].
        """
        return ["KernelCMBLensing", "lmax=3000"]


class KernelCMBISWConfig(BaseModel):
    """CMB Integrated Sachs-Wolfe (ISW) kernel configuration.

    The ISW effect is caused by the time-varying gravitational potential
    as photons traverse large-scale structures.

    :ivar lmax: Maximum multipole for noise power spectrum.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    lmax: Annotated[int, Field(gt=0)] = 3000

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelCMBISWConfig":
        """Create a KernelCMBISWConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for CMB ISW kernel.

        :return: List containing [model name, parameter description].
        """
        return ["KernelCMBISW", "lmax=3000"]


class KernelTSZConfig(BaseModel):
    """Thermal Sunyaev-Zeldovich (tSZ) kernel configuration.

    The tSZ effect is caused by inverse Compton scattering of CMB photons
    off hot electrons in galaxy clusters and large-scale structure.

    :ivar zmax: Maximum redshift for integration.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    zmax: Annotated[float, Field(gt=0.0)] = 6.0

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelTSZConfig":
        """Create a KernelTSZConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for tSZ kernel.

        :return: List containing [model name, parameter description].
        """
        return ["KernelTSZ", "zmax=6.0"]


class KernelNumberCountsConfig(BaseModel):
    """Galaxy number counts kernel configuration.

    This kernel represents the galaxy number density field for a specific
    redshift bin from various surveys.
    Galaxy clustering always uses lens bins.

    :ivar survey: Survey specification (e.g., 'LSST-Y1', 'LSST-Y10').
    :ivar bin_idx: Bin index within the selected survey.
    :ivar bias: Galaxy bias parameter.
    :ivar mag_bias: Magnification bias parameter.
    :ivar domagbias: Whether to include magnification bias.
    :ivar dorsd: Whether to include linear redshift-space distortions.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    survey: Annotated[str, Field()] = "LSST-Y1"
    bin_idx: Annotated[int, Field(ge=0)] = 0
    bias: Annotated[float, Field(gt=0.0)] = 1.5
    mag_bias: Annotated[float, Field()] = 0.0
    domagbias: Annotated[bool, Field()] = True
    dorsd: Annotated[bool, Field()] = False

    @property
    def bin_type(self) -> LSSTBinType:
        """Get the appropriate lens bin type based on survey specification."""
        match self.survey.upper():
            case "LSST-Y1":
                return LSSTBinType.Y1_LENS
            case "LSST-Y10":
                return LSSTBinType.Y10_LENS
            case _:
                raise ValueError(
                    f"Unknown survey '{self.survey}'. " f"Supported: LSST-Y1, LSST-Y10"
                )

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelNumberCountsConfig":
        """Create a KernelNumberCountsConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for number counts kernel.

        :return: List containing [model name, parameter description].
        """
        return [
            "KernelNumberCounts",
            "survey=LSST-Y1, bin_idx=0, bias=1.5, mag_bias=0.0, domagbias=True, "
            "dorsd=False",
        ]


class KernelWeakLensingConfig(BaseModel):
    """Weak lensing kernel configuration.

    This kernel represents the weak gravitational lensing shear field for
    a specific source redshift bin from various surveys.
    Weak lensing always uses source bins.

    :ivar survey: Survey specification (e.g., 'LSST-Y1', 'LSST-Y10').
    :ivar bin_idx: Bin index within the selected survey.
    :ivar nbar: Galaxy number density per square arcminute.
    :ivar intr_shear: Intrinsic shear dispersion.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    survey: Annotated[str, Field()] = "LSST-Y1"
    bin_idx: Annotated[int, Field(ge=0)] = 0
    nbar: Annotated[float, Field(gt=0.0)] = 3.0
    intr_shear: Annotated[float, Field(gt=0.0)] = 7.0

    @property
    def bin_type(self) -> LSSTBinType:
        """Get the appropriate source bin type based on survey specification."""
        match self.survey.upper():
            case "LSST-Y1":
                return LSSTBinType.Y1_SOURCE
            case "LSST-Y10":
                return LSSTBinType.Y10_SOURCE
            case _:
                raise ValueError(
                    f"Unknown survey '{self.survey}'. " f"Supported: LSST-Y1, LSST-Y10"
                )

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelWeakLensingConfig":
        """Create a KernelWeakLensingConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for weak lensing kernel.

        :return: List containing [model name, parameter description].
        """
        return [
            "KernelWeakLensing",
            "survey=LSST-Y1, bin_idx=0, nbar=3.0, intr_shear=7.0",
        ]


class KernelClusterTophatConfig(BaseModel):
    """Cluster number counts kernel configuration (thin-z approximation).

    This kernel represents cluster number counts in a single redshift bin
    using a simple top-hat window function. It implements the thin-z
    approximation where clusters are assumed to occupy a narrow redshift
    range.

    :ivar z_lower: Lower edge of the redshift bin.
    :ivar z_upper: Upper edge of the redshift bin.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    z_lower: Annotated[float, Field(ge=0.0)] = 0.2
    z_upper: Annotated[float, Field(gt=0.0)] = 0.8

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelClusterTophatConfig":
        """Create a KernelClusterTophatConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for cluster tophat kernel.

        :return: List containing [model name, parameter description].
        """
        return [
            "KernelClusterTophat",
            "z_lower=0.2, z_upper=0.8",
        ]


def _split_floats(value: object) -> object:
    """Accept ``a,b,c`` from the command line as a list of floats.

    The kernel options arrive as strings, one key=value at a time, so a
    vector-valued parameter has nowhere to be a list yet.

    :param value: Raw option value, a comma-separated string or an actual list.
    :return: A list of floats when given a string, otherwise the value unchanged.
    """
    if isinstance(value, str):
        return [float(part) for part in value.split(",")]

    return value


FloatList = Annotated[list[float], BeforeValidator(_split_floats)]


class _KernelRadialConfig(BaseModel):
    """Shared behaviour of the analytic radial windows.

    These are the closed-form shapes of ``NcXcorKernelRadial``, the ones certified
    against Arb (``data/truth_tables/xcor/``). Unlike every kernel above they carry
    no cosmology: the window is a function of comoving distance in Mpc, given
    directly. That is what makes them exactly known, and what makes them the right
    thing to plot a certified reference on top of.

    :ivar bessel_deriv: Derivative order of the spherical Bessel weight, 0, 1 or 2.
        Order 2 is the weight a redshift-space distortion term carries. A kernel
        with it set must be viewed non-Limber (``--l-limber -1``), since the
        redshift-space Limber methods do not implement a derivative component.
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    bessel_deriv: Annotated[int, Field(ge=0, le=2)] = 0

    @classmethod
    def from_args(cls, args: list[str]):
        """Create a configuration from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        return cls.model_validate(parse_options_strict(args))


class KernelRadialGaussConfig(_KernelRadialConfig):
    """Gaussian window in comoving distance, truncated at n_sigma.

    :ivar chi_mean: Window centre, in Mpc.
    :ivar chi_sigma: Window standard deviation, in Mpc.
    :ivar n_sigma: Truncation half-width, in units of sigma.
    """

    chi_mean: Annotated[float, Field(gt=0.0)] = 1500.0
    chi_sigma: Annotated[float, Field(gt=0.0)] = 300.0
    n_sigma: Annotated[float, Field(gt=0.0)] = 4.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the Gaussian radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticGauss",
            "chi_mean=1500, chi_sigma=300, n_sigma=4, bessel_deriv=0",
        ]


class KernelRadialTophatConfig(_KernelRadialConfig):
    """Top-hat window with hard edges.

    :ivar chi_lower: Lower edge, in Mpc.
    :ivar chi_upper: Upper edge, in Mpc.
    """

    chi_lower: Annotated[float, Field(ge=0.0)] = 500.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 2500.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the top-hat radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticTophat",
            "chi_lower=500, chi_upper=2500, bessel_deriv=0",
        ]


class KernelRadialTophatSmoothConfig(_KernelRadialConfig):
    """Top-hat convolved with a Gaussian, i.e. what a photometric bin looks like.

    :ivar chi_lower: Lower edge of the top-hat, in Mpc.
    :ivar chi_upper: Upper edge of the top-hat, in Mpc.
    :ivar chi_sigma: Smoothing scale, in Mpc.
    :ivar n_sigma: Truncation half-width beyond the edges, in units of sigma.
    """

    chi_lower: Annotated[float, Field(ge=0.0)] = 1000.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 2000.0
    chi_sigma: Annotated[float, Field(gt=0.0)] = 150.0
    n_sigma: Annotated[float, Field(gt=0.0)] = 6.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the smoothed top-hat radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticTophatSmooth",
            "chi_lower=1000, chi_upper=2000, chi_sigma=150, n_sigma=6",
        ]


class KernelRadialStudentTConfig(_KernelRadialConfig):
    """Student-t window: non-exponential tails, as physical kernels have.

    :ivar chi_mean: Window centre, in Mpc.
    :ivar chi_scale: Window scale, in Mpc.
    :ivar nu: Degrees of freedom, which tunes the tail exponent.
    :ivar n_scale: Truncation half-width, in units of the scale.
    """

    chi_mean: Annotated[float, Field(gt=0.0)] = 1500.0
    chi_scale: Annotated[float, Field(gt=0.0)] = 200.0
    nu: Annotated[float, Field(gt=0.0)] = 2.0
    n_scale: Annotated[float, Field(gt=0.0)] = 6.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the Student-t radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticStudentT",
            "chi_mean=1500, chi_scale=200, nu=2, n_scale=6",
        ]


class KernelRadialPowerExpConfig(_KernelRadialConfig):
    """Skewed, broad window: chi^alpha exp(-(chi/chi_scale)^beta).

    Covers the shapes an LSST-like dn/dz and an ISW kernel take.

    :ivar chi_scale: Scale, in Mpc.
    :ivar alpha: Power-law index.
    :ivar beta: Exponential index.
    :ivar chi_lower: Lower truncation, in Mpc.
    :ivar chi_upper: Upper truncation, in Mpc.
    """

    chi_scale: Annotated[float, Field(gt=0.0)] = 1200.0
    alpha: Annotated[float, Field(gt=0.0)] = 2.0
    beta: Annotated[float, Field(gt=0.0)] = 1.5
    chi_lower: Annotated[float, Field(ge=0.0)] = 50.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 4000.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the power-exponential radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticPowerExp",
            "chi_scale=1200, alpha=2, beta=1.5, chi_lower=50, chi_upper=4000",
        ]


class KernelRadialLensingConfig(_KernelRadialConfig):
    """Lensing-efficiency window over a top-hat source distribution.

    Broad and smoother than its source, with a hard edge at chi_lower where it
    does not vanish -- which is why its transform decays as 1/k.

    :ivar chi_lower: Observer-side edge, in Mpc.
    :ivar chi_source_lower: Lower edge of the source distribution, in Mpc.
    :ivar chi_source_upper: Upper edge of the source distribution, in Mpc.
    """

    chi_lower: Annotated[float, Field(ge=0.0)] = 50.0
    chi_source_lower: Annotated[float, Field(gt=0.0)] = 2000.0
    chi_source_upper: Annotated[float, Field(gt=0.0)] = 3000.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the lensing radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticLensing",
            "chi_lower=50, chi_source_lower=2000, chi_source_upper=3000",
        ]


class KernelRadialMultiConfig(_KernelRadialConfig):
    """Several Gaussian bumps, which may overlap or be disjoint.

    Bumps whose truncated supports meet form one component; disjoint ones become
    separate components, so this is the shape that exercises disconnected support.

    :ivar chi_mean: Bump centres, in Mpc.
    :ivar chi_sigma: Bump widths, in Mpc.
    :ivar weight: Relative bump amplitudes.
    :ivar n_sigma: Truncation half-width, in units of sigma.
    """

    chi_mean: FloatList = [1000.0, 1600.0]
    chi_sigma: FloatList = [300.0, 300.0]
    weight: FloatList = [1.0, 0.6]
    n_sigma: Annotated[float, Field(gt=0.0)] = 4.0

    @staticmethod
    def help_text() -> list[str]:
        """Return help text for the multi-bump radial window.

        :return: List containing [model name, parameter description].
        """
        return [
            "XcorKernelAnalyticMulti",
            "chi_mean=1000,1600 chi_sigma=300,300 weight=1,0.6 n_sigma=4",
        ]


# Type alias for all kernel configuration types
KernelConfigTypes = Union[
    KernelCMBLensingConfig,
    KernelCMBISWConfig,
    KernelTSZConfig,
    KernelNumberCountsConfig,
    KernelWeakLensingConfig,
    KernelClusterTophatConfig,
    KernelRadialGaussConfig,
    KernelRadialTophatConfig,
    KernelRadialTophatSmoothConfig,
    KernelRadialStudentTConfig,
    KernelRadialPowerExpConfig,
    KernelRadialLensingConfig,
    KernelRadialMultiConfig,
]


# Registry mapping CLI names to configuration classes
KERNEL_CONFIG_REGISTRY: dict[str, Type[BaseModel]] = {
    "cmb_lensing": KernelCMBLensingConfig,
    "cmb_isw": KernelCMBISWConfig,
    "tsz": KernelTSZConfig,
    "number-counts": KernelNumberCountsConfig,
    "weak-lensing": KernelWeakLensingConfig,
    "cluster-tophat": KernelClusterTophatConfig,
    "radial-gauss": KernelRadialGaussConfig,
    "radial-tophat": KernelRadialTophatConfig,
    "radial-tophat-smooth": KernelRadialTophatSmoothConfig,
    "radial-student-t": KernelRadialStudentTConfig,
    "radial-power-exp": KernelRadialPowerExpConfig,
    "radial-lensing": KernelRadialLensingConfig,
    "radial-multi": KernelRadialMultiConfig,
}


def get_kernel_registry_help_text() -> str:
    """Generate formatted help text for all available kernel types.

    Returns a formatted table showing kernel names, model names, and parameters.

    :return: Formatted help text as a string.
    """
    headers = ["Kernel Type", "Model", "Parameters"]
    rows = []

    for kernel_name, config_class in KERNEL_CONFIG_REGISTRY.items():
        help_info = config_class.help_text()  # type: ignore[attr-defined]
        rows.append([kernel_name] + help_info)

    return tabulate(rows, headers=headers, tablefmt="rounded_grid")


def parse_kernel_spec(spec: str) -> tuple[str, KernelConfigTypes]:
    """Parse a kernel specification string.

    The specification has the format:

        "<kernel_name> key=value key=value ..."

    The first token is the kernel name, which must be a key in
    KERNEL_CONFIG_REGISTRY. Remaining tokens are key=value pairs
    parsed and validated by the corresponding configuration class.

    :param spec: Kernel specification string.
    :return: A tuple of (kernel_name, config_object).
    :raises ValueError: If the kernel name is not recognized or if the
        specification is malformed.
    :raises ValidationError: If the configuration parameters are invalid.

    Examples::

        >>> parse_kernel_spec("cmb_lensing lmax=3000")
        ('cmb_lensing', KernelCMBLensingConfig(lmax=3000))

        >>> parse_kernel_spec("number-counts survey=LSST-Y1 bin_idx=0 bias=1.5")
        ('number-counts', KernelNumberCountsConfig(...))
    """
    # Split the specification using shell-like syntax
    tokens = shlex.split(spec)

    if not tokens:
        raise ValueError("Empty kernel specification")

    # First token is the kernel name
    kernel_name = tokens[0]

    # Look up the configuration class
    if kernel_name not in KERNEL_CONFIG_REGISTRY:
        available = ", ".join(KERNEL_CONFIG_REGISTRY.keys())
        raise ValueError(
            f"Unknown kernel type '{kernel_name}'. " f"Available types: {available}"
        )

    config_class = KERNEL_CONFIG_REGISTRY[kernel_name]

    # Parse remaining tokens as configuration options
    config_args = tokens[1:]

    # Create configuration object
    # All classes in the registry implement from_args classmethod
    config = config_class.from_args(config_args)  # type: ignore[attr-defined]

    return kernel_name, cast(KernelConfigTypes, config)
