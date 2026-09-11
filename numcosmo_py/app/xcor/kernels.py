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

import re
import shlex
from typing import (
    Annotated,
    Any,
    ClassVar,
    Union,
    Type,
    cast,
    get_args,
    get_origin,
)

from pydantic import BaseModel, BeforeValidator, Field, ConfigDict
from pydantic.fields import FieldInfo
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


class CMBLensingSource(GEnum):
    """Where the CMB photons are placed along the line of sight.

    :cvar THIN_SCREEN: one source plane at the decoupling redshift -
        CLI value: 'thin-screen'
    :cvar VISIBILITY: the recombination visibility function, normalized over the
        last-scattering shell - CLI value: 'visibility'
    :cvar VISIBILITY_REIONIZATION: the full visibility function, reionization
        bump included - CLI value: 'visibility-reionization'
    """

    # pylint: disable=no-member
    THIN_SCREEN = Nc.XcorKernelCMBLensingSource.THIN_SCREEN
    VISIBILITY = Nc.XcorKernelCMBLensingSource.VISIBILITY
    VISIBILITY_REIONIZATION = Nc.XcorKernelCMBLensingSource.VISIBILITY_REIONIZATION

    @classmethod
    def __get_pydantic_core_schema__(
        cls, _source_type: Any, _handler: Any
    ) -> core_schema.CoreSchema:
        """Get the Pydantic core schema for CMBLensingSource."""
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
    :ivar source: Placement of the CMB sources: a thin screen at decoupling, the
        recombination visibility function, or the full visibility including
        reionization.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelCMBLensing
    label: ClassVar[str] = "CMB Lensing"

    model_config = ConfigDict(extra="forbid", frozen=True)

    lmax: Annotated[int, Field(gt=0)] = 3000
    source: CMBLensingSource = CMBLensingSource.THIN_SCREEN

    @classmethod
    def from_args(cls, args: list[str]) -> "KernelCMBLensingConfig":
        """Create a KernelCMBLensingConfig from command line arguments.

        :param args: List of key=value strings.
        :return: Validated configuration object.
        :raises ValidationError: If arguments are invalid.
        """
        opts = parse_options_strict(args)
        return cls.model_validate(opts)


class KernelCMBISWConfig(BaseModel):
    """CMB Integrated Sachs-Wolfe (ISW) kernel configuration.

    The ISW effect is caused by the time-varying gravitational potential
    as photons traverse large-scale structures.

    :ivar lmax: Maximum multipole for noise power spectrum.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelCMBISW
    label: ClassVar[str] = "CMB ISW"

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


class KernelTSZConfig(BaseModel):
    """Thermal Sunyaev-Zeldovich (tSZ) kernel configuration.

    The tSZ effect is caused by inverse Compton scattering of CMB photons
    off hot electrons in galaxy clusters and large-scale structure.

    :ivar zmax: Maximum redshift for integration.
    """

    nc_type: ClassVar[type] = Nc.XcorKerneltSZ
    label: ClassVar[str] = "tSZ"

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

    nc_type: ClassVar[type] = Nc.XcorKernelGal
    label: ClassVar[str] = "Number Counts"

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

    nc_type: ClassVar[type] = Nc.XcorKernelWeakLensing
    label: ClassVar[str] = "Weak Lensing"

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


class KernelClusterTophatConfig(BaseModel):
    """Cluster number counts kernel configuration (thin-z approximation).

    This kernel represents cluster number counts in a single redshift bin
    using a simple top-hat window function. It implements the thin-z
    approximation where clusters are assumed to occupy a narrow redshift
    range.

    :ivar z_lower: Lower edge of the redshift bin.
    :ivar z_upper: Upper edge of the redshift bin.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelClusterTophat
    label: ClassVar[str] = "Cluster Tophat"

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

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticGauss
    label: ClassVar[str] = "Gaussian"

    chi_mean: Annotated[float, Field(gt=0.0)] = 1500.0
    chi_sigma: Annotated[float, Field(gt=0.0)] = 300.0
    n_sigma: Annotated[float, Field(gt=0.0)] = 4.0


class KernelRadialTophatConfig(_KernelRadialConfig):
    """Top-hat window with hard edges.

    :ivar chi_lower: Lower edge, in Mpc.
    :ivar chi_upper: Upper edge, in Mpc.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticTophat
    label: ClassVar[str] = "Top-hat"

    chi_lower: Annotated[float, Field(ge=0.0)] = 500.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 2500.0


class KernelRadialTophatSmoothConfig(_KernelRadialConfig):
    """Top-hat convolved with a Gaussian, i.e. what a photometric bin looks like.

    :ivar chi_lower: Lower edge of the top-hat, in Mpc.
    :ivar chi_upper: Upper edge of the top-hat, in Mpc.
    :ivar chi_sigma: Smoothing scale, in Mpc.
    :ivar n_sigma: Truncation half-width beyond the edges, in units of sigma.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticTophatSmooth
    label: ClassVar[str] = "Smoothed top-hat"

    chi_lower: Annotated[float, Field(ge=0.0)] = 1000.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 2000.0
    chi_sigma: Annotated[float, Field(gt=0.0)] = 150.0
    n_sigma: Annotated[float, Field(gt=0.0)] = 6.0


class KernelRadialStudentTConfig(_KernelRadialConfig):
    """Student-t window: non-exponential tails, as physical kernels have.

    :ivar chi_mean: Window centre, in Mpc.
    :ivar chi_scale: Window scale, in Mpc.
    :ivar nu: Degrees of freedom, which tunes the tail exponent.
    :ivar n_scale: Truncation half-width, in units of the scale.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticStudentT
    label: ClassVar[str] = "Student-t"

    chi_mean: Annotated[float, Field(gt=0.0)] = 1500.0
    chi_scale: Annotated[float, Field(gt=0.0)] = 200.0
    nu: Annotated[float, Field(gt=0.0)] = 2.0
    n_scale: Annotated[float, Field(gt=0.0)] = 6.0


class KernelRadialPowerExpConfig(_KernelRadialConfig):
    """Skewed, broad window: chi^alpha exp(-(chi/chi_scale)^beta).

    Covers the shapes an LSST-like dn/dz and an ISW kernel take.

    :ivar chi_scale: Scale, in Mpc.
    :ivar alpha: Power-law index.
    :ivar beta: Exponential index.
    :ivar chi_lower: Lower truncation, in Mpc.
    :ivar chi_upper: Upper truncation, in Mpc.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticPowerExp
    label: ClassVar[str] = "Power-exponential"

    chi_scale: Annotated[float, Field(gt=0.0)] = 1200.0
    alpha: Annotated[float, Field(gt=0.0)] = 2.0
    beta: Annotated[float, Field(gt=0.0)] = 1.5
    chi_lower: Annotated[float, Field(ge=0.0)] = 50.0
    chi_upper: Annotated[float, Field(gt=0.0)] = 4000.0


class KernelRadialLensingConfig(_KernelRadialConfig):
    """Lensing-efficiency window over a top-hat source distribution.

    Broad and smoother than its source, with a hard edge at chi_lower where it
    does not vanish -- which is why its transform decays as 1/k.

    :ivar chi_lower: Observer-side edge, in Mpc.
    :ivar chi_source_lower: Lower edge of the source distribution, in Mpc.
    :ivar chi_source_upper: Upper edge of the source distribution, in Mpc.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticLensing
    label: ClassVar[str] = "Lensing"

    chi_lower: Annotated[float, Field(ge=0.0)] = 50.0
    chi_source_lower: Annotated[float, Field(gt=0.0)] = 2000.0
    chi_source_upper: Annotated[float, Field(gt=0.0)] = 3000.0


class KernelRadialMultiConfig(_KernelRadialConfig):
    """Several Gaussian bumps, which may overlap or be disjoint.

    Bumps whose truncated supports meet form one component; disjoint ones become
    separate components, so this is the shape that exercises disconnected support.

    :ivar chi_mean: Bump centres, in Mpc.
    :ivar chi_sigma: Bump widths, in Mpc.
    :ivar weight: Relative bump amplitudes.
    :ivar n_sigma: Truncation half-width, in units of sigma.
    """

    nc_type: ClassVar[type] = Nc.XcorKernelAnalyticMulti
    label: ClassVar[str] = "Multi-bump"

    chi_mean: FloatList = [1000.0, 1600.0]
    chi_sigma: FloatList = [300.0, 300.0]
    weight: FloatList = [1.0, 0.6]
    n_sigma: Annotated[float, Field(gt=0.0)] = 4.0


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


_IVAR_RE = re.compile(r"^:ivar (\w+):\s*(.*)$")


def _field_descriptions(config_class: Type[BaseModel]) -> dict[str, str]:
    """Collect the ``:ivar:`` documentation of a config's fields.

    The class docstring is the single place a parameter is described, so the
    tables below read it rather than repeating it: a parameter documented there
    is documented on the command line too, and one that is not shows up as
    undocumented in the test that checks every field has a description.

    :param config_class: A kernel configuration class.
    :return: Mapping from field name to its documented description.
    """
    docs: dict[str, str] = {}

    for klass in reversed(config_class.__mro__):
        current: str | None = None
        for raw in (klass.__doc__ or "").splitlines():
            line = raw.strip()
            match = _IVAR_RE.match(line)
            if match is not None:
                current = match.group(1)
                docs[current] = match.group(2).strip()
            elif current is not None:
                # A blank line or another field directive ends the entry;
                # anything else is its continuation.
                if not line or line.startswith(":"):
                    current = None
                else:
                    docs[current] = f"{docs[current]} {line}"

    return docs


def _constraint_text(field: FieldInfo) -> str:
    """Describe the range a field is validated against.

    :param field: Pydantic field information.
    :return: A human-readable range, empty when the field is unconstrained.
    """
    bounds = []

    for meta in field.metadata:
        for attr, symbol in (("gt", ">"), ("ge", ">="), ("lt", "<"), ("le", "<=")):
            value = getattr(meta, attr, None)
            if value is not None:
                bounds.append(f"{symbol} {value:g}")

    return ", ".join(bounds)


def _type_text(annotation: Any) -> str:
    """Name the type a field takes on the command line.

    :param annotation: The field's type annotation.
    :return: A short type name.
    """
    if get_origin(annotation) is list:
        return f"{get_args(annotation)[0].__name__} list"

    return getattr(annotation, "__name__", str(annotation))


def _default_text(value: Any) -> str:
    """Render a default the way it would be typed back in.

    :param value: The field's default value.
    :return: The default as a command-line value.
    """
    if isinstance(value, list):
        return ",".join(f"{item:g}" for item in value)

    return str(value)


def _ordered_fields(config_class: Type[BaseModel]) -> list[str]:
    """List a config's fields, the ones it declares itself first.

    :param config_class: A kernel configuration class.
    :return: Field names, shape-specific ones before inherited ones.
    """
    fields = config_class.model_fields
    own = [
        name for name in getattr(config_class, "__annotations__", {}) if name in fields
    ]

    return own + [name for name in fields if name not in own]


def kernel_parameter_summary(config_class: Type[BaseModel]) -> str:
    """Summarise a kernel's parameters as the key=value list it is given as.

    :param config_class: A kernel configuration class.
    :return: A ``key=default`` list, comma separated.
    """
    fields = config_class.model_fields

    return ", ".join(
        f"{name}={_default_text(fields[name].default)}"
        for name in _ordered_fields(config_class)
    )


def get_kernel_registry_help_text() -> str:
    """Generate formatted help text for all available kernel types.

    Returns a formatted table showing kernel names, the NumCosmo object each
    builds, and its parameters with their defaults.

    :return: Formatted help text as a string.
    """
    headers = ["Kernel Type", "Builds", "Parameters (with defaults)"]
    rows = [
        [
            kernel_name,
            config_class.nc_type.__name__,  # type: ignore[attr-defined]
            kernel_parameter_summary(config_class),
        ]
        for kernel_name, config_class in KERNEL_CONFIG_REGISTRY.items()
    ]

    return tabulate(
        rows, headers=headers, tablefmt="rounded_grid", disable_numparse=True
    )


def get_kernel_help_text(kernel_name: str) -> str:
    """Generate formatted help text for one kernel type.

    Documents every parameter the kernel takes: its type, default, the range it
    is validated against and what it means.

    :param kernel_name: A key of :data:`KERNEL_CONFIG_REGISTRY`.
    :return: Formatted help text as a string.
    :raises ValueError: If the kernel name is not recognized.
    """
    if kernel_name not in KERNEL_CONFIG_REGISTRY:
        available = ", ".join(KERNEL_CONFIG_REGISTRY)
        raise ValueError(
            f"Unknown kernel type '{kernel_name}'. Available types: {available}"
        )

    config_class = KERNEL_CONFIG_REGISTRY[kernel_name]
    fields = config_class.model_fields
    descriptions = _field_descriptions(config_class)

    # The prose above the :ivar: block describes the shape itself.
    summary = []
    for raw in (config_class.__doc__ or "").splitlines():
        line = raw.strip()
        if line.startswith(":ivar"):
            break
        summary.append(line)

    rows = [
        [
            name,
            _type_text(fields[name].annotation),
            _default_text(fields[name].default),
            _constraint_text(fields[name]),
            descriptions.get(name, "").replace("``", ""),
        ]
        for name in _ordered_fields(config_class)
    ]
    table = tabulate(
        rows,
        headers=["Parameter", "Type", "Default", "Range", "Description"],
        tablefmt="rounded_grid",
        maxcolwidths=[None, None, None, None, 44],
        disable_numparse=True,
    )

    example = " ".join(
        f"{name}={_default_text(fields[name].default)}"
        for name in _ordered_fields(config_class)
    )

    nc_type_name = config_class.nc_type.__name__  # type: ignore[attr-defined]

    return "\n".join(
        [
            f"{kernel_name} -- builds {nc_type_name}",
            "",
            "\n".join(summary).strip(),
            "",
            table,
            "",
            f'Example: --kernel "{kernel_name} {example}"',
        ]
    )


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
