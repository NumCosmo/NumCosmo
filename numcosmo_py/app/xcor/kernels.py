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

from pydantic import BaseModel, Field, ConfigDict
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
    Y1_LENS = Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_LENS
    Y1_SOURCE = Nc.GalaxySDTrueRedshiftLSSTSRDType.Y1_SOURCE
    Y10_LENS = Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_LENS
    Y10_SOURCE = Nc.GalaxySDTrueRedshiftLSSTSRDType.Y10_SOURCE

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
    """

    model_config = ConfigDict(extra="forbid", frozen=True)

    survey: Annotated[str, Field()] = "LSST-Y1"
    bin_idx: Annotated[int, Field(ge=0)] = 0
    bias: Annotated[float, Field(gt=0.0)] = 1.5
    mag_bias: Annotated[float, Field()] = 0.0
    domagbias: Annotated[bool, Field()] = True

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
            "survey=LSST-Y1, bin_idx=0, bias=1.5, mag_bias=0.0, domagbias=True",
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


# Type alias for all kernel configuration types
KernelConfigTypes = Union[
    KernelCMBLensingConfig,
    KernelCMBISWConfig,
    KernelTSZConfig,
    KernelNumberCountsConfig,
    KernelWeakLensingConfig,
]


# Registry mapping CLI names to configuration classes
KERNEL_CONFIG_REGISTRY: dict[str, Type[BaseModel]] = {
    "cmb_lensing": KernelCMBLensingConfig,
    "cmb_isw": KernelCMBISWConfig,
    "tsz": KernelTSZConfig,
    "number-counts": KernelNumberCountsConfig,
    "weak-lensing": KernelWeakLensingConfig,
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
