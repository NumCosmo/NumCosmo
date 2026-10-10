#
# __init__.py
#
# Wed Mar 12 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# __init__.py
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

"""Cross-correlation kernel visualization and analysis tools."""

from .cls import ComputeCls, EllSpacing, sample_ells
from .common import XcorClosureOption, XcorKernelCommon
from .kernels import (
    KERNEL_CONFIG_REGISTRY,
    KernelCMBISWConfig,
    KernelCMBLensingConfig,
    KernelNumberCountsConfig,
    KernelTSZConfig,
    KernelWeakLensingConfig,
    LSSTBinType,
    get_kernel_registry_help_text,
    parse_kernel_spec,
)
from .view import ListKernels, ViewKernel

__all__ = [
    "KERNEL_CONFIG_REGISTRY",
    "ComputeCls",
    "EllSpacing",
    "KernelCMBISWConfig",
    "KernelCMBLensingConfig",
    "KernelNumberCountsConfig",
    "KernelTSZConfig",
    "KernelWeakLensingConfig",
    "LSSTBinType",
    "ListKernels",
    "ViewKernel",
    "XcorClosureOption",
    "XcorKernelCommon",
    "get_kernel_registry_help_text",
    "parse_kernel_spec",
    "sample_ells",
]
