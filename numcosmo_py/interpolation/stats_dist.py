#
# stats_dist.py
#
# Wed Feb 8 10:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# stats_dist.py
# Copyright (C) 2023 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Create a new ensemble sampler object."""

from numcosmo_py import Ncm, GEnum


class InterpolationMethod(GEnum):
    """Possible interpolation methods Ncm.StatsDist."""

    # pylint: disable=no-member
    KDE = Ncm.FitESMCMCWalkerAPESMethod.KDE
    VKDE = Ncm.FitESMCMCWalkerAPESMethod.VKDE


class InterpolationKernel(GEnum):
    """Possible interpolation kernels for Ncm.StatsDist."""

    # pylint: disable=no-member
    CAUCHY = Ncm.FitESMCMCWalkerAPESKType.CAUCHY
    ST3 = Ncm.FitESMCMCWalkerAPESKType.ST3
    GAUSS = Ncm.FitESMCMCWalkerAPESKType.GAUSS


def create_stats_dist(
    *,
    robust: bool = False,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    dim: int = 2,
    over_smooth: float = 1.0,
):
    """Create a new interpolation object."""

    if interpolation_kernel == InterpolationKernel.CAUCHY:
        kernel = Ncm.StatsDistKernelST.new(dim, 1.0)
    elif interpolation_kernel == InterpolationKernel.ST3:
        kernel = Ncm.StatsDistKernelST.new(dim, 3.0)
    elif interpolation_kernel == InterpolationKernel.GAUSS:
        kernel = Ncm.StatsDistKernelGauss.new(dim)
    else:
        raise RuntimeError(f"Kernel {interpolation_kernel} not supported")

    if interpolation_method == InterpolationMethod.KDE:
        sdist = Ncm.StatsDistKDE.new(kernel, Ncm.StatsDistCV.NONE)
    elif interpolation_method == InterpolationMethod.VKDE:
        sdist = Ncm.StatsDistVKDE.new(kernel, Ncm.StatsDistCV.NONE)

    sdist.set_over_smooth(over_smooth)
    if robust:
        sdist.set_cov_type(Ncm.StatsDistKDECovType.ROBUST)

    return sdist
