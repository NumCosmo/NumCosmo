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

from typing import Optional

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


class CrossValidationMethod(GEnum):
    """Cross validation methods for Ncm.StatsDist."""

    # pylint: disable=no-member
    NONE = Ncm.StatsDistCV.NONE
    SPLIT = Ncm.StatsDistCV.SPLIT
    SPLIT_NOFIT = Ncm.StatsDistCV.SPLIT_NOFIT


def create_stats_dist(
    *,
    robust: bool = False,
    interpolation_method: InterpolationMethod = InterpolationMethod.VKDE,
    interpolation_kernel: InterpolationKernel = InterpolationKernel.CAUCHY,
    cv_method: CrossValidationMethod = CrossValidationMethod.NONE,
    dim: int = 2,
    over_smooth: float = 1.0,
    split_fraction: Optional[float] = None,
    local_fraction: Optional[float] = None,
    verbose: bool = False,
):
    """Create a new interpolation object.

    :param robust: Use robust covariance matrix estimation.
    :param interpolation_method: Interpolation method.
    :param interpolation_kernel: Interpolation kernel.
    :param cv_method: Cross validation method.
    :param dim: Dimension of the interpolation.
    :param over_smooth: Oversmoothing factor.
    :param split_fraction: Split fraction.
    :param local_fraction: Local fraction.
    :param verbose: Verbose output.

    :return: A new Ncm.StatsDist object.
    """

    if interpolation_kernel == InterpolationKernel.CAUCHY:
        kernel = Ncm.StatsDistKernelST.new(dim, 1.0)
    elif interpolation_kernel == InterpolationKernel.ST3:
        kernel = Ncm.StatsDistKernelST.new(dim, 3.0)
    elif interpolation_kernel == InterpolationKernel.GAUSS:
        kernel = Ncm.StatsDistKernelGauss.new(dim)
    else:
        raise RuntimeError(f"Kernel {interpolation_kernel} not supported")

    if interpolation_method == InterpolationMethod.KDE:
        sdist = Ncm.StatsDistKDE.new(kernel, cv_method.genum)
        if local_fraction is not None:
            raise RuntimeError("local_fraction not supported for KDE")
    elif interpolation_method == InterpolationMethod.VKDE:
        sdist = Ncm.StatsDistVKDE.new(kernel, cv_method.genum)
        if local_fraction is not None:
            sdist.set_local_frac(local_fraction)

    sdist.set_over_smooth(over_smooth)
    if robust:
        sdist.set_cov_type(Ncm.StatsDistKDECovType.ROBUST)

    if split_fraction is not None:
        sdist.set_split_frac(split_fraction)

    sdist.set_print_fit(verbose)

    return sdist
