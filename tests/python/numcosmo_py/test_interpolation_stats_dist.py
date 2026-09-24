#!/usr/bin/env python
#
# test_interpolation_stats_dist.py
#
# Wed Sep 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_interpolation_stats_dist.py
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

"""Unit tests for the Ncm.StatsDist builder in numcosmo_py.interpolation.

create_stats_dist() is where the command line and the sampling helpers turn a
kernel name, an estimator name and a handful of settings into an Ncm.StatsDist.
Each name must reach the right object and each setting the right setter; the
combinations the estimators refuse must be refused here, before a run starts.
"""

import numpy as np
import pytest

from numcosmo_py import Ncm
from numcosmo_py.interpolation.stats_dist import (
    CrossValidationMethod,
    InterpolationKernel,
    InterpolationMethod,
    create_stats_dist,
)

Ncm.cfg_init()


@pytest.mark.parametrize("method", list(InterpolationMethod))
@pytest.mark.parametrize("kernel", list(InterpolationKernel))
def test_create_stats_dist_kernel_and_method(method, kernel):
    """Every kernel and estimator pair builds, takes its settings and prepares."""
    dim = 2
    local_fraction = 0.5 if method == InterpolationMethod.VKDE else None
    sdist = create_stats_dist(
        robust=True,
        interpolation_method=method,
        interpolation_kernel=kernel,
        cv_method=CrossValidationMethod.NONE,
        dim=dim,
        over_smooth=1.2,
        split_fraction=0.6,
        local_fraction=local_fraction,
    )

    if method == InterpolationMethod.KDE:
        assert type(sdist) is Ncm.StatsDistKDE  # pylint: disable=unidiomatic-typecheck
    else:
        assert isinstance(sdist, Ncm.StatsDistVKDE)
        assert sdist.get_local_frac() == local_fraction

    sdist_kernel = sdist.get_kernel()
    if kernel == InterpolationKernel.GAUSS:
        assert isinstance(sdist_kernel, Ncm.StatsDistKernelGauss)
    else:
        assert isinstance(sdist_kernel, Ncm.StatsDistKernelST)
        nu = {
            InterpolationKernel.CAUCHY: 1.0,
            InterpolationKernel.ST3: 3.0,
            InterpolationKernel.AUTO: 10.0,
        }[kernel]
        assert sdist_kernel.get_nu() == nu

    assert sdist.get_auto_kernel() == (kernel == InterpolationKernel.AUTO)
    assert sdist.get_over_smooth() == 1.2
    assert sdist.get_split_frac() == 0.6
    assert sdist.get_cov_type() == Ncm.StatsDistKDECovType.ROBUST
    assert not sdist.get_center_shrink()
    assert not sdist.get_print_fit()

    rng = np.random.default_rng(20260923)
    for point in rng.normal(size=(80, dim)):
        sdist.add_obs(Ncm.Vector.new_array(point.tolist()))

    sdist.prepare(None)
    assert np.isfinite(sdist.eval(Ncm.Vector.new_array([0.0] * dim)))


def test_create_stats_dist_kde_refuses_local_fraction():
    """The local fraction belongs to the VKDE; the KDE has no use for it."""
    with pytest.raises(RuntimeError, match="local_fraction not supported"):
        create_stats_dist(
            interpolation_method=InterpolationMethod.KDE,
            interpolation_kernel=InterpolationKernel.GAUSS,
            local_fraction=0.5,
        )


def test_create_stats_dist_auto_kernel_needs_student_t():
    """The deprecated auto_kernel flag tunes a Student-t kernel; a Gaussian is refused."""
    with pytest.raises(ValueError, match="Student-t"):
        create_stats_dist(
            interpolation_kernel=InterpolationKernel.GAUSS,
            auto_kernel=True,
        )


def test_create_stats_dist_unknown_kernel():
    """A kernel name outside the enumeration is refused, not guessed."""
    with pytest.raises(RuntimeError, match="not supported"):
        create_stats_dist(
            interpolation_kernel="not-a-kernel",  # type: ignore[arg-type]
        )
