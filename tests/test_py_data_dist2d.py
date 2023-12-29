#!/usr/bin/env python
#
# test_py_data_dist2d.py
#
# thu Dec 29 15:03:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_dist2d.py
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

"""Tests on NcmDataDist2D class."""

import math
from typing import Tuple

from numpy.testing import assert_allclose
from scipy.stats import norm

from numcosmo_py import Ncm

Ncm.cfg_init()


class DataDist2dTest(Ncm.DataDist2d):
    """Test class for NcmDataDist2D."""

    def __init__(self, n_points: int):
        """Constructor for DataDist2dTest."""
        super().__init__(n_points=n_points)

    def do_dist2d_m2lnL_val(  # pylint: disable=arguments-differ
        self,
        _: Ncm.MSet,
        x: float,
        y: float,
    ) -> float:
        """Test function for NcmDataDist2D.do_m2lnL_val."""
        return -2.0 * norm.logpdf(x) + -2.0 * norm.logpdf(y)

    def do_inv_pdf(  # pylint: disable=arguments-differ
        self,
        _: Ncm.MSet,
        u: float,
        v: float,
    ) -> Tuple[float, float]:
        """Test function for NcmDataDist2D.do_inv_pdf."""
        return (norm.ppf(u), norm.ppf(v))


def test_data_dist2d_set_get_size():
    """Test NcmDataDist2D."""

    n_points = 10
    data_dist = DataDist2dTest(n_points=n_points)
    assert data_dist.get_size() == n_points

    data_dist.set_size(2 * n_points)
    assert data_dist.get_size() == 2 * n_points


def test_data_dist2d_get_matrix():
    """Test NcmDataDist2D."""

    n_points = 10
    data_dist = DataDist2dTest(n_points=n_points)

    matrix = data_dist.get_data()

    assert matrix.nrows() == n_points
    assert matrix.ncols() == 2


def test_data_dist2d_resample():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()
    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    n_points = 100
    n_runs = 100
    data_dist = DataDist2dTest(n_points=n_points)
    lnnorm = -n_points * math.log(2.0 * math.pi)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        sv.set(0, data_dist.m2lnL_val(mset) + 2.0 * lnnorm)
        sv.update()

    assert_allclose(
        sv.get_mean(0), 2.0 * n_points, atol=3.0 * math.sqrt(4.0 * n_points / n_runs)
    )

    assert_allclose(data_dist.inv_pdf(mset, 0.5, 0.5), [0.0, 0.0])


def test_data_dist2d_bootstrap():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()
    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    n_points = 100
    n_runs = 100
    data_dist = DataDist2dTest(n_points=n_points)

    data_dist.resample(mset, rng)
    data_dist.bootstrap_create()

    assert data_dist.bootstrap_enabled()
    assert isinstance(data_dist.peek_bootstrap(), Ncm.Bootstrap)

    lnnorm = -n_points * math.log(2.0 * math.pi)

    for _ in range(n_runs):
        data_dist.bootstrap_resample(rng)
        sv.set(0, data_dist.m2lnL_val(mset) + 2.0 * lnnorm)
        sv.update()

    assert_allclose(
        sv.get_mean(0), 2.0 * n_points, atol=20.0 * math.sqrt(4.0 * n_points / n_runs)
    )


def test_data_dist2d_serialize():
    """Test NcmDataDist1D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()

    n_points = 100
    data_dist = DataDist2dTest(n_points=n_points)

    data_dist.resample(mset, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    data_dist_dup = ser.dup_obj(data_dist)

    assert data_dist_dup.get_size() == data_dist.get_size()
    assert data_dist_dup.get_data().nrows() == data_dist.get_data().nrows()

    matrix = data_dist.get_data()
    matrix_dup = data_dist_dup.get_data()

    assert_allclose(matrix.dup_array(), matrix_dup.dup_array())


if __name__ == "__main__":
    test_data_dist2d_set_get_size()
    test_data_dist2d_get_matrix()
    test_data_dist2d_resample()
    test_data_dist2d_bootstrap()
