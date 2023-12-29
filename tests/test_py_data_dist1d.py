#!/usr/bin/env python
#
# test_py_data_dist1d.py
#
# thu Dec 29 12:23:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_dist1d.py
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

"""Tests on NcmDataDist1D class."""

import math
from numpy.testing import assert_allclose
from scipy.stats import norm

from numcosmo_py import Ncm

Ncm.cfg_init()


class DataDist1dTest(Ncm.DataDist1d):
    """Test class for NcmDataDist1D."""

    def __init__(self, n_points: int):
        """Constructor for DataDist1dTest."""
        super().__init__(n_points=n_points)

    def do_dist1d_m2lnL_val(  # pylint: disable=arguments-differ
        self, _: Ncm.MSet, x: float
    ) -> float:
        """Test function for NcmDataDist1D.do_m2lnL_val."""
        return -2.0 * norm.logpdf(x)

    def do_inv_pdf(  # pylint: disable=arguments-differ
        self, _: Ncm.MSet, u: float
    ) -> float:
        """Test function for NcmDataDist1D.do_inv_pdf."""
        return norm.ppf(u)


def test_data_dist1d_set_get_size():
    """Test NcmDataDist1D."""

    n_points = 10
    data_dist = DataDist1dTest(n_points=n_points)
    assert data_dist.get_size() == n_points

    data_dist.set_size(2 * n_points)
    assert data_dist.get_size() == 2 * n_points


def test_data_dist1d_get_vector():
    """Test NcmDataDist1D."""

    n_points = 10
    data_dist = DataDist1dTest(n_points=n_points)

    vector = data_dist.get_data()

    assert vector.len() == n_points


def test_data_dist1d_resample():
    """Test NcmDataDist1D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()
    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    n_points = 100
    n_runs = 100
    data_dist = DataDist1dTest(n_points=n_points)
    lnnorm = -n_points * math.log(2.0 * math.pi)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        sv.set(0, data_dist.m2lnL_val(mset) + lnnorm)
        sv.update()

    assert_allclose(
        sv.get_mean(0), n_points, atol=3.0 * math.sqrt(2.0 * n_points / n_runs)
    )


def test_data_dist1d_bootstrap():
    """Test NcmDataDist1D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()
    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    n_points = 100
    n_runs = 100
    data_dist = DataDist1dTest(n_points=n_points)

    data_dist.resample(mset, rng)
    data_dist.bootstrap_create()

    assert data_dist.bootstrap_enabled()
    assert isinstance(data_dist.peek_bootstrap(), Ncm.Bootstrap)

    lnnorm = -n_points * math.log(2.0 * math.pi)

    for _ in range(n_runs):
        data_dist.bootstrap_resample(rng)
        sv.set(0, data_dist.m2lnL_val(mset) + lnnorm)
        sv.update()

    assert_allclose(
        sv.get_mean(0), n_points, atol=20.0 * math.sqrt(2.0 * n_points / n_runs)
    )


def test_data_dist1d_serialize():
    """Test NcmDataDist1D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()

    n_points = 100
    data_dist = DataDist1dTest(n_points=n_points)

    data_dist.resample(mset, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    data_dist_dup = ser.dup_obj(data_dist)

    assert data_dist_dup.get_size() == data_dist.get_size()
    assert data_dist_dup.get_data().len() == data_dist.get_data().len()

    vector = data_dist.get_data()
    vector_dup = data_dist_dup.get_data()

    assert_allclose(vector.dup_array(), vector_dup.dup_array())


if __name__ == "__main__":
    test_data_dist1d_set_get_size()
    test_data_dist1d_get_vector()
    test_data_dist1d_resample()
    test_data_dist1d_bootstrap()
    test_data_dist1d_serialize()
