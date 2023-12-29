#!/usr/bin/env python
#
# test_py_data_gauss.py
#
# thu Dec 29 15:57:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_gauss.py
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

"""Tests on NcmDataGauss class."""

import math
from numpy.testing import assert_allclose
import numpy as np
from numcosmo_py import Ncm

Ncm.cfg_init()


class DataGaussTest(Ncm.DataGauss):
    """Test class for NcmDataGauss."""

    def __init__(self, corr: float = 0.0, sigma1: float = 1.0, sigma2: float = 1.0):
        """Constructor for DataGaussTest."""

        oneminuscorr = 1.0 - corr * corr
        inv_cov = [
            1.0 / (sigma1 * sigma1 * oneminuscorr),
            -corr / (sigma1 * sigma2 * oneminuscorr),
            -corr / (sigma1 * sigma2 * oneminuscorr),
            1.0 / (sigma2 * sigma2 * oneminuscorr),
        ]

        inv_cov = Ncm.Matrix.new_array(inv_cov, 2)
        mean = Ncm.Vector.new_array([0.0, 0.0])
        super().__init__(n_points=2, inv_cov=inv_cov, mean=mean, init=True)

    def do_mean_func(  # pylint: disable=arguments-differ
        self, _: Ncm.MSet, vp: Ncm.Vector
    ) -> None:
        """Do mean function."""

        vp.set(0, 0.0)
        vp.set(1, 0.0)


def test_data_gauss_set_get_size():
    """Test NcmDataGauss."""

    data_dist = DataGaussTest()
    assert data_dist.get_size() == 2

    data_dist.set_size(10)
    assert data_dist.get_size() == 10


def test_data_gauss_get_inv_cov():
    """Test NcmDataGauss."""

    data_dist = DataGaussTest()

    inv_cov = data_dist.peek_inv_cov()

    assert inv_cov.nrows() == 2
    assert inv_cov.ncols() == 2

    assert_allclose(inv_cov.get(0, 0), 1.0)
    assert_allclose(inv_cov.get(0, 1), 0.0)
    assert_allclose(inv_cov.get(1, 0), 0.0)
    assert_allclose(inv_cov.get(1, 1), 1.0)


def test_data_gauss_resample():
    """Test NcmDataGauss."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()
    sv = Ncm.StatsVec.new(2, Ncm.StatsVecType.COV, False)

    n_runs = 100
    data_dist = DataGaussTest(corr=0.5, sigma1=1.0, sigma2=2.0)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        sv.append(data_dist.peek_mean(), True)

    assert_allclose(sv.get_mean(0), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_mean(1), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_sd(0), 1.0, atol=0.1)
    assert_allclose(sv.get_sd(1), 2.0, atol=0.1)
    assert_allclose(sv.get_cor(0, 1), 0.5, atol=0.1)


def test_data_gauss_bootstrap():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()

    n_runs = 100
    data_dist = DataGaussTest()

    data_dist.resample(mset, rng)
    data_dist.bootstrap_create()

    assert data_dist.bootstrap_enabled()
    assert isinstance(data_dist.peek_bootstrap(), Ncm.Bootstrap)

    bootstrap = data_dist.peek_bootstrap()
    for _ in range(n_runs):
        data_dist.bootstrap_resample(rng)
        index_freq = np.array(bootstrap.get_sortncomp())
        freq = index_freq[1::2]

        assert np.sum(freq) == 2


def test_data_gauss_serialize():
    """Test NcmDataGauss."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.empty_new()

    data_dist = DataGaussTest(corr=0.55, sigma1=1.234, sigma2=2.345)

    data_dist.resample(mset, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    data_dist_dup = ser.dup_obj(data_dist)

    assert data_dist_dup.get_size() == data_dist.get_size()

    assert_allclose(
        data_dist.peek_inv_cov().dup_array(), data_dist_dup.peek_inv_cov().dup_array()
    )
    assert_allclose(
        data_dist.peek_mean().dup_array(), data_dist_dup.peek_mean().dup_array()
    )


if __name__ == "__main__":
    test_data_gauss_set_get_size()
    test_data_gauss_get_inv_cov()
    test_data_gauss_resample()
    test_data_gauss_bootstrap()
    test_data_gauss_serialize()
