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
from numcosmo_py.ncm import MSet, Matrix

Ncm.cfg_init()


class DataGaussTest(Ncm.DataGauss):
    """Test class for NcmDataGauss."""

    def __init__(self, corr: float = 0.0, sigma1: float = 1.0, sigma2: float = 1.0):
        """Constructor for DataGaussTest."""

        oneminuscorr = 1.0 - corr * corr
        inv_cov_array = [
            1.0 / (sigma1 * sigma1 * oneminuscorr),
            -corr / (sigma1 * sigma2 * oneminuscorr),
            -corr / (sigma1 * sigma2 * oneminuscorr),
            1.0 / (sigma2 * sigma2 * oneminuscorr),
        ]

        inv_cov = Ncm.Matrix.new_array(inv_cov_array, 2)
        mean = Ncm.Vector.new_array([0.0, 0.0])
        super().__init__(n_points=2, inv_cov=inv_cov, mean=mean, init=True)

    def do_mean_func(  # pylint: disable=arguments-differ
        self, mset: Ncm.MSet, vp: Ncm.Vector
    ) -> None:
        """Do mean function."""

        mvnd = mset.peek(Ncm.ModelMVND.id())
        assert mvnd is not None
        assert mvnd.vparam_len(0) == 2

        vp.set(0, mvnd.orig_vparam_get(0, 0))
        vp.set(1, mvnd.orig_vparam_get(0, 1))


class DataGaussTestUpdateCov(DataGaussTest):
    """Test class for NcmDataGauss with update covariance."""

    def do_inv_cov_func(  # pylint: disable=arguments-differ
        self, _: MSet, inv_cov: Matrix
    ) -> bool:
        """Do inverse covariance function."""

        inv_cov.memcpy(self.peek_inv_cov())

        return True


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
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    sv = Ncm.StatsVec.new(5, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(2)

    n_runs = 100
    data_dist = DataGaussTest(corr=0.5, sigma1=1.0, sigma2=2.0)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)

        sv.set(0, data_dist.peek_mean().get(0))
        sv.set(1, data_dist.peek_mean().get(1))
        sv.set(2, data_dist.m2lnL_val(mset))
        data_dist.leastsquares_f(mset, residual)
        sv.set(3, residual.get(0))
        sv.set(4, residual.get(1))
        sv.update()

    assert_allclose(sv.get_mean(0), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_mean(1), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_sd(0), 1.0, atol=0.1)
    assert_allclose(sv.get_sd(1), 2.0, atol=0.1)
    assert_allclose(sv.get_cor(0, 1), 0.5, atol=0.1)

    assert_allclose(sv.get_mean(2), 2.0, atol=0.4)

    # Residuals should have mean zero and standard deviation one and
    # be uncorrelated
    assert_allclose(sv.get_mean(3), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_mean(4), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_sd(3), 1.0, atol=0.1)
    assert_allclose(sv.get_sd(4), 1.0, atol=0.1)
    assert_allclose(sv.get_cor(3, 4), 0.0, atol=0.1)


def test_data_gauss_bootstrap():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    sv = Ncm.StatsVec.new(3, Ncm.StatsVecType.COV, False)

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
        sv.set(0, data_dist.peek_mean().get(0))
        sv.set(1, data_dist.peek_mean().get(1))
        sv.set(2, data_dist.m2lnL_val(mset))
        sv.update()

        assert np.sum(freq) == 2

    assert_allclose(sv.get_mean(2), 0.0, atol=0.4)

    data_dist.set_size(10)
    assert bootstrap.get_bsize() == 10
    assert bootstrap.get_fsize() == 10


def test_data_gauss_serialize():
    """Test NcmDataGauss."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])

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


def test_data_gauss_update_cov_resample():
    """Test NcmDataGauss."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    sv = Ncm.StatsVec.new(5, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(2)

    n_runs = 100
    data_dist = DataGaussTestUpdateCov(corr=0.5, sigma1=1.0, sigma2=2.0)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)

        sv.set(0, data_dist.peek_mean().get(0))
        sv.set(1, data_dist.peek_mean().get(1))
        sv.set(2, data_dist.m2lnL_val(mset))
        data_dist.leastsquares_f(mset, residual)
        sv.set(3, residual.get(0))
        sv.set(4, residual.get(1))
        sv.update()

    assert_allclose(sv.get_mean(0), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_mean(1), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_sd(0), 1.0, atol=0.1)
    assert_allclose(sv.get_sd(1), 2.0, atol=0.1)
    assert_allclose(sv.get_cor(0, 1), 0.5, atol=0.1)

    assert_allclose(sv.get_mean(2), 2.0, atol=0.4)

    # Residuals should have mean zero and standard deviation one and
    # be uncorrelated
    assert_allclose(sv.get_mean(3), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_mean(4), 0.0, atol=3.0 / math.sqrt(n_runs))
    assert_allclose(sv.get_sd(3), 1.0, atol=0.1)
    assert_allclose(sv.get_sd(4), 1.0, atol=0.1)
    assert_allclose(sv.get_cor(3, 4), 0.0, atol=0.1)


def test_data_gauss_fisher():
    """Test NcmDataGauss fisher matrix."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    data_dist = DataGaussTest(corr=0.5, sigma1=1.0, sigma2=2.0)

    data_dist.resample(mset, rng)
    fisher = data_dist.fisher_matrix(mset)

    fisher.cholesky_decomp(ord("U"))
    fisher.cholesky_inverse(ord("U"))
    fisher.copy_triangle(ord("U"))

    assert_allclose(fisher.dup_array(), [1.0, 1.0, 1.0, 4.0])


def test_data_gauss_fisher_bias():
    """Test NcmDataGauss fisher matrix and bias."""

    true_theta = [1.0, 2.0]
    theta_shift = [0.1, -0.02]

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    mset.fparam_set(0, true_theta[0] - theta_shift[0])
    mset.fparam_set(1, true_theta[1] - theta_shift[1])

    data_dist = DataGaussTest(corr=0.5, sigma1=1.0, sigma2=2.0)

    data_dist.resample(mset, rng)
    fisher0 = data_dist.fisher_matrix(mset)
    fisher1, delta_theta = data_dist.fisher_matrix_bias(
        mset,
        Ncm.Vector.new_array(true_theta),
    )

    assert_allclose(fisher0.dup_array(), fisher1.dup_array(), atol=1.0e-6)

    fisher0.cholesky_solve(delta_theta, ord("U"))
    assert_allclose(delta_theta.dup_array(), theta_shift, atol=1.0e-1)


if __name__ == "__main__":
    test_data_gauss_set_get_size()
    test_data_gauss_get_inv_cov()
    test_data_gauss_resample()
    test_data_gauss_bootstrap()
    test_data_gauss_serialize()
    test_data_gauss_update_cov_resample()
    test_data_gauss_fisher()
    test_data_gauss_fisher_bias()
