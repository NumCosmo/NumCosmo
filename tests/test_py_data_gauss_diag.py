#!/usr/bin/env python
#
# test_py_data_gauss_diag.py
#
# thu Dec 29 22:09:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_gauss_diag.py
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

"""Tests on NcmDataGaussDiag class."""

from numpy.testing import assert_allclose
import numpy as np
from numcosmo_py import Ncm
from numcosmo_py.ncm import MSet

Ncm.cfg_init()


class DataGaussDiagTest(Ncm.DataGaussDiag):
    """Test class for NcmDataGaussDiag."""

    def __init__(self, n_points: int = 200):
        """Constructor for DataGaussTest."""

        mean_array = np.linspace(-1.0, 1.0, n_points)
        sigma_array = np.linspace(1.0e-1, 2.0e-1, n_points)

        mean = Ncm.Vector.new_array(mean_array.tolist())
        sigma = Ncm.Vector.new_array(sigma_array.tolist())
        super().__init__(n_points=n_points, sigma=sigma, mean=mean, init=True)

    def do_mean_func(  # pylint: disable=arguments-differ
        self, mset: Ncm.MSet, vp: Ncm.Vector
    ) -> None:
        """Do mean function."""

        mvnd = mset.peek(Ncm.ModelMVND.id())
        assert mvnd is not None
        assert mvnd.vparam_len(0) == 2

        a = mvnd.orig_vparam_get(0, 0)
        b = mvnd.orig_vparam_get(0, 1)

        vp.set_array(np.linspace(a, b, self.get_size()).tolist())


class DataGaussDiagTestUpdateSigma(DataGaussDiagTest):
    """Test class for NcmDataGauss with update covariance."""

    def do_sigma_func(  # pylint: disable=arguments-differ
        self, _: MSet, sigma: Ncm.Vector
    ) -> bool:
        """Do inverse covariance function."""
        if self.peek_std() != sigma:
            sigma.memcpy(self.peek_std())

        return True


def test_data_gauss_set_get_size():
    """Test NcmDataGauss."""

    n_points = 200

    data_dist = DataGaussDiagTest(n_points=n_points)
    assert data_dist.get_size() == n_points
    assert data_dist.get_dof() == n_points

    data_dist.set_size(n_points * 2)
    assert data_dist.get_size() == n_points * 2
    assert data_dist.get_dof() == n_points * 2

    data_dist.set_property("w-mean", True)
    assert data_dist.get_dof() == n_points * 2 - 1


def test_data_gauss_get_inv_cov():
    """Test NcmDataGauss."""

    n_points = 200
    data_dist = DataGaussDiagTest(n_points=n_points)

    sigma = data_dist.peek_std()

    assert sigma.len() == n_points


def test_data_gauss_resample():
    """Test NcmDataGauss."""

    n_points = 20
    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    sv = Ncm.StatsVec.new(2 * n_points + 1, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(n_points)

    mset.fparam_set(0, 5.0)
    mset.fparam_set(1, 15.0)

    n_runs = 100
    data_dist = DataGaussDiagTest(n_points=n_points)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        data_dist.leastsquares_f(mset, residual)

        mean = data_dist.peek_mean().dup_array()
        ls_f = residual.dup_array()
        m2lnL = [data_dist.m2lnL_val(mset)]
        sv.append(
            Ncm.Vector.new_array(mean + ls_f + m2lnL),
            True,
        )

    mean = sv.peek_mean().dup_array()

    assert_allclose(mean[:n_points], np.linspace(5.0, 15.0, n_points), atol=0.6)
    assert_allclose(mean[n_points : 2 * n_points], 0.0, atol=0.6)
    assert_allclose(mean[2 * n_points], n_points, atol=n_points * 0.2)


def test_data_gauss_bootstrap():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    n_runs = 100
    n_points = 200
    data_dist = DataGaussDiagTest(n_points=n_points)

    data_dist.resample(mset, rng)
    data_dist.bootstrap_create()

    assert data_dist.bootstrap_enabled()
    assert isinstance(data_dist.peek_bootstrap(), Ncm.Bootstrap)

    bootstrap = data_dist.peek_bootstrap()
    for _ in range(n_runs):
        data_dist.bootstrap_resample(rng)
        index_freq = np.array(bootstrap.get_sortncomp())
        freq = index_freq[1::2]
        sv.set(0, data_dist.m2lnL_val(mset))
        sv.update()

        assert np.sum(freq) == n_points

    assert_allclose(sv.get_mean(0), n_points, atol=10.0)

    data_dist.set_size(10)
    assert bootstrap.get_bsize() == 10
    assert bootstrap.get_fsize() == 10


def test_data_gauss_serialize():
    """Test NcmDataGauss."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])

    n_points = 200
    data_dist = DataGaussDiagTest(n_points=n_points)

    data_dist.resample(mset, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    data_dist_dup = ser.dup_obj(data_dist)

    assert data_dist_dup.get_size() == data_dist.get_size()

    assert_allclose(
        data_dist.peek_std().dup_array(), data_dist_dup.peek_std().dup_array()
    )
    assert_allclose(
        data_dist.peek_mean().dup_array(), data_dist_dup.peek_mean().dup_array()
    )


def test_data_gauss_update_cov_resample():
    """Test NcmDataGauss."""

    n_points = 20
    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    sv = Ncm.StatsVec.new(2 * n_points + 1, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(n_points)

    mset.fparam_set(0, -2.0)
    mset.fparam_set(1, 2.0)

    n_runs = 100
    data_dist = DataGaussDiagTestUpdateSigma(n_points=n_points)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        data_dist.leastsquares_f(mset, residual)

        mean_vec = data_dist.peek_mean().dup_array()
        ls_f = residual.dup_array()
        m2lnL = [data_dist.m2lnL_val(mset)]
        sv.append(
            Ncm.Vector.new_array(mean_vec + ls_f + m2lnL),
            True,
        )

    mean = sv.peek_mean().dup_array()

    assert_allclose(mean[:n_points], np.linspace(-2.0, 2.0, n_points), atol=0.6)
    assert_allclose(mean[n_points : 2 * n_points], 0.0, atol=0.6)
    assert_allclose(mean[2 * n_points], n_points, atol=n_points * 0.2)


def test_data_gauss_fisher():
    """Test NcmDataGauss fisher matrix."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    mset.fparam_set(0, 1.0)
    mset.fparam_set(1, 2.0)

    n_points = 100
    data_dist = DataGaussDiagTest(n_points=n_points)

    data_dist.resample(mset, rng)
    fisher = data_dist.fisher_matrix(mset)

    fisher.cholesky_decomp(ord("U"))
    fisher.cholesky_inverse(ord("U"))
    fisher.copy_triangle(ord("U"))

    assert np.all(np.isfinite(fisher.dup_array()))


def test_data_gauss_resample_wmean():
    """Test NcmDataGauss."""

    n_points = 20
    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    sv = Ncm.StatsVec.new(2 * n_points + 1, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(n_points)

    mset.fparam_set(0, 5.0)
    mset.fparam_set(1, 15.0)

    n_runs = 100
    data_dist = DataGaussDiagTest(n_points=n_points)
    data_dist.set_property("w-mean", True)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        data_dist.leastsquares_f(mset, residual)

        mean = data_dist.peek_mean().dup_array()
        ls_f = residual.dup_array()
        m2lnL = [data_dist.m2lnL_val(mset)]
        sv.append(
            Ncm.Vector.new_array(mean + ls_f + m2lnL),
            True,
        )

    mean = sv.peek_mean().dup_array()

    assert_allclose(mean[:n_points], np.linspace(5.0, 15.0, n_points), atol=0.6)
    assert_allclose(mean[n_points : 2 * n_points], 0.0, atol=0.6)
    assert_allclose(mean[2 * n_points], n_points, atol=n_points * 0.2)


if __name__ == "__main__":
    test_data_gauss_set_get_size()
    test_data_gauss_get_inv_cov()
    test_data_gauss_resample()
    test_data_gauss_bootstrap()
    test_data_gauss_serialize()
    test_data_gauss_update_cov_resample()
    test_data_gauss_fisher()
    test_data_gauss_resample_wmean()
