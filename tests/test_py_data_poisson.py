#!/usr/bin/env python
#
# test_py_data_poisson.py
#
# thu Dec 29 22:09:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_data_poisson.py
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

"""Tests on NcmDataPoisson class."""

from numpy.testing import assert_allclose
import numpy as np
from numcosmo_py import Ncm

Ncm.cfg_init()


class DataPoissonTest(Ncm.DataPoisson):
    """Test class for NcmDataPoisson."""

    def __init__(self, n_bins: int = 200):
        """Constructor for DataGaussTest."""

        super().__init__(n_bins=n_bins, init=True)
        self.n_bins = n_bins

    def do_mean_func(  # pylint: disable=arguments-differ
        self, mset: Ncm.MSet, n: int
    ) -> float:
        """Do mean function."""

        mvnd = mset.peek(Ncm.ModelMVND.id())
        assert mvnd is not None
        assert mvnd.vparam_len(0) == 2

        a = mvnd.orig_vparam_get(0, 0)
        b = mvnd.orig_vparam_get(0, 1)
        return a + b / (self.n_bins - 1.0) * n


def test_data_poisson_set_get_size():
    """Test NcmDataPoisson."""

    n_bins = 200

    data_dist = DataPoissonTest(n_bins=n_bins)
    assert data_dist.get_size() == n_bins
    assert data_dist.get_dof() == n_bins

    data_dist.set_size(n_bins * 2)
    assert data_dist.get_size() == n_bins * 2
    assert data_dist.get_dof() == n_bins * 2


def test_data_poisson_init_from_vector():
    """Test NcmDataPoission."""

    nodes = Ncm.Vector.new_array(np.linspace(0.0, 1.0, 100 + 1))
    values = Ncm.Vector.new_array(np.linspace(100.0, 200.0, 100))

    data_dist = DataPoissonTest()

    data_dist.init_from_vector(nodes, values)

    assert data_dist.get_size() == 100
    assert data_dist.get_dof() == 100
    assert_allclose(data_dist.get_sum(), 100.0 * 150.0)

    values = data_dist.get_hist_vals()

    assert_allclose(values.dup_array(), values.dup_array())


def test_data_poisson_init_zero():
    """Test NcmDataPoisson."""

    data_dist = DataPoissonTest()

    nodes = Ncm.Vector.new_array(np.linspace(0.0, 1.0, 100 + 1))

    data_dist.init_zero(nodes)

    assert data_dist.get_size() == 100
    assert data_dist.get_dof() == 100
    assert_allclose(data_dist.get_sum(), 0.0)

    values = data_dist.get_hist_vals()

    assert_allclose(values.dup_array(), np.zeros(100))


def test_data_poisson_init_from_binning():
    """Test NcmDataPoisson."""

    data_dist = DataPoissonTest()

    nodes = Ncm.Vector.new_array(np.linspace(0.0, 1.0, 100 + 1))
    values = Ncm.Vector.new_array(np.linspace(0.0001, 0.9999, 1000))

    data_dist.init_from_binning(nodes, values)

    assert data_dist.get_size() == 100
    assert data_dist.get_dof() == 100
    assert_allclose(data_dist.get_sum(), 1000.0)

    values = data_dist.get_hist_vals()

    assert_allclose(values.dup_array(), np.ones(100) * 10.0)


def test_data_poisson_resample():
    """Test NcmDataPoisson."""

    n_bins = 20
    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    sv = Ncm.StatsVec.new(2 * n_bins + 1, Ncm.StatsVecType.COV, False)
    residual = Ncm.Vector.new(n_bins)

    mset.fparam_set(0, 5.0)
    mset.fparam_set(1, 15.0)

    n_runs = 100
    data_dist = DataPoissonTest(n_bins=n_bins)

    for _ in range(n_runs):
        data_dist.resample(mset, rng)
        data_dist.leastsquares_f(mset, residual)

        mean = data_dist.get_hist_vals().dup_array()
        ls_f = residual.dup_array()
        m2lnL = [data_dist.m2lnL_val(mset)]
        sv.append(
            Ncm.Vector.new_array(mean + ls_f + m2lnL),
            True,
        )

    mean = sv.peek_mean().dup_array()

    assert_allclose(
        mean[:n_bins], np.linspace(5.0, 15.0, n_bins, endpoint=True), atol=10.0
    )
    assert_allclose(mean[n_bins : 2 * n_bins], 1.0, atol=1.0)
    assert_allclose(mean[2 * n_bins], n_bins, atol=n_bins * 0.2)


def test_data_poisson_bootstrap():
    """Test NcmDataDist2D."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)

    mset.fparam_set(0, 5.0)
    mset.fparam_set(1, 15.0)

    n_runs = 100
    n_bins = 200
    data_dist = DataPoissonTest(n_bins=n_bins)

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

        assert np.sum(freq) == n_bins

    assert_allclose(sv.get_mean(0), n_bins, atol=50.0)

    data_dist.set_size(10)
    assert bootstrap.get_bsize() == 10
    assert bootstrap.get_fsize() == 10


def test_data_poisson_serialize():
    """Test NcmDataPoisson."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    mset.fparam_set(0, 5.0)
    mset.fparam_set(1, 15.0)

    n_bins = 200
    data_dist = DataPoissonTest(n_bins=n_bins)

    data_dist.resample(mset, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)

    data_dist_dup = ser.dup_obj(data_dist)

    assert data_dist_dup.get_size() == data_dist.get_size()

    assert_allclose(
        data_dist.get_hist_vals().dup_array(), data_dist_dup.get_hist_vals().dup_array()
    )
    assert_allclose(
        data_dist.get_hist_means(mset).dup_array(),
        data_dist_dup.get_hist_means(mset).dup_array(),
    )


def test_data_poisson_fisher():
    """Test NcmDataPoisson fisher matrix."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()

    mset.fparam_set(0, 1.0)
    mset.fparam_set(1, 2.0)

    n_bins = 100
    data_dist = DataPoissonTest(n_bins=n_bins)

    data_dist.resample(mset, rng)
    fisher = data_dist.fisher_matrix(mset)

    fisher.cholesky_decomp(ord("U"))
    fisher.cholesky_inverse(ord("U"))
    fisher.copy_triangle(ord("U"))

    assert np.all(np.isfinite(fisher.dup_array()))


if __name__ == "__main__":
    test_data_poisson_set_get_size()
    test_data_poisson_init_from_vector()
    test_data_poisson_init_zero()
    test_data_poisson_resample()
    test_data_poisson_bootstrap()
    test_data_poisson_serialize()
    test_data_poisson_fisher()
