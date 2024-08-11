#!/usr/bin/env python
#
# test_py_dataset.py
#
# thu Dec 31 11:29:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_dataset.py
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

"""Tests on NcmDataset class."""

import math
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm

from .test_py_data_gauss_diag import DataGaussDiagTest
from .test_py_data_gauss import DataGaussTest

Ncm.cfg_init()


def test_dataset_constructor():
    """Test constructor."""
    dset = Ncm.Dataset.new()
    assert dset is not None
    assert isinstance(dset, Ncm.Dataset)

    dset2 = dset.ref()
    assert dset2 == dset

    dset.append_data(Ncm.DataRosenbrock.new())
    dset.append_data(Ncm.DataRosenbrock.new())
    dset.append_data(Ncm.DataFunnel.new())

    assert dset.get_n() > 0
    assert dset.get_dof() > 0
    assert dset.get_length() > 0

    assert dset.all_init()


def test_dataset_new_array_eval():
    """Test eval."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    assert dset.get_ndata() == 3
    rosenbrock_model = Ncm.ModelRosenbrock.new()
    funnel_model = Ncm.ModelFunnel.new(5)

    mset = Ncm.MSet.new_array([rosenbrock_model, funnel_model])

    assert math.isfinite(dset.m2lnL_val(mset))


def test_dataset_new_array_dup():
    """Test eval."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    assert dset.get_ndata() == 3
    rosenbrock_model = Ncm.ModelRosenbrock.new()
    funnel_model = Ncm.ModelFunnel.new(5)

    mset = Ncm.MSet.new_array([rosenbrock_model, funnel_model])

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dset2 = dset.dup(ser)

    assert math.isfinite(dset.m2lnL_val(mset))
    assert math.isfinite(dset2.m2lnL_val(mset))

    assert_allclose(dset.m2lnL_val(mset), dset2.m2lnL_val(mset))


def test_dataset_new_array_copy():
    """Test eval."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    assert dset.get_ndata() == 3
    rosenbrock_model = Ncm.ModelRosenbrock.new()
    funnel_model = Ncm.ModelFunnel.new(5)

    mset = Ncm.MSet.new_array([rosenbrock_model, funnel_model])

    dset2 = dset.copy()

    assert math.isfinite(dset.m2lnL_val(mset))
    assert math.isfinite(dset2.m2lnL_val(mset))

    assert_allclose(dset.m2lnL_val(mset), dset2.m2lnL_val(mset))
    assert dset.get_ndata() == dset2.get_ndata()

    for i in range(dset.get_ndata()):
        assert dset.peek_data(i) == dset2.peek_data(i)
        assert dset.get_data(i) == dset2.get_data(i)


def test_dataset_new_array_get_data_array():
    """Test eval."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    assert dset.get_ndata() == 3

    data_array = dset.get_data_array()

    assert data_array.len() == 3
    for i in range(data_array.len()):
        assert data_array.get(i) == dset.get_data(i)


def test_dataset_resample():
    """Test resample."""
    rng = Ncm.RNG.new()
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])

    n_runs = 1000

    assert dset.has_leastsquares_f()
    assert dset.has_m2lnL_val()

    ls_f = Ncm.Vector.new(dset.get_n())

    sv = Ncm.StatsVec.new(1 + dset.get_n(), Ncm.StatsVecType.COV, False)
    for _ in range(n_runs):
        dset.resample(mset, rng)
        sv.set(0, dset.m2lnL_val(mset))

        assert_allclose(
            np.sum([dset.m2lnL_i_val(mset, i) for i in range(dset.get_length())]),
            dset.m2lnL_val(mset),
        )

        sv.update()
        dset.leastsquares_f(mset, ls_f)
        for i in range(dset.get_n()):
            sv.set(1 + i, ls_f.get(i))

    mean = np.array(sv.peek_mean().dup_array())
    sd = np.array([sv.get_sd(i) for i in range(sv.len())])
    assert_allclose(mean[0], dset.get_n(), atol=0.1)
    assert_allclose(mean[1:], 0, atol=0.15)
    assert_allclose(sd[1:], 1.0, atol=0.1)


def test_dataset_bootstrap_partial_resample():
    """Test resample."""
    rng = Ncm.RNG.new()
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    dset.resample(mset, rng)
    dset.bootstrap_set(Ncm.DatasetBStrapType.PARTIAL)
    dset.append_data(DataGaussTest(corr=0.3, sigma1=10.0, sigma2=2.0))

    n_runs = 1000

    assert dset.has_leastsquares_f()
    assert dset.has_m2lnL_val()

    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)
    for _ in range(n_runs):
        dset.bootstrap_resample(rng)
        sv.set(0, dset.m2lnL_val(mset))
        sv.update()

    assert_allclose(sv.get_mean(0), dset.get_n(), atol=20.0)


def test_dataset_bootstrap_total_resample():
    """Test resample."""
    rng = Ncm.RNG.new()
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    dset.resample(mset, rng)
    dset.bootstrap_set(Ncm.DatasetBStrapType.TOTAL)
    dset.append_data(DataGaussTest(corr=0.3, sigma1=10.0, sigma2=2.0))

    n_runs = 1000

    assert dset.has_leastsquares_f()
    assert dset.has_m2lnL_val()

    sv = Ncm.StatsVec.new(1, Ncm.StatsVecType.VAR, False)
    for _ in range(n_runs):
        dset.bootstrap_resample(rng)
        sv.set(0, dset.m2lnL_val(mset))
        sv.update()

    assert_allclose(sv.get_mean(0), dset.get_n(), atol=10.0)


def test_dataset_log_info():
    """Test log_info."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    dset.log_info()
    assert dset.get_info() is not None


def test_dataset_no_mean_vector():
    """Test mean vector."""
    dset = Ncm.Dataset.new_array(
        [Ncm.DataRosenbrock.new(), Ncm.DataRosenbrock.new(), Ncm.DataFunnel.new()]
    )

    assert not dset.has_mean_vector()


def test_dataset_mixed_mean_vector():
    """Test mean vector."""
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
            Ncm.DataFunnel.new(),
        ]
    )

    assert not dset.has_mean_vector()


def test_dataset_has_mean_vector():
    """Test mean vector."""
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )

    assert dset.has_mean_vector()


def test_dataset_mean_vector():
    """Test mean vector."""
    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)

    mset.fparams_set_array(np.array([1.0, 2.0]))

    assert dset.has_mean_vector()

    mean_vector = Ncm.Vector.new(dset.get_n())
    dset.mean_vector(mset, mean_vector)

    local_mean_vector = []
    for i in range(dset.get_ndata()):
        data = dset.get_data(i)
        mean_vector_i = Ncm.Vector.new(data.get_length())
        data.mean_vector(mset, mean_vector_i)
        local_mean_vector += mean_vector_i.dup_array()

    assert_allclose(mean_vector.dup_array(), local_mean_vector)


def test_dataset_fisher():
    """Test resample."""
    theta_true = np.array([1.0, 2.0])
    delta_theta = np.array([0.1, 0.2])

    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.fparams_set_array(theta_true - delta_theta)
    mset.prepare_fparam_map()

    fisher = dset.fisher_matrix(mset)
    fisher.cholesky_decomp(ord("U"))
    fisher.cholesky_inverse(ord("U"))
    fisher.copy_triangle(ord("U"))

    assert np.all(np.isfinite(fisher.dup_array()))


def test_dataset_fisher_bias():
    """Test resample."""
    true_theta = np.array([1.0, 2.0])
    theta_shift = np.array([0.1, 0.2])

    dset = Ncm.Dataset.new_array(
        [
            DataGaussDiagTest(n_points=200),
            DataGaussTest(corr=0.3, sigma1=1.0, sigma2=2.0),
        ]
    )
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(2)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    mset.fparams_set_array(true_theta - theta_shift)

    fisher0 = dset.fisher_matrix(mset)
    fisher1, delta_theta = dset.fisher_matrix_bias(
        mset,
        Ncm.Vector.new_array(
            np.concatenate((np.linspace(true_theta[0], true_theta[1], 200), true_theta))
        ),
    )

    assert_allclose(fisher0.dup_array(), fisher1.dup_array(), atol=1.0e-6)

    fisher0.cholesky_solve(delta_theta, ord("U"))
    assert_allclose(delta_theta.dup_array(), theta_shift, atol=1.0e-1)


if __name__ == "__main__":
    test_dataset_constructor()
    test_dataset_new_array_eval()
    test_dataset_new_array_dup()
    test_dataset_new_array_copy()
    test_dataset_resample()
    test_dataset_bootstrap_partial_resample()
    test_dataset_bootstrap_total_resample()
    test_dataset_log_info()
    test_dataset_mean_vector()
    test_dataset_fisher()
    test_dataset_fisher_bias()
