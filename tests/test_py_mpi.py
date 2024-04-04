#!/usr/bin/env python
#
# test_py_mpi.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_mpi.py
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

"""Unit tests for NumCosmo MPI objects. """

import sys
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm

sys.argv = Ncm.cfg_init_full(sys.argv)


def test_mpi_job_fit_run_array() -> None:
    """Testing MPI fit using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobFit.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 4))

    mj.run_array(param_a, m2lnL_a)

    true_mean = data_mvnd.peek_mean().dup_array()

    for m2lnL in m2lnL_a:
        assert_allclose(np.array(m2lnL.dup_array())[1:], true_mean, rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_fit_funcs_run_array() -> None:
    """Testing MPI fit using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobFit.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 4 + 1))

    mj.run_array(param_a, m2lnL_a)

    true_mean = data_mvnd.peek_mean().dup_array()

    for m2lnL in m2lnL_a:
        assert_allclose(np.array(m2lnL.dup_array())[1:5], true_mean, rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_fit_run_array_async() -> None:
    """Testing MPI fit using run_array_async."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobFit.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 4))

    mj.run_array_async(param_a, m2lnL_a)

    true_mean = data_mvnd.peek_mean().dup_array()

    for m2lnL in m2lnL_a:
        assert_allclose(np.array(m2lnL.dup_array())[1:], true_mean, rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_fit_funcs_run_array_async() -> None:
    """Testing MPI fit using run_array_async."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobFit.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 4 + 1))

    mj.run_array_async(param_a, m2lnL_a)

    true_mean = data_mvnd.peek_mean().dup_array()

    for m2lnL in m2lnL_a:
        assert_allclose(np.array(m2lnL.dup_array())[1:5], true_mean, rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_test_run_array() -> None:
    """Testing MPI job test using run_array."""

    rng = Ncm.RNG.new(None)
    rng.set_random_seed(True)

    mj = Ncm.MPIJobTest.new()
    mj.set_rand_vector(12, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    a = []
    b = []
    for t in np.arange(12.0):
        a.append(Ncm.Vector.new_array([t]))
        b.append(Ncm.Vector.new(1))

    mj.run_array(a, b)

    for a_i, b_i in zip(a, b):
        i = int(a_i.get(0))
        assert_allclose(b_i.get(0), mj.props.vector.get(i), rtol=1.0e-15)


def test_mpi_job_test_run_array_async() -> None:
    """Testing MPI job test using run_array_async."""

    rng = Ncm.RNG.new(None)
    rng.set_random_seed(True)

    mj = Ncm.MPIJobTest.new()
    mj.set_rand_vector(12, rng)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    a = []
    b = []
    for t in np.arange(12.0):
        a.append(Ncm.Vector.new_array([t]))
        b.append(Ncm.Vector.new(1))

    mj.run_array(a, b)

    for a_i, b_i in zip(a, b):
        i = int(a_i.get(0))
        assert_allclose(b_i.get(0), mj.props.vector.get(i), rtol=1.0e-15)


def test_mpi_job_feval_run_array() -> None:
    """Testing MPI job feval using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobFEval.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 0))

    mj.run_array(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        mset.param_set_vector(param)
        assert_allclose(np.array(m2lnL.dup_array())[0], fit.m2lnL_val(), rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_feval_funcs_run_array() -> None:
    """Testing MPI job feval using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobFEval.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 1))

    mj.run_array(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        mset.param_set_vector(param)
        assert_allclose(np.array(m2lnL.dup_array())[0], fit.m2lnL_val(), rtol=1.0e-5)
        assert_allclose(np.array(m2lnL.dup_array())[1], param.get(0), rtol=1.0e-5)

    mj.free_all_slaves()


def test_mpi_job_feval_run_array_async() -> None:
    """Testing MPI job feval using run_array_async."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobFEval.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 0))

    mj.run_array_async(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        mset.param_set_vector(param)
        assert_allclose(np.array(m2lnL.dup_array())[0], fit.m2lnL_val(), rtol=1.0e-2)

    mj.free_all_slaves()


def test_mpi_job_feval_funcs_run_array_async() -> None:
    """Testing MPI job feval using run_array_async."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=1.0,
        sigma_max=2.0,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )

    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobFEval.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []
    for _ in range(24):
        param_a.append(
            Ncm.Vector.new_array([rng.gaussian_gen(-1.0, 1.0) for _ in range(4)])
        )
        m2lnL_a.append(Ncm.Vector.new(1 + 1))

    mj.run_array_async(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        mset.param_set_vector(param)
        assert_allclose(np.array(m2lnL.dup_array())[0], fit.m2lnL_val(), rtol=1.0e-5)
        assert_allclose(np.array(m2lnL.dup_array())[1], param.get(0), rtol=1.0e-5)

    mj.free_all_slaves()


def test_mpi_job_mcmc_run_array():
    """Testing MPI job mcmc using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=0.1,
        sigma_max=0.2,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )
    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobMCMC.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []

    for _ in range(128):
        param = [rng.uniform_gen(0.0, 1.0) for _ in range(4)]
        mset.fparams_set_array(param)
        m2lnL_star = fit.m2lnL_val()
        prob = rng.uniform_gen(0.0, 1.0)
        jump = rng.uniform_gen(-0.1, 1.0)
        param_a.append(Ncm.Vector.new_array(param + [m2lnL_star, np.log(prob), jump]))
        m2lnL_a.append(Ncm.Vector.new(1 + 1 + 0))

    mj.run_array(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        assert_allclose(
            m2lnL.get(0),
            1.0
            if param.get(6) < 0.0
            else 1.0
            if param.get(6) < np.exp(param.get(5))
            else 0.0,
        )
        assert_allclose(m2lnL.get(1), param.get(4))

    mj.free_all_slaves()


def test_mpi_job_mcmc_funcs_run_array():
    """Testing MPI job mcmc using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=0.1,
        sigma_max=0.2,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )
    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobMCMC.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []

    for _ in range(128):
        param = [rng.uniform_gen(0.0, 1.0) for _ in range(4)]
        mset.fparams_set_array(param)
        m2lnL_star = fit.m2lnL_val()
        prob = rng.uniform_gen(0.0, 1.0)
        jump = rng.uniform_gen(-0.1, 1.0)
        param_a.append(Ncm.Vector.new_array(param + [m2lnL_star, np.log(prob), jump]))
        m2lnL_a.append(Ncm.Vector.new(1 + 1 + 1))

    mj.run_array(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        accepted = True if param.get(6) < np.exp(param.get(5)) else False

        assert_allclose(m2lnL.get(0), 1.0 if accepted else 0.0)
        assert_allclose(m2lnL.get(1), param.get(4))
        if accepted:
            assert_allclose(m2lnL.get(2), param.get(0))

    mj.free_all_slaves()


def test_mpi_job_mcmc_run_array_async():
    """Testing MPI job mcmc using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=0.1,
        sigma_max=0.2,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )
    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    mj = Ncm.MPIJobMCMC.new(fit, None)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []

    for _ in range(128):
        param = [rng.uniform_gen(0.0, 1.0) for _ in range(4)]
        mset.fparams_set_array(param)
        m2lnL_star = fit.m2lnL_val()
        prob = rng.uniform_gen(0.0, 1.0)
        jump = rng.uniform_gen(-0.1, 1.0)
        param_a.append(Ncm.Vector.new_array(param + [m2lnL_star, np.log(prob), jump]))
        m2lnL_a.append(Ncm.Vector.new(1 + 1 + 0))

    mj.run_array_async(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        assert_allclose(
            m2lnL.get(0),
            1.0
            if param.get(6) < 0.0
            else 1.0
            if param.get(6) < np.exp(param.get(5))
            else 0.0,
        )
        assert_allclose(m2lnL.get(1), param.get(4))

    mj.free_all_slaves()


def test_mpi_job_mcmc_funcs_run_array_async():
    """Testing MPI job mcmc using run_array."""

    rng = Ncm.RNG.new()
    mset = Ncm.MSet.new_array([Ncm.ModelMVND.new(4)])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(
        dim=4,
        sigma_min=0.1,
        sigma_max=0.2,
        cor_level=50.0,
        mean_min=0.0,
        mean_max=1.0,
        rng=rng,
    )
    dset = Ncm.Dataset.new_array([data_mvnd])
    lh = Ncm.Likelihood(dataset=dset)
    fit = Ncm.Fit.factory(
        Ncm.FitType.NLOPT, "ln-neldermead", lh, mset, Ncm.FitGradType.NUMDIFF_FORWARD
    )

    funcs = Ncm.ObjArray.new()
    funcs.add(Ncm.PriorGaussParam.new_name("NcmModelMVND:mu_0", 0.0, 1.0))
    mj = Ncm.MPIJobMCMC.new(fit, funcs)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    mj.init_all_slaves(ser)

    param_a = []
    m2lnL_a = []

    for _ in range(128):
        param = [rng.uniform_gen(0.0, 1.0) for _ in range(4)]
        mset.fparams_set_array(param)
        m2lnL_star = fit.m2lnL_val()
        prob = rng.uniform_gen(0.0, 1.0)
        jump = rng.uniform_gen(-0.1, 1.0)
        param_a.append(Ncm.Vector.new_array(param + [m2lnL_star, np.log(prob), jump]))
        m2lnL_a.append(Ncm.Vector.new(1 + 1 + 1))

    mj.run_array_async(param_a, m2lnL_a)

    for param, m2lnL in zip(param_a, m2lnL_a):
        accepted = True if param.get(6) < np.exp(param.get(5)) else False

        assert_allclose(m2lnL.get(0), 1.0 if accepted else 0.0)
        assert_allclose(m2lnL.get(1), param.get(4))
        if accepted:
            assert_allclose(m2lnL.get(2), param.get(0))

    mj.free_all_slaves()


if __name__ == "__main__":
    test_mpi_job_fit_run_array()
    test_mpi_job_fit_funcs_run_array()
    test_mpi_job_fit_run_array_async()
    test_mpi_job_fit_funcs_run_array_async()
    test_mpi_job_test_run_array()
    test_mpi_job_test_run_array_async()
    test_mpi_job_feval_run_array()
    test_mpi_job_feval_funcs_run_array()
    test_mpi_job_feval_run_array_async()
    test_mpi_job_feval_funcs_run_array_async()
    test_mpi_job_mcmc_run_array()
    test_mpi_job_mcmc_funcs_run_array()
    test_mpi_job_mcmc_run_array_async()
    test_mpi_job_mcmc_funcs_run_array_async()
