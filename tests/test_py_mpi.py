#!/usr/bin/env python
#
# example_mpi_fit.py
#
# Mon May 22 16:00:00 2023
# Copyright  2023  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# example_mpi_fit.py
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

"""Example testing MPI fit."""

import sys
from numpy.testing import assert_allclose
import numpy as np

from numcosmo_py import Ncm

#
#  Initializing the library objects, this must be called before
#  any other library function.
#
sys.argv = Ncm.cfg_init_full(sys.argv)


def test_mpi_fit_run_array() -> None:
    """Example testing MPI fit."""

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


def test_mpi_fit_run_array_async() -> None:
    """Example testing MPI fit."""

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


def test_mpi_job_test_run_array() -> None:
    """Example testing MPI objects."""

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
    """Example testing MPI objects."""

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


if __name__ == "__main__":
    test_mpi_fit_run_array()
    test_mpi_fit_run_array_async()
    test_mpi_job_test_run_array()
    test_mpi_job_test_run_array_async()
