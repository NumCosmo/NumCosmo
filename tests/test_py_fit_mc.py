#!/usr/bin/env python
#
# test_py_fit_mc.py
#
# Sat Sep 13 12:53:40 2025
# Copyright  2025  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_fit_mc.py
# Copyright (C) 2025 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for NcmFit object."""

import re
import pytest
import numpy as np
import numpy.testing as npt

from numcosmo_py import Ncm

Ncm.cfg_init()

MEAN = 1.2543


@pytest.fixture(
    name="fit", params=[1, 2, 5, 10], ids=["dim=1", "dim=2", "dim=5", "dim=10"]
)
def fixture_fit(request) -> Ncm.Fit:
    """Fixture for NcmFit object."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    problem_size: int = request.param
    model_mvnd = Ncm.ModelMVND.new(problem_size)
    data_mvnd = [
        Ncm.DataGaussCovMVND.new_full(
            problem_size, 1.0e-1, 2.0e-1, 30.0, MEAN, MEAN, rng
        )
        for _ in range(4)
    ]

    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    mset.fparams_set_array(np.ones(problem_size) * MEAN)
    likelihood = Ncm.Likelihood.new(Ncm.Dataset.new_array(data_mvnd))

    return Ncm.Fit.factory(
        Ncm.FitType.GSL_MMS, None, likelihood, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )


@pytest.fixture(
    name="mc",
    params=list(Ncm.FitMCResampleType),
    ids=[x.name for x in list(Ncm.FitMCResampleType)],
)
def fixture_mc(fit: Ncm.Fit, request) -> Ncm.FitMC:
    """Fixture for NcmFitMC object."""
    return Ncm.FitMC.new(fit, request.param, Ncm.FitRunMsgs.SIMPLE)


@pytest.mark.parametrize("nthreads", [1, 4], ids=["threads=1", "threads=4"])
def test_fit_mc_run(capfd, mc: Ncm.FitMC, nthreads: int):
    """Test NcmFit factory."""
    n_runs = 100

    mc.set_nthreads(nthreads)

    mc.start_run()
    mc.run(n_runs)
    mc.end_run()

    out, _ = capfd.readouterr()

    assert re.search(r"Task:NcmFitMC, completed: 100 of 100", out)

    mcat = mc.peek_catalog()

    assert mcat.nchains() == 1
    assert mcat.len() == n_runs

    npt.assert_allclose(
        mcat.get_mean().dup_array(), np.ones(mcat.ncols() - 1) * MEAN, rtol=1.0e-1
    )


def test_serialize_deserialize(mc: Ncm.FitMC):
    """Test NcmFit serialize and deserialize."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)

    mc2: Ncm.FitMC = ser.dup_obj(mc)

    assert mc.props.nthreads == mc2.props.nthreads
    assert mc.props.rtype == mc2.props.rtype
    assert mc.props.mtype == mc2.props.mtype
