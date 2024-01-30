#!/usr/bin/env python
#
# test_py_fit.py
#
# Min Jan 29 20:26:43 2024
# Copyright  2024  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_py_fit.py
# Copyright (C) 2024 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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
import sys
import pytest

from numcosmo_py import Ncm

Ncm.cfg_init()


@pytest.fixture(name="fit")
def fixture_fit() -> Ncm.Fit:
    """Fixture for NcmFit object."""

    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(5)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(5, 1.0, 2.0, 30.0, 0.0, 0.0, rng)

    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    likelihood = Ncm.Likelihood.new(Ncm.Dataset.new_array([data_mvnd]))

    return Ncm.Fit.factory(
        Ncm.FitType.GSL_MMS, None, likelihood, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )


def test_fit_factory(capfd, fit: Ncm.Fit):
    """Test NcmFit factory."""

    mset = fit.peek_mset()
    mset.param_set(Ncm.ModelMVND.id(), 0, 1.0e300)

    fit.run(Ncm.FitRunMsgs.NONE)

    _, err = capfd.readouterr()

    assert re.search(r"initial point provides m2lnL", err)


def test_fit_log_state_non_computed(capfd, fit: Ncm.Fit):
    """Test NcmFit log_state."""

    fit.set_messages(Ncm.FitRunMsgs.SIMPLE)
    fit.log_state()

    out, _ = capfd.readouterr()
    log_lines = out.splitlines()

    assert re.fullmatch(r"#\s*m2lnL\s*=\s*0\s*", log_lines[5])


def test_fit_log_state_computed(capfd, fit: Ncm.Fit):
    """Test NcmFit log_state."""

    fit.set_messages(Ncm.FitRunMsgs.SIMPLE)
    fit.m2lnL_val()
    fit.log_state()

    out, _ = capfd.readouterr()
    log_lines = out.splitlines()

    assert re.fullmatch(r"#\s*m2lnL\s*=\s*0\s+\(.*\)\s*", log_lines[5])


def test_fit_obs_fisher_log(capfd, fit: Ncm.Fit):
    """Test NcmFit obs_fisher_log."""

    fit.set_messages(Ncm.FitRunMsgs.SIMPLE)
    fit.obs_fisher()

    out, _ = capfd.readouterr()

    assert re.search(r"Computing Hessian matrix using numerical differentiation", out)
    assert re.search(r"Computing Hessian matrix", out)
    assert re.search(r"trying", out)


def test_fit_obs_set_logger(capsys, fit: Ncm.Fit):
    """Test NcmFit obs_set_logger."""

    fit.set_messages(Ncm.FitRunMsgs.SIMPLE)

    def writer(_fit: Ncm.Fit, msg: str) -> None:
        sys.stdout.write(msg)

    def updater(_fit: Ncm.Fit, _n: int) -> None:
        sys.stdout.write(".")

    def start_update(_fit: Ncm.Fit, _msg: str) -> None:
        sys.stdout.write("#")

    def end_update(_fit: Ncm.Fit, _msg: str) -> None:
        sys.stdout.write("\n")

    fit.set_logger(writer, updater, start_update, end_update)

    fit.run(Ncm.FitRunMsgs.SIMPLE)

    out, _ = capsys.readouterr()

    assert re.search(r"GSL Multidimensional Minimization", out)
    assert re.search(r"\.{3,}", out)
    assert re.search(r"Fit parameters", out)
