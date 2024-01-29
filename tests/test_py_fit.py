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

from numcosmo_py import Ncm

Ncm.cfg_init()


def test_fit_factory(capfd):
    """Test NcmFit factory."""

    rng = Ncm.RNG.seeded_new(None, 1234)
    model_mvnd = Ncm.ModelMVND.new(5)
    data_mvnd = Ncm.DataGaussCovMVND.new_full(5, 1.0, 2.0, 30.0, 0.0, 0.0, rng)

    mset = Ncm.MSet.new_array([model_mvnd])
    mset.param_set_all_ftype(Ncm.ParamType.FREE)
    likelihood = Ncm.Likelihood.new(Ncm.Dataset.new_array([data_mvnd]))

    fit = Ncm.Fit.factory(
        Ncm.FitType.GSL_MMS, None, likelihood, mset, Ncm.FitGradType.NUMDIFF_CENTRAL
    )

    model_mvnd.param_set(0, 1.0e300)
    print(model_mvnd.param_get(0))
    fit.run(Ncm.FitRunMsgs.NONE)

    _, err = capfd.readouterr()

    assert re.search(r"initial point provides m2lnL", err)
