#!/usr/bin/env python
#
# test_planck18_native.py
#
# Fri August 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck18_native.py
# Copyright (C) 2026 Sandro Dias Pinto Vitenti <vitenti@uel.br>
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

"""Tests for the native (clik-free) Planck 2018 experiment assembly.

``generate_planck18_native`` glues the four native likelihood blocks onto one
shared CBE. The blocks themselves are covered by the per-likelihood tests; here
they are stood in for by synthetic ones (see tests/python/fixtures_planck.py) so
the assembly, the priors and the derived-function wiring are exercised without a
local ``plc_3.0`` tree.
"""

import numpy as np
import pytest

from python.fixtures_planck import (
    make_commander_cldf,
    make_lensing_cldf,
    make_simall_cldf,
    make_smica_cldf,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import HIPrimModel
import numcosmo_py.experiments.planck_native_release as pnr
from numcosmo_py.experiments.planck18 import (
    Planck18Types,
    generate_planck18_native,
    mset_set_parameters,
)
from numcosmo_py.experiments.planck_commander import build_commander
from numcosmo_py.experiments.planck_lensing import build_lensing
from numcosmo_py.experiments.planck_simall import build_simall
from numcosmo_py.experiments.planck_smica import build_smica_tt, build_smica_ttteee

Ncm.cfg_init()


@pytest.fixture(name="synthetic_blocks", scope="module")
def fixture_synthetic_blocks(tmp_path_factory):
    """Builders for the synthetic stand-ins, keyed by release id."""
    root = tmp_path_factory.mktemp("planck_native")
    trees = {
        pnr.PlanckReleaseId.PR3_SIMALL_EE: (make_simall_cldf(root), build_simall),
        pnr.PlanckReleaseId.PR3_COMMANDER: (make_commander_cldf(root), build_commander),
        pnr.PlanckReleaseId.PR3_PLIK_TT: (make_smica_cldf(root), build_smica_tt),
        pnr.PlanckReleaseId.PR3_PLIK_TTTEEE: (
            make_smica_cldf(root, pol=True),
            build_smica_ttteee,
        ),
        pnr.PlanckReleaseId.PR3_LENSING: (make_lensing_cldf(root), build_lensing),
    }

    def load(rid, pb=None, cache_dir=None):  # signature of load_planck_release
        del cache_dir
        clik, builder = trees[rid]
        return builder(clik, pb)

    return load


@pytest.mark.parametrize(
    "data_type,model_name",
    [
        (Planck18Types.TT, "PlanckFICorTT"),
        (Planck18Types.TTTEEE, "PlanckFICorTTTEEE"),
    ],
)
@pytest.mark.parametrize("use_lensing", [False, True])
def test_generate_native_from_release(
    monkeypatch, synthetic_blocks, data_type, model_name, use_lensing
):
    """The native experiment assembles, evaluates and carries its derived functions.

    ``from_release=True`` is the path that needs no local clik tree; the release
    loader is replaced by synthetic blocks so the assembly runs anywhere.
    """
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    experiment, mfunc_array = generate_planck18_native(
        data_type=data_type,
        use_lensing_likelihood=use_lensing,
        from_release=True,
    )

    likelihood = experiment.peek("likelihood")
    mset = experiment.peek("model-set")
    assert isinstance(likelihood, Ncm.Likelihood)
    assert isinstance(mset, Ncm.MSet)
    assert mset.peek_by_name("NcPlanckFI").__class__.__name__ == model_name

    dset = likelihood.peek_dataset()
    assert dset.get_ndata() == (4 if use_lensing else 3)
    # Native blocks only: the experiment must reload without clik or the PLC.
    names = {dset.get_data(i).__class__.__name__ for i in range(dset.get_ndata())}
    assert "DataPlanckLKL" not in names

    assert experiment.peek("distance") is not None
    assert experiment.peek("ps-ml") is not None
    assert experiment.peek("ps-ml-filter") is not None
    assert mfunc_array.len() == 4

    mset_set_parameters(mset, data_type, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()
    for i in range(dset.get_ndata()):
        dset.get_data(i).prepare(mset)
    assert np.isfinite(dset.m2lnL_val(mset))


def test_generate_native_two_fluids_raises_precision(monkeypatch, synthetic_blocks):
    """The two-fluids primordial model raises the CBE primordial sampling."""
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    experiment, _ = generate_planck18_native(
        data_type=Planck18Types.TT,
        prim_model=HIPrimModel.TWO_FLUIDS,
        from_release=True,
    )

    mset = experiment.peek("model-set")
    assert isinstance(mset, Ncm.MSet)
    cosmo = mset.peek_by_name("NcHICosmo")
    assert isinstance(cosmo.peek_prim(), Nc.HIPrimTwoFluids)


def test_generate_native_invalid_data_type(monkeypatch, synthetic_blocks):
    """An unknown data combination is rejected."""
    monkeypatch.setattr(pnr, "load_planck_release", synthetic_blocks)

    with pytest.raises(ValueError, match="Invalid data type"):
        generate_planck18_native(data_type="TTTEEE_lensing", from_release=True)


def test_generate_native_without_local_data(monkeypatch):
    """Building from a local clik tree that is absent fails with a clear message."""
    import numcosmo_py.experiments.planck_lite as pl  # pylint: disable=C0415

    monkeypatch.setattr(pl, "find_baseline_file", lambda relpath: None)

    with pytest.raises(FileNotFoundError, match="Planck baseline data not found"):
        generate_planck18_native(data_type=Planck18Types.TT)
