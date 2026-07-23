#!/usr/bin/env python
#
# test_planck_release.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_release.py
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

"""Tests for the curated native Planck likelihood release (build + local load)."""

import os

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import find_baseline_file
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH
from numcosmo_py.experiments.planck_native_release import (
    PlanckReleaseId,
    build_release,
    load_planck_release,
    release_filename,
)

Ncm.cfg_init()

needs_data = pytest.mark.skipif(
    find_baseline_file(COMMANDER_RELPATH) is None,
    reason="Planck baseline (plc_3.0) clik data not found",
)

_EXPECTED_TYPE = {
    PlanckReleaseId.COMMANDER: "DataPlanckCommander",
    PlanckReleaseId.SIMALL_EE: "DataPlanckSimall",
    PlanckReleaseId.SIMALL_BB: "DataPlanckSimall",
    PlanckReleaseId.SIMALL_EEBB: "DataPlanckSimall",
    PlanckReleaseId.PLIK_TT: "DataPlanckSmica",
    PlanckReleaseId.PLIK_TTTEEE: "DataPlanckSmica",
    PlanckReleaseId.PLIK_LITE_TT: "DataPlanckPlikLite",
    PlanckReleaseId.PLIK_LITE_TTTEEE: "DataPlanckPlikLite",
    PlanckReleaseId.LENSING: "DataPlanckLensing",
    PlanckReleaseId.LENSING_MARGED: "DataPlanckLensing",
}


@pytest.mark.app
@needs_data
def test_build_release_and_load(tmp_path):
    """build_release writes serialized objects that reload to the right types."""
    out = str(tmp_path)
    written = build_release(out_dir=out)
    assert written

    for rid in PlanckReleaseId:
        path = os.path.join(out, release_filename(rid))
        if not os.path.exists(path):
            continue  # source clik data for this id was missing
        data = load_planck_release(rid, cache_dir=out)
        assert data.__class__.__name__ == _EXPECTED_TYPE[rid]


@pytest.mark.app
@needs_data
def test_release_block_evaluates(tmp_path):
    """A block loaded from the release self-configures its CBE and evaluates."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    out = str(tmp_path)
    build_release(out_dir=out)

    cbe = Nc.HIPertBoltzmannCBE.new()  # bare: the loaded block configures it
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    data = load_planck_release(PlanckReleaseId.COMMANDER, cbe, cache_dir=out)
    assert isinstance(data, Nc.DataPlanckCommander)
    data.prepare(mset)
    assert np.isfinite(data.m2lnL_val(mset))
