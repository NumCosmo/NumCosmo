#!/usr/bin/env python
#
# test_planck_golden.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_golden.py
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

"""Golden snapshot pinning the native Planck likelihood m2lnL values.

The per-likelihood clik-match tests are *relative* (native vs the live clik
wrapper at a runtime cosmology): they prove the native code reproduces clik, but
a change that drifts both sides together (e.g. reverting the lensing cors fix on
both the embedded clik loader and the converter) would pass silently.

This test pins the *absolute* native m2lnL of every ported Planck likelihood at a
fixed fiducial cosmology against a stored reference vector under
``data/truth_tables/planck_m2lnl_golden.bin`` (a serialized #NcmVector). The
comparison is tolerance-based so it survives sub-ULP cross-stack (libm/GSL/BLAS/
CLASS) drift while flagging real regressions, which shift m2lnL by O(0.1) or more
(the lensing bug shifted the CMB-marginalized value from 8.95 to 15.16).

Regenerate the reference (only when an intentional change moves the values)::

    python tests/python/nc/data/test_planck_golden.py
"""

import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import (
    find_baseline_file,
    PLIK_LITE_TTTEEE_RELPATH,
    build_plik_lite,
)
from numcosmo_py.experiments.planck_smica import (
    PLIK_TT_RELPATH,
    PLIK_TTTEEE_RELPATH,
    build_smica_tt,
    build_smica_ttteee,
)
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH, build_commander
from numcosmo_py.experiments.planck_simall import SIMALL_EEBB_RELPATH, build_simall
from numcosmo_py.experiments.planck_lensing import (
    LENSING_FULL_RELPATH,
    LENSING_MARGED_RELPATH,
    build_lensing,
)

Ncm.cfg_init()

GOLDEN_FILE = "truth_tables/planck_m2lnl_golden.bin"
# CLASS/libm/BLAS rounding drifts the absolute m2lnL by a few ULP across builds;
# this tolerance absorbs that, while real code regressions are O(0.1) or larger.
GOLDEN_RTOL = 1.0e-6
GOLDEN_ATOL = 1.0e-6

# (name, relpath, builder). Fixed order defines the golden vector layout.
_CASES = [
    ("plik_lite_ttteee", PLIK_LITE_TTTEEE_RELPATH, build_plik_lite),
    ("smica_tt", PLIK_TT_RELPATH, build_smica_tt),
    ("smica_ttteee", PLIK_TTTEEE_RELPATH, build_smica_ttteee),
    ("commander", COMMANDER_RELPATH, build_commander),
    ("simall_eebb", SIMALL_EEBB_RELPATH, build_simall),
    ("lensing_full", LENSING_FULL_RELPATH, build_lensing),
    ("lensing_marged", LENSING_MARGED_RELPATH, build_lensing),
]
_KEYS = [c[0] for c in _CASES]


def _make(name, relpath, builder):
    """Build one native likelihood + its mset at the fixed fiducial cosmology."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cbe = Nc.HIPertBoltzmannCBE.new()
    # Unlike NcDataPlanckLensing, the plik_lite/smica/commander/simall classes do
    # not raise the Boltzmann targets/lmax themselves (in their own tests the clik
    # wrapper does it); configure the CBE explicitly for a native-only setup.
    lmax = 2508
    cbe.set_target_Cls(
        Nc.DataCMBDataType.TT
        | Nc.DataCMBDataType.EE
        | Nc.DataCMBDataType.TE
        | Nc.DataCMBDataType.BB
        | Nc.DataCMBDataType.PHIPHI
    )
    cbe.set_TT_lmax(lmax)
    cbe.set_EE_lmax(lmax)
    cbe.set_TE_lmax(lmax)
    cbe.set_BB_lmax(lmax)
    cbe.set_PHIPHI_lmax(lmax)

    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)

    if name == "smica_ttteee":
        planck = Nc.PlanckFICorTTTEEE()
        ptype = Planck18Types.TTTEEE
    else:
        planck = Nc.PlanckFICorTT()
        ptype = Planck18Types.TT

    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, ptype, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    native = builder(find_baseline_file(relpath), cbe)
    native.prepare(mset)
    return native, mset


def _compute_all():
    """Return the m2lnL of every case, in the fixed _KEYS order."""
    values = []
    for name, relpath, builder in _CASES:
        native, mset = _make(name, relpath, builder)
        values.append(native.m2lnL_val(mset))
    return values


@pytest.mark.app
def test_planck_m2lnl_golden():
    """Native Planck m2lnL values match the stored fixed-cosmology reference."""
    for name, relpath, _ in _CASES:
        if find_baseline_file(relpath) is None:
            pytest.skip(f"baseline data not found ({relpath})")

    path = Ncm.cfg_get_data_filename(GOLDEN_FILE, True)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    golden = ser.from_binfile(path)
    assert isinstance(golden, Ncm.Vector)
    assert golden.len() == len(_KEYS)

    values = _compute_all()
    for i, name in enumerate(_KEYS):
        assert values[i] == pytest.approx(
            golden.get(i), rel=GOLDEN_RTOL, abs=GOLDEN_ATOL
        ), f"{name}: {values[i]} vs golden {golden.get(i)}"


if __name__ == "__main__":
    # Regenerate data/truth_tables/planck_m2lnl_golden.bin.
    import os

    vals = _compute_all()
    vec = Ncm.Vector.new_array(vals)
    out = os.path.join(
        os.path.dirname(__file__), "..", "..", "..", "..", "data", GOLDEN_FILE
    )
    out = os.path.abspath(out)
    ser_out = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    ser_out.to_binfile(vec, out)
    for k, v in zip(_KEYS, vals):
        print(f"  {k:20s} {v:.15g}")
    print(f"wrote {out}")
