#!/usr/bin/env python
#
# test_planck_synthetic_golden.py
#
# Fri August 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_synthetic_golden.py
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

"""Golden snapshot pinning the native Planck likelihood assembly, data-free.

``test_planck_golden.py`` pins the same quantity against the real ``plc_3.0``
data, which makes it a stronger physical statement and a much weaker regression
test: it needs data no CI lane has, so it only ever runs on a developer machine,
and with a Boltzmann solve in the loop its value tracks the Boltzmann code's
grids and the compiler's floating-point choices as much as the code under test.

This snapshot removes both limitations. The likelihoods are built from the
synthetic cldf trees (deterministic, seeded) and driven by #FixedClBoltzmann,
which hands out stored spectra instead of solving anything, so ``-2\\ln L`` is a
pure function of the committed inputs. It therefore runs in every lane on every
platform, and its tolerance can be tight enough to catch a real change in the
assembly rather than merely a large one.

What it does *not* check is the coupling to the Boltzmann code or agreement with
clik; that is ``test_planck_golden.py``'s job.

Regenerate the reference (only when an intentional change moves the values)::

    python tests/python/nc/data/test_planck_synthetic_golden.py
"""

import os
import sys

import pytest

# Run as a script (to regenerate the reference) this file's own directory is what
# lands on sys.path, not the tests root the shared fixtures are imported from.
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "..", ".."))

# flake8: noqa: E402
# pylint: disable=wrong-import-position
from python.fixtures_planck import (
    FixedClBoltzmann,
    make_commander_cldf,
    make_lensing_cldf,
    make_plik_lite_cldf,
    make_simall_cldf,
    make_smica_cldf,
    PLIK_LITE_LMAX,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import HIPrimModel, create_cosmo
from numcosmo_py.experiments.planck_commander import build_commander
from numcosmo_py.experiments.planck_lensing import build_lensing
from numcosmo_py.experiments.planck_lite import build_plik_lite
from numcosmo_py.experiments.planck_simall import build_simall
from numcosmo_py.experiments.planck_smica import build_smica_tt, build_smica_ttteee
from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

Ncm.cfg_init()

GOLDEN_FILE = "truth_tables/planck_synthetic_m2lnl_golden.bin"
# With no Boltzmann solve in the loop the only spread left is how a compiler
# arranges the same arithmetic: all seven values come out bit-identical between
# the -O2 and the -O0 (coverage) builds. The bound is kept a few orders above
# that measurement to allow for other compilers, and is still seven orders
# tighter than the real-data snapshot's -- tight enough to fail on a genuine
# change in the assembly rather than only on a large one.
GOLDEN_RTOL = 1.0e-9


def _commander(root, pb):
    return build_commander(make_commander_cldf(root), pb)


def _simall(root, pb):
    return build_simall(make_simall_cldf(root, spectra=("EE", "BB")), pb)


def _plik_lite(root, pb):
    clik = make_plik_lite_cldf(root, spectra=("TT", "TE", "EE"))
    return build_plik_lite(clik, pb, lmax=PLIK_LITE_LMAX)


def _smica_tt(root, pb):
    return build_smica_tt(make_smica_cldf(root), pb)


def _smica_ttteee(root, pb):
    return build_smica_ttteee(make_smica_cldf(root, pol=True), pb)


def _lensing(root, pb):
    return build_lensing(make_lensing_cldf(root), pb)


def _lensing_marged(root, pb):
    return build_lensing(make_lensing_cldf(root, marged=True), pb)


# (name, builder, needs the TTTEEE nuisance model). Fixed order defines the
# golden vector layout; append, never reorder.
_CASES = [
    ("commander", _commander, False),
    ("simall_eebb", _simall, False),
    ("plik_lite_ttteee", _plik_lite, False),
    ("smica_tt", _smica_tt, False),
    ("smica_ttteee", _smica_ttteee, True),
    ("lensing_full", _lensing, False),
    ("lensing_marged", _lensing_marged, False),
]
_KEYS = [c[0] for c in _CASES]


def _compute_all(root):
    """Return the m2lnL of every case, in the fixed _KEYS order."""
    values = []
    for name, builder, pol in _CASES:
        pb = FixedClBoltzmann()
        planck = Nc.PlanckFICorTTTEEE() if pol else Nc.PlanckFICorTT()
        planck.params_set_default_ftype()
        # The cosmology only supplies the parameter slots here -- the spectra
        # come from the fixed Boltzmann, so nothing is solved from it.
        cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
        mset = Ncm.MSet.new_array([planck, cosmo])
        mset_set_parameters(
            mset,
            Planck18Types.TTTEEE if pol else Planck18Types.TT,
            HIPrimModel.POWER_LAW,
        )
        mset.prepare_fparam_map()

        data = builder(root, pb)
        data.prepare(mset)
        values.append(data.m2lnL_val(mset))

    return values


def test_planck_synthetic_m2lnl_golden(tmp_path):
    """The synthetic Planck m2lnL values match the stored reference."""
    path = Ncm.cfg_get_data_filename(GOLDEN_FILE, True)
    golden = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).from_binfile(path)
    assert isinstance(golden, Ncm.Vector)
    assert golden.len() == len(_KEYS)

    values = _compute_all(tmp_path)
    for i, name in enumerate(_KEYS):
        assert values[i] == pytest.approx(
            golden.get(i), rel=GOLDEN_RTOL
        ), f"{name}: {values[i]} vs golden {golden.get(i)}"


if __name__ == "__main__":
    # Regenerate data/truth_tables/planck_synthetic_m2lnl_golden.bin.
    import os
    import tempfile

    with tempfile.TemporaryDirectory() as scratch:
        vals = _compute_all(scratch)

    vec = Ncm.Vector.new_array(vals)
    out = os.path.abspath(
        os.path.join(
            os.path.dirname(__file__), "..", "..", "..", "..", "data", GOLDEN_FILE
        )
    )
    Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(vec, out)
    for k, v in zip(_KEYS, vals):
        print(f"  {k:20s} {v:.15g}")
    print(f"wrote {out}")
