#!/usr/bin/env python
#
# test_hipert_boltzmann_cbe.py
#
# Sat August 29 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_hipert_boltzmann_cbe.py
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

"""Tests for the CLASS-backed Boltzmann, #NcHIPertBoltzmannCBE.

The CMB likelihoods are tested against stored spectra (see
tests/python/fixtures_planck.py), which deliberately takes CLASS out of their
loop. This module is the other half of that split: it checks that the CBE
actually produces sane spectra, on its own and cheaply.

The assertions are physical invariants rather than reference values -- positive
auto-spectra, the Cauchy-Schwarz bound on TE, the first acoustic peak, and
lensing moving the damping tail -- so this stays a fast sanity check that does
not need a golden table or a tight tolerance. A multipole range well below
Planck's is enough for all of them, which is what keeps it cheap.
"""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc
from numcosmo_py.cosmology import create_cosmo, HIPrimModel

Ncm.cfg_init()

LMAX = 1000
_SPECTRA = Nc.DataCMBDataType.TT | Nc.DataCMBDataType.EE | Nc.DataCMBDataType.TE


def _cbe(lensed: bool, phiphi: bool = True):
    """A CBE targeting TT/EE/TE (+PHIPHI) up to LMAX."""
    pb = Nc.HIPertBoltzmannCBE.new()
    target = _SPECTRA | Nc.DataCMBDataType.PHIPHI if phiphi else _SPECTRA
    pb.set_target_Cls(target)
    for name in ("TT", "EE", "TE"):
        getattr(pb, f"set_{name}_lmax")(LMAX)
    if phiphi:
        pb.set_PHIPHI_lmax(LMAX)
    pb.set_lensed_Cls(lensed)

    return pb


def _cls(pb, name: str) -> np.ndarray:
    """The prepared $C_\\ell$ of @name as a numpy array."""
    vec = Ncm.Vector.new(LMAX + 1)
    getattr(pb, f"get_{name}_Cls")(vec)

    return np.array([vec.get(i) for i in range(LMAX + 1)])


@pytest.fixture(name="lensed_cls", scope="module")
def fixture_lensed_cls():
    """One lensed solve shared by the whole module (one CLASS run)."""
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    pb = _cbe(lensed=True)
    pb.prepare(cosmo)

    return {name: _cls(pb, name) for name in ("TT", "EE", "TE", "PHIPHI")}


def test_auto_spectra_are_positive(lensed_cls):
    """TT, EE and phi-phi are positive definite over the whole range."""
    for name in ("TT", "EE", "PHIPHI"):
        cl = lensed_cls[name][2:]
        assert np.all(np.isfinite(cl)), f"{name} has non-finite entries"
        assert np.all(cl > 0.0), f"{name} is not positive"


def test_te_obeys_cauchy_schwarz(lensed_cls):
    """|C_l^TE| <= sqrt(C_l^TT C_l^EE): TE is a cross-spectrum of the same fields.

    A tolerance-free invariant -- it holds for any cosmology, so it catches a
    mis-scaled or mis-indexed TE without pinning a reference value.
    """
    tt, ee, te = (lensed_cls[n][2:] for n in ("TT", "EE", "TE"))

    assert np.all(np.abs(te) <= np.sqrt(tt * ee))


def test_te_changes_sign(lensed_cls):
    """TE oscillates about zero, unlike the auto-spectra."""
    te = lensed_cls["TE"][2:]

    assert np.any(te > 0.0) and np.any(te < 0.0)


def test_first_acoustic_peak(lensed_cls):
    """D_l^TT peaks near l ~ 220, the standard LCDM first acoustic peak."""
    ell = np.arange(LMAX + 1, dtype=float)
    dl = np.zeros(LMAX + 1)
    dl[2:] = lensed_cls["TT"][2:] * ell[2:] * (ell[2:] + 1.0) / (2.0 * np.pi)

    assert 180 <= int(np.argmax(dl[2:400])) + 2 <= 260


def test_lensing_smooths_the_damping_tail():
    """Lensing leaves the large scales alone and fills in the high-l tail.

    Two solves of the same cosmology differing only in lensed_Cls, which is the
    switch the CMB likelihoods rely on.
    """
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)

    unlensed_pb = _cbe(lensed=False)
    unlensed_pb.prepare(cosmo)
    unlensed = _cls(unlensed_pb, "TT")

    lensed_pb = _cbe(lensed=True)
    lensed_pb.prepare(cosmo)
    lensed = _cls(lensed_pb, "TT")

    # Large scales are essentially untouched by lensing.
    assert lensed[2:30] == pytest.approx(unlensed[2:30], rel=1.0e-2)
    # The damping tail is not: lensing transfers power into the troughs.
    assert np.max(np.abs(lensed[800:] / unlensed[800:] - 1.0)) > 1.0e-2
