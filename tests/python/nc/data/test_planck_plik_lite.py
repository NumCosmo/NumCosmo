#!/usr/bin/env python
#
# test_planck_plik_lite.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_plik_lite.py
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

"""Tests on the native NcDataPlanckPlikLite likelihood."""

import numpy as np
import pytest

from python.fixtures_planck import (
    PLIK_LITE_LMAX,
    make_plik_lite_cldf,
    model_vector,
    planck_mset,
    plik_lite_tables,
    theory_cls,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import (
    PLIK_LITE_TT_RELPATH,
    PLIK_LITE_TTTEEE_RELPATH,
    NBIN_TT,
    NBIN_TOTAL,
    PLMIN,
    NBIN_TE,
    NBIN_EE,
    find_baseline_file,
    build_plik_lite,
    build_plik_lite_tt,
)

Ncm.cfg_init()

_CLIK = find_baseline_file(PLIK_LITE_TT_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"plik_lite baseline data not found ({PLIK_LITE_TT_RELPATH})"
)


def test_type_is_resampleable_gausscov():
    """The native type is an NcmDataGaussCov (so it inherits resample)."""
    plik = Nc.DataPlanckPlikLite()
    assert isinstance(plik, Ncm.DataGaussCov)
    assert isinstance(plik, Ncm.Data)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from clik data, serialize the dataset to binary, reload intact."""
    plik = build_plik_lite_tt(_CLIK)  # no Boltzmann: data-only
    assert plik.get_size() == NBIN_TT

    mean = Ncm.DataGaussCov.peek_mean(plik)
    cov = Ncm.DataGaussCov.peek_cov(plik)

    dset = Ncm.Dataset.new_array([plik])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    variant = ser.to_variant(dset)
    ser.clear_instances(False)
    dset2 = ser.from_variant(variant)
    plik2 = dset2.peek_data(0)

    assert isinstance(plik2, Nc.DataPlanckPlikLite)
    assert plik2.get_size() == NBIN_TT
    assert plik2.get_property("calib-name") == "A_planck"

    mean2 = Ncm.DataGaussCov.peek_mean(plik2)
    cov2 = Ncm.DataGaussCov.peek_cov(plik2)
    for i in range(NBIN_TT):
        assert mean2.get(i) == mean.get(i)
        for j in range(0, NBIN_TT, 13):
            assert cov2.get(i, j) == cov.get(i, j)


@pytest.mark.app
@needs_data
def test_matches_clik_reference():
    """Native m2lnL equals the clik plik_lite reference on a CLASS cosmology.

    Both likelihoods share one CBE Boltzmann + cosmology + PlanckFI, so they see
    identical theory Cl. The clik reference reproduces the file's check_value, so
    an exact match validates the native reimplementation. Also checks resample.

    Uses the TTTEEE file: it exercises the full TT+TE+EE code path (TT is its
    first block). Only one clik plik_lite can be built per process (the Fortran
    ``plik_cmbonly`` keeps global module state), so a single superset case is used.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import (
        mset_set_parameters,
        Planck18Types,
    )

    clik_path = find_baseline_file(PLIK_LITE_TTTEEE_RELPATH)
    if clik_path is None:
        pytest.skip(f"baseline data not found ({PLIK_LITE_TTTEEE_RELPATH})")

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()  # plik_lite's only nuisance param is A_planck
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TTTEEE_LITE, cbe
    )
    native = build_plik_lite(clik_path, cbe)
    assert native.get_size() == NBIN_TOTAL

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, abs=1.0e-6)

    # resample: the capability the clik wrapper lacks.
    rng = Ncm.RNG.seeded_new(None, 123)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


# Bins per spectrum, in the file's data-vector order (TT, TE, EE).
_NBIN = {"TT": NBIN_TT, "TE": NBIN_TE, "EE": NBIN_EE}


def _expected_bandpowers(cbe, spectra, calib):
    """Closed-form bandpower model: binned raw Cl over each bin / A_planck^2."""
    _, _, blmin, blmax, bweight = plik_lite_tables()
    expected = []
    for name in spectra:
        nbin = _NBIN[name]
        cl = theory_cls(cbe, name, PLIK_LITE_LMAX)
        for b in range(nbin):
            lo, hi = int(blmin[b]), int(blmax[b])
            window = bweight[lo : hi + 1]
            expected.append(np.sum(cl[lo + PLMIN : hi + PLMIN + 1] * window) / calib**2)

    return np.array(expected)


@pytest.mark.parametrize("spectra", [["TT"], ["TT", "TE", "EE"]])
def test_synthetic_matches_closed_form(tmp_path, spectra):
    """The binned model bandpowers match the closed form, spectrum by spectrum.

    Exercises both the TT-only and the full TT/TE/EE layouts: the converter
    picks the active blocks out of the 613-long data vector and covariance from
    ``has_cl``, and the native mean_func must bin the matching theory spectrum.
    """
    clik = make_plik_lite_cldf(tmp_path, spectra=spectra)
    cbe = Nc.HIPertBoltzmannCBE.new()
    calib = 1.004
    mset, _ = planck_mset(A_planck=calib)

    plik = build_plik_lite(clik, cbe, lmax=PLIK_LITE_LMAX)
    assert plik.get_size() == sum(_NBIN[name] for name in spectra)

    plik.prepare(mset)

    assert model_vector(plik, mset) == pytest.approx(
        _expected_bandpowers(cbe, spectra, calib), rel=1.0e-13
    )


def test_synthetic_spectra_override(tmp_path):
    """An explicit @spectra list overrides the file's has_cl selection."""
    clik = make_plik_lite_cldf(tmp_path, spectra=["TT", "TE", "EE"])
    cbe = Nc.HIPertBoltzmannCBE.new()
    mset, _ = planck_mset()

    plik = build_plik_lite(clik, cbe, spectra=["EE"], lmax=PLIK_LITE_LMAX)
    assert plik.get_size() == NBIN_EE

    plik.prepare(mset)
    assert model_vector(plik, mset) == pytest.approx(
        _expected_bandpowers(cbe, ["EE"], 1.0), rel=1.0e-13
    )


def test_synthetic_tt_helper_and_resample(tmp_path):
    """build_plik_lite_tt selects TT only, and the result resamples."""
    clik = make_plik_lite_cldf(tmp_path, spectra=["TT", "TE", "EE"])
    cbe = Nc.HIPertBoltzmannCBE.new()
    mset, _ = planck_mset()

    plik = build_plik_lite_tt(clik, cbe, lmax=PLIK_LITE_LMAX)
    assert plik.get_size() == NBIN_TT

    plik.prepare(mset)
    assert np.isfinite(plik.m2lnL_val(mset))

    rng = Ncm.RNG.seeded_new(None, 17)
    plik.resample(mset, rng)
    assert np.isfinite(plik.m2lnL_val(mset))


def test_synthetic_serialize_roundtrip(tmp_path):
    """A synthetic plik_lite survives a serialization round trip unchanged."""
    clik = make_plik_lite_cldf(tmp_path)
    cbe = Nc.HIPertBoltzmannCBE.new()
    mset, _ = planck_mset()

    plik = build_plik_lite(clik, cbe, lmax=PLIK_LITE_LMAX)
    plik.prepare(mset)
    m2lnl = plik.m2lnL_val(mset)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    plik2 = ser.from_variant(ser.to_variant(plik))
    assert isinstance(plik2, Nc.DataPlanckPlikLite)
    assert plik2.peek_hipert_boltzmann() is not None

    plik2.prepare(mset)
    assert plik2.m2lnL_val(mset) == m2lnl
