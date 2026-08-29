#!/usr/bin/env python
#
# test_planck_lensing.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_lensing.py
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

"""Tests on the native NcDataPlanckLensing (Planck 2018 CMB lensing) likelihood."""

import numpy as np
import pytest

from python.fixtures_planck import (
    LENSING_LMAX,
    LENSING_NBINS,
    lensing_tables,
    make_lensing_cldf,
    model_vector,
    planck_mset,
    theory_cls,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lensing import (
    LENSING_FULL_RELPATH,
    LENSING_MARGED_RELPATH,
    build_lensing,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

# (relpath, DataPlanckLKL enum, whether the file needs CMB spectra) per flavor.
_FLAVORS = {
    "marged": (
        LENSING_MARGED_RELPATH,
        Nc.DataPlanckLKLType.BASELINE_18_LENSING_CMB_MARGINALIZED,
    ),
    "full": (
        LENSING_FULL_RELPATH,
        Nc.DataPlanckLKLType.BASELINE_18_LENSING,
    ),
}


def _clik(relpath):
    return find_baseline_file(relpath)


def _mset():
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()
    return mset


def test_type_is_gauss():
    """The native type is an NcmDataGauss subclass."""
    assert issubclass(Nc.DataPlanckLensing, Ncm.DataGauss)


@pytest.mark.parametrize("flavor", list(_FLAVORS))
def test_construct_and_serialize_roundtrip(flavor):
    """Build from the lensing clik data, serialize and reload intact."""
    relpath, _ = _FLAVORS[flavor]
    if _clik(relpath) is None:
        pytest.skip(f"lensing baseline data not found ({relpath})")

    lens = build_lensing(_clik(relpath))  # data-only (no Boltzmann)
    assert lens.get_length() == 9

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([lens])
    lens2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(lens2, Nc.DataPlanckLensing)
    assert lens2.get_length() == 9


@pytest.mark.planck_data
@pytest.mark.parametrize("flavor", list(_FLAVORS))
def test_matches_clik_reference(flavor):
    """Native lensing m2lnL agrees with the clik lensing reference.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT (for A_planck).
    The native inverse-covariance Gaussian assembly of the band-power model
    reproduces the clik dts_lensing likelihood to machine precision.
    """
    relpath, enum = _FLAVORS[flavor]
    if _clik(relpath) is None:
        pytest.skip(f"lensing baseline data not found ({relpath})")

    cbe = Nc.HIPertBoltzmannCBE.new()
    mset = _mset()

    ref = Nc.DataPlanckLKL.full_new_id(enum, cbe)
    native = build_lensing(_clik(relpath), cbe)

    ref.prepare(mset)
    native.prepare(mset)

    assert native.m2lnL_val(mset) == pytest.approx(ref.m2lnL_val(mset), rel=1.0e-11)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@pytest.mark.planck_data
def test_prepare_self_configures_boltzmann():
    """Lensing raises the Boltzmann targets/lmax itself on a bare CBE.

    Unlike plik_lite/smica/commander/simall (which rely on the caller to set the
    targets and lmax), NcDataPlanckLensing self-configures in prepare(). Drive a
    fresh CBE that only has the default TT target at a low lmax and check the full
    file lifts PHIPHI + TT/EE/TE up to its own lmax and still evaluates.
    """
    if _clik(LENSING_FULL_RELPATH) is None:
        pytest.skip(f"lensing baseline data not found ({LENSING_FULL_RELPATH})")

    cbe = Nc.HIPertBoltzmannCBE.new()  # default: TT only, lmax ~30
    mset = _mset()

    native = build_lensing(_clik(LENSING_FULL_RELPATH), cbe)
    native.prepare(mset)

    assert cbe.get_PHIPHI_lmax() >= 2500
    assert cbe.get_TT_lmax() >= 2500
    assert cbe.get_EE_lmax() >= 2500
    assert cbe.get_TE_lmax() >= 2500
    assert np.isfinite(native.m2lnL_val(mset))


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


def _expected_bandpowers(cbe, marged, calib, has_calib=True):
    """Closed-form band-power model of the synthetic lensing file."""
    hascl, bins, cor0, cors, _, _ = lensing_tables(marged)
    nl = LENSING_LMAX + 1
    ell = np.arange(nl)
    w_pp = ell**2 * (ell + 1.0) ** 2 / (2.0 * np.pi)
    w_cmb = ell * (ell + 1.0) / (2.0 * np.pi)

    phi = theory_cls(cbe, "PHIPHI", LENSING_LMAX)
    expected = bins @ (phi * w_pp) - cor0 + cors[:, :nl] @ (phi * w_pp)

    # hascl flags [TT, EE, BB, TE, TB, EB]; blocks follow in that order, each
    # rescaled by 1/A_planck^2 when the file declares a calibration.
    scale = 1.0 / calib**2 if has_calib else 1.0
    names = [n for n, f in zip(("TT", "EE", "BB", "TE", "TB", "EB"), hascl) if f]
    for k, name in enumerate(names):
        block = cors[:, (k + 1) * nl : (k + 2) * nl]
        expected += block @ (theory_cls(cbe, name, LENSING_LMAX) * w_cmb) * scale

    return expected


@pytest.mark.parametrize("marged", [False, True])
def test_synthetic_matches_closed_form(tmp_path, marged):
    """The assembled band-powers match the closed-form projection.

    Covers both flavors: the full file renormalizes with the CMB spectra its
    ``hascl`` selects (scaled by the calibration), the CMB-marginalized one
    carries the phi block alone.
    """
    clik = make_lensing_cldf(tmp_path, marged=marged)
    cbe = Nc.HIPertBoltzmannCBE.new()
    calib = 1.0  # exercised separately below
    mset, _ = planck_mset(A_planck=calib)

    lens = build_lensing(clik, cbe)
    assert lens.get_length() == LENSING_NBINS

    lens.prepare(mset)
    model = model_vector(lens, mset)

    assert model == pytest.approx(_expected_bandpowers(cbe, marged, calib), rel=1.0e-12)


def test_synthetic_calibration_scales_cmb_block(tmp_path):
    """A_planck rescales the CMB renormalization block, not the phi projection."""
    clik = make_lensing_cldf(tmp_path)
    cbe = Nc.HIPertBoltzmannCBE.new()

    mset_one, _ = planck_mset(A_planck=1.0)
    mset_cal, _ = planck_mset(A_planck=1.05)

    lens = build_lensing(clik, cbe)
    lens.prepare(mset_one)
    model_one = model_vector(lens, mset_one)
    model_cal = model_vector(lens, mset_cal)

    assert model_cal != pytest.approx(model_one, rel=1.0e-10)
    assert model_cal == pytest.approx(
        _expected_bandpowers(cbe, False, 1.05), rel=1.0e-12
    )


def test_synthetic_no_calibration(tmp_path):
    """With has_calib off the CMB block is not rescaled by A_planck."""
    clik = make_lensing_cldf(tmp_path, has_calib=False)
    cbe = Nc.HIPertBoltzmannCBE.new()

    lens = build_lensing(clik, cbe)
    assert not lens.get_property("has-calib")

    mset_one, _ = planck_mset(A_planck=1.0)
    mset_cal, _ = planck_mset(A_planck=1.05)
    lens.prepare(mset_one)

    assert model_vector(lens, mset_cal) == pytest.approx(
        model_vector(lens, mset_one), rel=1.0e-14
    )


def test_synthetic_resample_and_serialize(tmp_path):
    """A synthetic lensing block resamples and survives serialization."""
    clik = make_lensing_cldf(tmp_path)
    cbe = Nc.HIPertBoltzmannCBE.new()
    mset, _ = planck_mset()

    lens = build_lensing(clik, cbe)
    lens.prepare(mset)
    assert np.isfinite(lens.m2lnL_val(mset))

    rng = Ncm.RNG.seeded_new(None, 321)
    lens.resample(mset, rng)
    assert np.isfinite(lens.m2lnL_val(mset))

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    lens2 = ser.from_variant(ser.to_variant(lens))
    assert isinstance(lens2, Nc.DataPlanckLensing)

    lens2.set_hipert_boltzmann(cbe)
    lens2.prepare(mset)
    assert lens2.m2lnL_val(mset) == lens.m2lnL_val(mset)


def test_synthetic_new_from_file(tmp_path):
    """The `new_from_file` loader reads back a serialized lensing.

    Only the compact likelihoods are checked this way: the loader is shared
    across all five, and the text serialization of the high-l blocks is large.
    """
    clik = make_lensing_cldf(tmp_path)
    data = build_lensing(clik)

    path = str(tmp_path / "lensing.obj")
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).to_file(data, path)

    loaded = Nc.DataPlanckLensing.new_from_file(path)
    assert isinstance(loaded, Nc.DataPlanckLensing)
    assert loaded.get_length() == data.get_length()
