#!/usr/bin/env python
#
# test_planck_smica.py
#
# Fri June 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_smica.py
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

"""Tests on the native NcDataPlanckSmica (plik high-ell SMICA) likelihood."""

import numpy as np
import pytest

from python.fixtures_planck import (
    SMICA_LMIN,
    SMICA_LMAX,
    SMICA_NBINS,
    make_smica_cldf,
    model_vector,
    planck_mset,
    smica_binning,
    smica_mask_and_ordering,
    theory_cls,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_smica import (
    PLIK_TT_RELPATH,
    PLIK_TTTEEE_RELPATH,
    build_smica_tt,
    build_smica_ttteee,
)
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_CLIK = find_baseline_file(PLIK_TT_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"plik baseline data not found ({PLIK_TT_RELPATH})"
)

_CLIK_TTTEEE = find_baseline_file(PLIK_TTTEEE_RELPATH)
needs_data_ttteee = pytest.mark.skipif(
    _CLIK_TTTEEE is None,
    reason=f"plik baseline data not found ({PLIK_TTTEEE_RELPATH})",
)


def test_type_is_resampleable_gauss():
    """The native type is an NcmDataGauss (so it inherits resample).

    Asserted on the type, not an instance: NcDataPlanckSmica is fully specified
    at construction (all structural properties are CONSTRUCT_ONLY and validated
    in constructed()), so it is deliberately not default-constructible.
    """
    assert issubclass(Nc.DataPlanckSmica, Ncm.DataGauss)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from the plik clik data, serialize the dataset, reload intact."""
    smica = build_smica_tt(_CLIK)  # data-only (no Boltzmann)
    npt = smica.get_size()
    assert npt == 765

    mean = Ncm.DataGauss.peek_mean(smica)
    dset = Ncm.Dataset.new_array([smica])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.get_size() == npt
    mean2 = Ncm.DataGauss.peek_mean(smica2)
    for i in range(0, npt, 37):
        assert mean2.get(i) == mean.get(i)


@pytest.mark.planck_data
@needs_data
def test_matches_clik_reference():
    """Native SMICA m2lnL agrees with the clik plik TT reference on a CLASS cosmology.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT. The native
    reimplementation reproduces the full R_q assembly (CMB + 10 foreground/
    calibration components) to machine precision.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TT, cbe)
    native = build_smica_tt(_CLIK, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-12)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@needs_data_ttteee
def test_ttteee_construct_and_serialize_roundtrip():
    """Build the TTTEEE (m=6) SMICA from clik data; serialize/reload intact."""
    smica = build_smica_ttteee(_CLIK_TTTEEE)
    npt = smica.get_size()
    assert npt == 2289

    mean = Ncm.DataGauss.peek_mean(smica)
    dset = Ncm.Dataset.new_array([smica])
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.get_size() == npt
    mean2 = Ncm.DataGauss.peek_mean(smica2)
    for i in range(0, npt, 53):
        assert mean2.get(i) == mean.get(i)


@pytest.mark.planck_data
@needs_data_ttteee
def test_ttteee_matches_clik_reference():
    """Native TTTEEE (m=6) SMICA m2lnL agrees with the clik reference.

    Same NcDataPlanckSmica class as TT, now with the polarization channels and
    components (CMB TT/EE/TE, pwfe galactic dust, EE end-to-end noise, the
    icalTP two-term calibration mixing and totcalP). Machine precision.
    """
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    cbe = Nc.HIPertBoltzmannCBE.new()
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTTTEEE()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TTTEEE, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    ref = Nc.DataPlanckLKL.full_new_id(
        Nc.DataPlanckLKLType.BASELINE_18_HIGHL_TTTEEE, cbe
    )
    native = build_smica_ttteee(_CLIK_TTTEEE, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-12)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


@pytest.fixture(name="synthetic_cbe", scope="module")
def fixture_synthetic_cbe():
    """One CBE shared by the synthetic tests (they all use the same cosmology)."""
    return Nc.HIPertBoltzmannCBE.new()


@pytest.fixture(name="synthetic_tt_clik", scope="module")
def fixture_synthetic_tt_clik(tmp_path_factory):
    """A synthetic TT plik tree, written once for the whole module."""
    return make_smica_cldf(tmp_path_factory.mktemp("plik_tt"))


@pytest.fixture(name="synthetic_ttteee_clik", scope="module")
def fixture_synthetic_ttteee_clik(tmp_path_factory):
    """A synthetic TTTEEE plik tree, written once for the whole module."""
    return make_smica_cldf(tmp_path_factory.mktemp("plik_ttteee"), pol=True)


# Nuisance parameters that switch every foreground component off. The point
# sources have no such setting -- clik defaults every temperature pair to one --
# so they are given known values and kept in the expectation below; the beam
# leakage has no amplitude either, and is switched off through a zero template
# in the synthetic file.
_CMB_ONLY = {
    "A_cib_217": 0.0,
    "cib_index": -1.3,
    "xi_sz_cib": 0.0,
    "A_sz": 0.0,
    "ksz_norm": 0.0,
    "gal545_A_100": 0.0,
    "gal545_A_143": 0.0,
    "gal545_A_143_217": 0.0,
    "gal545_A_217": 0.0,
    "A_sbpx_100_100_TT": 0.0,
    "A_sbpx_143_143_TT": 0.0,
    "A_sbpx_143_217_TT": 0.0,
    "A_sbpx_217_217_TT": 0.0,
    "ps_A_100_100": 10.0,
    "ps_A_143_143": 20.0,
    "ps_A_143_217": 30.0,
    "ps_A_217_217": 40.0,
    "calib_100T": 1.0,
    "calib_217T": 1.0,
    "A_planck": 1.0,
}


def _expected_cmb_only(cbe, m):
    """Closed-form R_q with every foreground off: binned CMB plus point sources."""
    bin_lmin, bin_lmax, _ = smica_binning()
    _, _, quad = smica_mask_and_ordering(m)

    cl = theory_cls(cbe, "TT", SMICA_LMAX)
    band = np.array(
        [
            np.mean(cl[SMICA_LMIN + lo : SMICA_LMIN + hi + 1])
            for lo, hi in zip(bin_lmin, bin_lmax)
        ]
    )

    # clik normalizes the flat point-source amplitudes at ell = 3000.
    nrm = 3000.0 * 3001.0 / (2.0 * np.pi)
    ps = np.full((m, m), 1.0 / nrm)
    ps[0, 0] = _CMB_ONLY["ps_A_100_100"] / nrm
    ps[1, 1] = _CMB_ONLY["ps_A_143_143"] / nrm
    ps[1, 2] = ps[2, 1] = _CMB_ONLY["ps_A_143_217"] / nrm
    ps[2, 2] = _CMB_ONLY["ps_A_217_217"] / nrm

    return (band[:, None, None] + ps[None, :, :]).ravel()[quad]


def test_synthetic_cmb_only_closed_form(synthetic_tt_clik, synthetic_cbe):
    """With every foreground off, R_q reduces to the binned CMB plus point sources.

    This pins the parts of the assembly that have a closed form -- the CMB
    mixing, the bandpower binning operator, the point-source normalization, the
    calibration factors and the masked-entry extraction -- to machine precision,
    which the real (opaque) plik data cannot do.
    """
    clik, cbe = synthetic_tt_clik, synthetic_cbe
    mset, _ = planck_mset(**_CMB_ONLY)

    smica = build_smica_tt(clik, cbe)
    smica.prepare(mset)

    assert model_vector(smica, mset) == pytest.approx(
        _expected_cmb_only(cbe, 3), rel=1.0e-13
    )


def test_synthetic_calibration_rescales_model(synthetic_tt_clik, synthetic_cbe):
    """calib_100T/calib_217T/A_planck rescale R_q by the documented factors.

    R_q[i,j] *= cal_i cal_j / A_planck^2 with cal = 1/sqrt(calib_fT) (143 is the
    reference channel), applied after the whole component assembly.
    """
    clik, cbe = synthetic_tt_clik, synthetic_cbe
    mset, _ = planck_mset(**_CMB_ONLY)

    smica = build_smica_tt(clik, cbe)
    smica.prepare(mset)
    plain = model_vector(smica, mset)

    calib_100t, calib_217t, a_planck = 1.02, 0.97, 1.01
    mset_cal, _ = planck_mset(
        **dict(
            _CMB_ONLY,
            calib_100T=calib_100t,
            calib_217T=calib_217t,
            A_planck=a_planck,
        )
    )
    scaled = model_vector(smica, mset_cal)

    cal = np.array([1.0 / np.sqrt(calib_100t), 1.0, 1.0 / np.sqrt(calib_217t)])
    factor = (np.outer(cal, cal) / a_planck**2).ravel()
    _, _, quad = smica_mask_and_ordering(3)
    # quad indexes a flat (bin, channel, channel) array, so the channel pair of
    # each masked entry is its index modulo m*m.
    per_entry = factor[quad % 9]

    assert scaled == pytest.approx(plain * per_entry, rel=1.0e-13)


def test_synthetic_all_foregrounds_evaluate(synthetic_tt_clik, synthetic_cbe):
    """With every foreground amplitude on, the full assembly stays finite.

    Drives the components with no closed form here (CIB, tSZ, kSZ, CIBxtSZ,
    galactic dust with its pivot normalization, subpixel and beam leakage) and
    checks they actually move the model.
    """
    clik, cbe = synthetic_tt_clik, synthetic_cbe

    mset_off, _ = planck_mset(**_CMB_ONLY)
    mset_on, _ = planck_mset(
        **dict(
            _CMB_ONLY,
            A_cib_217=48.5,
            cib_index=-1.3,
            xi_sz_cib=0.32,
            A_sz=7.03,
            ksz_norm=0.5,
            gal545_A_100=8.86,
            gal545_A_143=10.8,
            gal545_A_143_217=19.43,
            gal545_A_217=94.8,
            A_sbpx_100_100_TT=1.0,
            A_sbpx_143_143_TT=1.0,
            A_sbpx_143_217_TT=1.0,
            A_sbpx_217_217_TT=1.0,
        )
    )

    smica = build_smica_tt(clik, cbe)
    smica.prepare(mset_off)

    model_off = model_vector(smica, mset_off)
    model_on = model_vector(smica, mset_on)

    assert np.all(np.isfinite(model_on))
    assert np.all(model_on != model_off)
    assert np.isfinite(smica.m2lnL_val(mset_on))


def test_synthetic_ttteee(synthetic_ttteee_clik, synthetic_cbe):
    """The TTTEEE (six-channel) build evaluates, resamples and serializes.

    Covers the polarization-only paths: the T/E field map, the TE/EE CMB mixing,
    the galactic power-law dust components, the end-to-end EE correlated noise
    and the two-term icalTP calibration mixing.
    """
    clik, cbe = synthetic_ttteee_clik, synthetic_cbe
    mset, _ = planck_mset(pol=True)

    smica = build_smica_ttteee(clik, cbe)
    assert smica.get_length() == len(smica_mask_and_ordering(6)[2])

    smica.prepare(mset)
    assert np.all(np.isfinite(model_vector(smica, mset)))

    m2lnl = smica.m2lnL_val(mset)
    assert np.isfinite(m2lnl)

    rng = Ncm.RNG.seeded_new(None, 99)
    smica.resample(mset, rng)
    assert np.isfinite(smica.m2lnL_val(mset))

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(smica))
    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.get_length() == smica.get_length()


def test_synthetic_serialize_roundtrip(synthetic_tt_clik, synthetic_cbe):
    """A synthetic TT SMICA survives a serialization round trip unchanged."""
    clik, cbe = synthetic_tt_clik, synthetic_cbe
    mset, _ = planck_mset()

    smica = build_smica_tt(clik, cbe)
    smica.prepare(mset)
    m2lnl = smica.m2lnL_val(mset)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    smica2 = ser.from_variant(ser.to_variant(smica))
    assert isinstance(smica2, Nc.DataPlanckSmica)
    assert smica2.peek_hipert_boltzmann() is not None

    smica2.prepare(mset)
    assert smica2.m2lnL_val(mset) == m2lnl
