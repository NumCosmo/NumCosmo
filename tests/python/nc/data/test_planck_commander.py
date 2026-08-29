#!/usr/bin/env python
#
# test_planck_commander.py
#
# Wed July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_commander.py
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

"""Tests on the native NcDataPlanckCommander (low-ell gibbs_gauss) likelihood."""

import numpy as np
import pytest

from python.fixtures_planck import (
    COMMANDER_NL,
    FixedClBoltzmann,
    commander_tables,
    make_commander_cldf,
    planck_mset,
    theory_cls,
)
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH, build_commander
from numcosmo_py.experiments.planck_lite import find_baseline_file

Ncm.cfg_init()

_CLIK = find_baseline_file(COMMANDER_RELPATH)
needs_data = pytest.mark.skipif(
    _CLIK is None, reason=f"commander baseline data not found ({COMMANDER_RELPATH})"
)


def test_type_is_ncmdata():
    """The native type is an NcmData subclass."""
    assert issubclass(Nc.DataPlanckCommander, Ncm.Data)


@needs_data
def test_construct_and_serialize_roundtrip():
    """Build from the commander clik data, serialize and reload intact."""
    cmd = build_commander(_CLIK)  # data-only (no Boltzmann)
    assert cmd.get_length() == 28

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    dset = Ncm.Dataset.new_array([cmd])
    cmd2 = ser.from_variant(ser.to_variant(dset)).peek_data(0)

    assert isinstance(cmd2, Nc.DataPlanckCommander)
    assert cmd2.get_length() == 28


@pytest.mark.planck_data
@needs_data
def test_matches_clik_reference():
    """Native commander m2lnL agrees with the clik low-ell TT reference.

    Both share one CBE Boltzmann + cosmology + PlanckFICorTT (for A_planck).
    The native Gaussianized-Blackwell-Rao assembly (spline transform, Gaussian
    in x-space, spline Jacobian and offset) reproduces clik to ~1e-7.
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

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe)
    native = build_commander(_CLIK, cbe)

    ref.prepare(mset)
    native.prepare(mset)
    m2_ref = ref.m2lnL_val(mset)
    m2_nat = native.m2lnL_val(mset)

    assert m2_nat == pytest.approx(m2_ref, rel=1.0e-6)

    rng = Ncm.RNG.seeded_new(None, 42)
    native.resample(mset, rng)
    assert np.isfinite(native.m2lnL_val(mset))


@pytest.mark.planck_data
@needs_data
def test_clik_pi_compat_matches_clik():
    """clik_pi_compat lands the native m2lnL on clik, to the last few ULP.

    clik's gibbs uses a single-precision pi in the Dl = Cl*l(l+1)/2pi conversion.
    Reproducing it closes the ~1e-7 gap the (more accurate) default leaves: the
    compat path agrees with clik to ~5e-15, five orders of magnitude closer, which
    is what this pins. Equality is not asserted -- the two implementations reach
    the same value by different floating-point routes, so the last bits are a
    property of the build, not of the code.
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

    ref = Nc.DataPlanckLKL.full_new_id(Nc.DataPlanckLKLType.BASELINE_18_LOWL_TT, cbe)
    native = build_commander(_CLIK, cbe, clik_pi_compat=True)

    ref.prepare(mset)
    native.prepare(mset)

    m2_ref = ref.m2lnL_val(mset)
    m2_default = build_commander(_CLIK, cbe).m2lnL_val(mset)

    assert native.m2lnL_val(mset) == pytest.approx(m2_ref, rel=1.0e-12)
    # The default (double-pi) path stays measurably further out, so the tolerance
    # above cannot be passing merely by being loose.
    assert abs(m2_default / m2_ref - 1.0) > 1.0e-9


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


def _commander_lnl(dl_vec, xa, ya, mu, cov):
    """Closed-form gibbs_gauss lnL for the synthetic (exactly linear) transform."""
    slope = (ya[0, -1] - ya[0, 0]) / (xa[0, -1] - xa[0, 0])
    x = ya[:, 0] + (dl_vec - xa[:, 0]) * slope
    d = x - mu
    chi2 = d @ np.linalg.inv(cov) @ d

    return -0.5 * chi2 + dl_vec.size * np.log(slope)


def test_synthetic_matches_closed_form(tmp_path):
    """The native m2lnL reproduces the Gaussianized Blackwell-Rao closed form.

    The synthetic ``sigma.fits`` tabulates a straight Cl->x line (zero second
    derivatives), so the spline transform, its Jacobian and the mu_sigma offset
    all have closed forms; matching them checks the whole chain, including the
    covariance inversion and the A_planck rescaling.
    """
    clik = make_commander_cldf(tmp_path)
    pb = FixedClBoltzmann()
    calib = 1.003
    mset, _ = planck_mset(A_planck=calib)

    cmd = build_commander(clik, pb)
    assert cmd.get_length() == COMMANDER_NL

    cmd.prepare(mset)

    xa, ya, _, mu, cov, mu_sigma = commander_tables()
    ell = np.arange(2, 2 + COMMANDER_NL)
    cl = theory_cls(pb, "TT", 1 + COMMANDER_NL)[2:]
    dl = cl / calib**2 * ell * (ell + 1.0) / (2.0 * np.pi)

    expected = -2.0 * (
        _commander_lnl(dl, xa, ya, mu, cov) - _commander_lnl(mu_sigma, xa, ya, mu, cov)
    )
    assert cmd.m2lnL_val(mset) == pytest.approx(expected, rel=1.0e-12)


def test_synthetic_clik_pi_compat_shifts_result(tmp_path):
    """clik_pi_compat changes the Dl conversion, and so the answer, slightly."""
    clik = make_commander_cldf(tmp_path)
    pb = FixedClBoltzmann()
    mset, _ = planck_mset()

    exact = build_commander(clik, pb)
    compat = build_commander(clik, pb, clik_pi_compat=True)
    assert compat.get_property("clik-pi-compat")
    assert not exact.get_property("clik-pi-compat")

    exact.prepare(mset)
    compat.prepare(mset)

    m2_exact = exact.m2lnL_val(mset)
    m2_compat = compat.m2lnL_val(mset)
    assert m2_exact != m2_compat
    assert m2_compat == pytest.approx(m2_exact, rel=1.0e-5)


def test_synthetic_out_of_prior_range(tmp_path):
    """A theory Dl outside the tabulated range is rejected, not extrapolated."""
    clik = make_commander_cldf(tmp_path)
    pb = FixedClBoltzmann()
    # A tiny calibration blows Dl = Cl/A^2 * l(l+1)/2pi past the table's top end.
    mset, _ = planck_mset(A_planck=0.01)

    cmd = build_commander(clik, pb)
    cmd.prepare(mset)

    assert cmd.m2lnL_val(mset) == 1.0e30


def test_synthetic_serialize_roundtrip(tmp_path):
    """A synthetic commander survives a serialization round trip unchanged."""
    clik = make_commander_cldf(tmp_path)
    pb = FixedClBoltzmann()
    mset, _ = planck_mset()

    cmd = build_commander(clik, pb)
    cmd.prepare(mset)
    m2lnl = cmd.m2lnL_val(mset)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    cmd2 = ser.from_variant(ser.to_variant(cmd))
    assert isinstance(cmd2, Nc.DataPlanckCommander)

    cmd2.set_hipert_boltzmann(pb)
    cmd2.prepare(mset)
    assert cmd2.m2lnL_val(mset) == m2lnl


def test_synthetic_new_from_file(tmp_path):
    """The `new_from_file` loader reads back a serialized commander.

    Only the compact likelihoods are checked this way: the loader is shared
    across all five, and the text serialization of the high-l blocks is large.
    """
    clik = make_commander_cldf(tmp_path)
    data = build_commander(clik)

    path = str(tmp_path / "commander.obj")
    Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP).to_file(data, path)

    loaded = Nc.DataPlanckCommander.new_from_file(path)
    assert isinstance(loaded, Nc.DataPlanckCommander)
    assert loaded.get_length() == data.get_length()
