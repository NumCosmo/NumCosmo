#!/usr/bin/env python
#
# test_powspec_filter.py
#
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

"""Tests for NcmPowspecFilter."""

import numpy as np
import pytest

from numcosmo_py import Ncm, Nc

# Required for this file to run at all: conftest.py skips anything whose keywords
# contain "powspec" -- which the directory name alone supplies -- unless
# --run-powspec is given, while the py-powspec lane selects with `-m powspec`,
# which matches markers only. Without this the file is skipped in one lane and
# deselected in the other.
pytestmark = pytest.mark.powspec

Ncm.cfg_init()

# A power law is the one case with a closed-form answer: for P(k) = A k^n the
# filtered variance is exactly sigma^2(R) ~ R^-(n+3), independent of the window,
# so dlnvar/dlnr = -(n+3) and every higher log-derivative vanishes identically.
# That makes it an external check on the transform rather than a comparison of
# the filter against itself.
POWER_LAW_INDICES = [-2.0, -1.5]
FILTER_TYPES = [Ncm.PowspecFilterType.TOPHAT, Ncm.PowspecFilterType.GAUSS]


def _power_law_filter(n_index: float, ftype: Ncm.PowspecFilterType, nderivs: int):
    """Build a filter over a tabulated pure power law P(k) = A k^n."""
    ka = np.geomspace(1.0e-6, 1.0e3, 3000)
    za = np.linspace(0.0, 1.0, 10, dtype=np.float64)
    za_v, ka_v = np.meshgrid(za, ka)
    lnPk = np.log(1.0e-9 * ka_v**n_index) + 0.0 * za_v

    Pk2d = Ncm.Spline2dBicubic(
        spline=Ncm.SplineCubicNotaknot.new(),
        x_vector=Ncm.Vector.new_array(za.tolist()),
        y_vector=Ncm.Vector.new_array(np.log(ka).tolist()),
        z_matrix=Ncm.Matrix.new_array(lnPk.flatten().tolist(), len(za)),
    )
    ps = Ncm.PowspecSpline2d.new(Pk2d)
    ps.prepare()

    psf = Ncm.PowspecFilter.new(ps, ftype)
    psf.require_nderivs(nderivs)
    psf.set_reltol(1.0e-6)
    psf.set_best_lnr0()
    psf.prepare(Nc.HICosmoLCDM.new())

    return psf


def test_nderivs_default() -> None:
    """The filter computes the first derivative unless more is required."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)

    assert psf.get_nderivs() == 1


def test_require_nderivs_never_lowers() -> None:
    """require_nderivs is a ratchet: it raises the order and never lowers it."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)
    assert psf.get_nderivs() == 1

    psf.require_nderivs(3)
    assert psf.get_nderivs() == 3

    # A smaller requirement must not discard the order somebody else needs.
    psf.require_nderivs(1)
    assert psf.get_nderivs() == 3

    psf.require_nderivs(3)
    assert psf.get_nderivs() == 3


def test_set_nderivs_lowers_and_stays_usable() -> None:
    """set_nderivs may lower the order, discarding the tables above it.

    This is the one way down: require_nderivs is a ratchet. The discarded
    orders must actually be released and the filter must still evaluate the
    orders it kept, after re-preparing.
    """
    cosmo = Nc.HICosmoLCDM.new()
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 3)
    assert psf.get_nderivs() == 3

    lnr = np.log(5.0)
    var_before = psf.eval_dnvar_dlnrn(0.0, lnr, 0)

    psf.set_nderivs(1)
    assert psf.get_nderivs() == 1

    # Lowering invalidates the calibration, so the filter must be prepared again.
    psf.prepare(cosmo)

    assert psf.eval_dnvar_dlnrn(0.0, lnr, 0) == pytest.approx(var_before, rel=1.0e-12)
    assert psf.eval_dnlnvar_dlnrn(0.0, lnr, 1) == pytest.approx(-1.5, abs=1.0e-4)


def test_nderivs_property_round_trip() -> None:
    """The order is readable and writable through the GObject property."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)
    assert psf.props.nderivs == 1

    psf.props.nderivs = 2
    assert psf.props.nderivs == 2
    assert psf.get_nderivs() == 2


@pytest.mark.parametrize("ftype", FILTER_TYPES)
@pytest.mark.parametrize("n_index", POWER_LAW_INDICES)
def test_power_law_first_log_derivative(
    n_index: float, ftype: Ncm.PowspecFilterType
) -> None:
    """dlnvar/dlnr must equal -(n+3) for a power-law spectrum."""
    psf = _power_law_filter(n_index, ftype, 1)
    lnr = np.linspace(np.log(1.0), np.log(20.0), 9)

    d1 = np.array([psf.eval_dnlnvar_dlnrn(0.0, x, 1) for x in lnr])

    assert_allclose_msg = f"filter {ftype}, n = {n_index}"
    np.testing.assert_allclose(
        d1, -(n_index + 3.0), atol=1.0e-4, err_msg=assert_allclose_msg
    )


@pytest.mark.parametrize("ftype", FILTER_TYPES)
@pytest.mark.parametrize("n_index", POWER_LAW_INDICES)
def test_power_law_second_log_derivative(
    n_index: float, ftype: Ncm.PowspecFilterType
) -> None:
    """d2lnvar/dlnr2 must vanish for a power-law spectrum.

    This is the order that used to be obtained by differentiating the bicubic
    interpolation of the first derivative; taking it from the transform itself
    is worth one to two orders of magnitude here.
    """
    psf = _power_law_filter(n_index, ftype, 2)
    lnr = np.linspace(np.log(1.0), np.log(20.0), 9)

    d2 = np.array([psf.eval_dnlnvar_dlnrn(0.0, x, 2) for x in lnr])

    np.testing.assert_allclose(
        d2, 0.0, atol=1.0e-4, err_msg=f"filter {ftype}, n = {n_index}"
    )


def test_eval_r_forms_match_lnr_forms() -> None:
    """The r-argument evaluators are the lnr ones at the same point."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 2)

    for r in (0.5, 2.0, 9.0):
        lnr = np.log(r)

        assert psf.eval_var(0.0, r) == pytest.approx(
            psf.eval_var_lnr(0.0, lnr), rel=1.0e-14
        )
        assert psf.eval_sigma(0.0, r) == pytest.approx(
            psf.eval_sigma_lnr(0.0, lnr), rel=1.0e-14
        )
        assert psf.eval_lnvar_lnr(0.0, lnr) == pytest.approx(
            np.log(psf.eval_var_lnr(0.0, lnr)), rel=1.0e-14
        )
        # sigma is the square root of the variance, by definition.
        assert psf.eval_sigma_lnr(0.0, lnr) == pytest.approx(
            np.sqrt(psf.eval_var_lnr(0.0, lnr)), rel=1.0e-14
        )


def test_dvar_and_dlnvar_relations() -> None:
    """The raw, log and per-r derivative forms are consistent at a point."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 2)

    for r in (0.7, 3.0, 11.0):
        lnr = np.log(r)
        var = psf.eval_var_lnr(0.0, lnr)
        dvar = psf.eval_dvar_dlnr(0.0, lnr)

        assert dvar == pytest.approx(psf.eval_dnvar_dlnrn(0.0, lnr, 1), rel=1.0e-14)
        assert psf.eval_dlnvar_dlnr(0.0, lnr) == pytest.approx(dvar / var, rel=1.0e-14)
        # d/dr = (1/r) d/dlnr.
        assert psf.eval_dlnvar_dr(0.0, lnr) == pytest.approx(
            psf.eval_dlnvar_dlnr(0.0, lnr) / r, rel=1.0e-14
        )
        # n = 0 of the log-derivative ladder is ln(var) itself.
        assert psf.eval_dnlnvar_dlnrn(0.0, lnr, 0) == pytest.approx(
            np.log(var), rel=1.0e-14
        )


@pytest.mark.parametrize("ftype", FILTER_TYPES)
def test_volume_rm3_is_filter_specific(ftype: Ncm.PowspecFilterType) -> None:
    """Each window has its own V/r^3, positive and independent of the spectrum."""
    psf = _power_law_filter(-1.5, ftype, 1)

    volume = psf.volume_rm3()
    assert volume > 0.0

    other = _power_law_filter(-2.0, ftype, 1).volume_rm3()
    assert volume == pytest.approx(other, rel=1.0e-14)


def test_volume_rm3_differs_between_filters() -> None:
    """A top-hat and a Gaussian of the same radius enclose different volumes."""
    tophat = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1).volume_rm3()
    gauss = _power_law_filter(-1.5, Ncm.PowspecFilterType.GAUSS, 1).volume_rm3()

    assert tophat != pytest.approx(gauss, rel=1.0e-6)


def test_peek_powspec_returns_the_source() -> None:
    """The filter hands back the spectrum it was built on."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)

    ps = psf.peek_powspec()
    assert isinstance(ps, Ncm.Powspec)
    assert psf.get_filter_type() == Ncm.PowspecFilterType.TOPHAT


def test_reltol_accessors_round_trip() -> None:
    """Both calibration tolerances are settable and readable."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)

    psf.set_reltol(1.0e-8)
    psf.set_reltol_z(1.0e-7)

    assert psf.get_reltol() == pytest.approx(1.0e-8)
    assert psf.get_reltol_z() == pytest.approx(1.0e-7)


def test_zi_zf_require_is_a_ratchet() -> None:
    """require_zi/require_zf widen the range and never narrow it."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)

    psf.set_zi(0.5)
    psf.set_zf(3.0)
    assert psf.props.zi == pytest.approx(0.5)
    assert psf.props.zf == pytest.approx(3.0)

    # Widening in either direction takes effect.
    psf.require_zi(0.1)
    psf.require_zf(6.0)
    assert psf.props.zi == pytest.approx(0.1)
    assert psf.props.zf == pytest.approx(6.0)

    # Narrowing requests are ignored, so a shared filter stays usable.
    psf.require_zi(2.0)
    psf.require_zf(1.0)
    assert psf.props.zi == pytest.approx(0.1)
    assert psf.props.zf == pytest.approx(6.0)


def test_lnr0_moves_the_output_range() -> None:
    """Re-centring the output shifts the r range it can be evaluated on."""
    psf = _power_law_filter(-1.5, Ncm.PowspecFilterType.TOPHAT, 1)

    r_min_before, r_max_before = psf.get_r_min(), psf.get_r_max()
    assert 0.0 < r_min_before < r_max_before

    psf.set_lnr0(psf.props.lnr0 + np.log(10.0))

    assert psf.props.lnr0 == pytest.approx(
        np.log(10.0) + np.log(r_min_before * r_max_before) / 2.0, rel=1.0e-6
    )
    assert psf.get_r_min() == pytest.approx(10.0 * r_min_before, rel=1.0e-10)
    assert psf.get_r_max() == pytest.approx(10.0 * r_max_before, rel=1.0e-10)

    # set_best_lnr0 puts it back on the range the spectrum actually supports.
    psf.set_best_lnr0()
    assert psf.get_r_min() == pytest.approx(r_min_before, rel=1.0e-10)


@pytest.mark.parametrize("ftype", FILTER_TYPES)
def test_raw_derivatives_match_log_derivatives(ftype: Ncm.PowspecFilterType) -> None:
    """The log-derivatives must follow from the raw ones at the same point."""
    psf = _power_law_filter(-1.5, ftype, 2)
    lnr = np.linspace(np.log(1.0), np.log(20.0), 9)

    for x in lnr:
        var = psf.eval_dnvar_dlnrn(0.0, x, 0)
        dvar = psf.eval_dnvar_dlnrn(0.0, x, 1)
        d2var = psf.eval_dnvar_dlnrn(0.0, x, 2)

        assert var == pytest.approx(psf.eval_var_lnr(0.0, x), rel=1.0e-14)
        assert dvar / var == pytest.approx(
            psf.eval_dnlnvar_dlnrn(0.0, x, 1), rel=1.0e-12
        )
        assert d2var / var - (dvar / var) ** 2 == pytest.approx(
            psf.eval_dnlnvar_dlnrn(0.0, x, 2), rel=1.0e-12
        )
