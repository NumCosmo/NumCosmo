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
