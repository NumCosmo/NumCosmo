#!/usr/bin/env python
#
# test_multiplicity_func_base.py
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

"""Tests for the NcMultiplicityFunc base class.

These cover the behaviour the abstract class implements itself, rather than any
particular fit: the mass-definition conversion of Delta, and the filtered
power-spectrum slot that non-universal models read.
"""

import pytest

from numcosmo_py import Ncm, Nc

Ncm.cfg_init()


@pytest.fixture(name="cosmo")
def fixture_cosmo() -> Nc.HICosmo:
    """Return a simple flat cosmology."""
    cosmo = Nc.HICosmoDEXcdm.new()
    cosmo.omega_x2omega_k()
    cosmo.param_set_by_name("H0", 70.0)
    cosmo.param_set_by_name("Omegac", 0.25)
    cosmo.param_set_by_name("Omegab", 0.05)
    cosmo.param_set_by_name("Omegak", 0.0)
    return cosmo


def _tinker(mdef: Nc.MultiplicityFuncMassDef, delta: float) -> Nc.MultiplicityFunc:
    """Return a Tinker multiplicity function with the requested mass definition."""
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(mdef)
    mulf.set_Delta(delta)
    return mulf


def test_matter_delta_mean_is_delta_itself(cosmo: Nc.HICosmo) -> None:
    """A mean-overdensity Delta is already referred to the matter density."""
    mulf = _tinker(Nc.MultiplicityFuncMassDef.MEAN, 200.0)

    for z in (0.0, 0.5, 2.0):
        assert mulf.get_matter_Delta(cosmo, z) == pytest.approx(200.0, rel=1.0e-14)


def test_matter_delta_critical_divides_by_omega_m(cosmo: Nc.HICosmo) -> None:
    """A critical-overdensity Delta converts with a factor 1 / Omega_m(z)."""
    mulf = _tinker(Nc.MultiplicityFuncMassDef.CRITICAL, 200.0)

    for z in (0.0, 0.5, 2.0):
        omega_m = cosmo.E2Omega_m(z) / cosmo.E2(z)
        expected = 200.0 / omega_m

        assert mulf.get_matter_Delta(cosmo, z) == pytest.approx(expected, rel=1.0e-14)
        # Omega_m(z) < 1 here, so the matter-referred Delta is always the larger.
        assert mulf.get_matter_Delta(cosmo, z) > 200.0


def test_matter_delta_critical_approaches_mean_as_matter_dominates(
    cosmo: Nc.HICosmo,
) -> None:
    """The two definitions converge as Omega_m(z) grows towards one.

    They never coincide exactly here: radiation still holds a few percent of the
    budget at z ~ 100, so this checks the gap shrinks monotonically rather than
    asserting a limit the cosmology does not actually reach.
    """
    crit = _tinker(Nc.MultiplicityFuncMassDef.CRITICAL, 200.0)
    mean = _tinker(Nc.MultiplicityFuncMassDef.MEAN, 200.0)

    gaps = [
        crit.get_matter_Delta(cosmo, z) / mean.get_matter_Delta(cosmo, z) - 1.0
        for z in (0.0, 1.0, 5.0, 20.0)
    ]

    for prev, cur in zip(gaps, gaps[1:]):
        assert cur < prev
    assert gaps[-1] < 0.05


def test_psf_slot_defaults_to_unset() -> None:
    """Universal models never set a filter, so the slot starts empty."""
    mulf = Nc.MultiplicityFuncTinker.new()

    assert mulf.peek_psf() is None
    assert mulf.get_property("powerspectrum-filtered") is None


def test_psf_slot_round_trip() -> None:
    """The filter can be set, read back, and cleared again."""
    ps = Nc.PowspecMLCBE.new()
    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    mulf = Nc.MultiplicityFuncTinker.new()

    mulf.set_psf(psf)
    assert mulf.peek_psf() == psf

    # Setting the same object again is a no-op rather than a double ref.
    mulf.set_psf(psf)
    assert mulf.peek_psf() == psf

    mulf.set_psf(None)
    assert mulf.peek_psf() is None


def test_psf_slot_via_property() -> None:
    """The filter is also reachable through the GObject property."""
    ps = Nc.PowspecMLCBE.new()
    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.GAUSS)
    mulf = Nc.MultiplicityFuncTinker.new()

    mulf.set_property("powerspectrum-filtered", psf)
    assert mulf.get_property("powerspectrum-filtered") == psf
    assert mulf.peek_psf() == psf


def test_mass_function_injects_its_own_filter() -> None:
    """NcHaloMassFunction wires its filter into the multiplicity function.

    This is what keeps sigma_R and the slope read by a non-universal model on the
    same power spectrum, instead of relying on the caller to set both.
    """
    ps = Nc.PowspecMLCBE.new()
    psf = Ncm.PowspecFilter.new(ps, Ncm.PowspecFilterType.TOPHAT)
    mulf = Nc.MultiplicityFuncCastro.new()
    dist = Nc.Distance.new(5.0)

    assert mulf.peek_psf() is None

    mfp = Nc.HaloMassFunction.new(dist, psf, mulf)

    assert mulf.peek_psf() == psf
    assert mfp.peek_psf() == psf
