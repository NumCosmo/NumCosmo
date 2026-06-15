#!/usr/bin/env python
#
# test_galaxy_hod.py
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
#

"""Tests for the NcGalaxyHOD halo occupation distribution (Zheng07)."""

import math

from numcosmo_py import Nc, Ncm

Ncm.cfg_init()

LOG_MMIN = 12.72
SIGMA_LOG_M = 0.26
LOG_M0 = 12.7
LOG_M1 = 13.93
ALPHA = 1.15


def _ref_mean_central(log10_m: float) -> float:
    """Reference Zheng07 mean central occupation."""
    arg = (log10_m - LOG_MMIN) / SIGMA_LOG_M
    return 0.5 + Ncm.util_normal_gaussian_integral(0.0, arg * math.sqrt(2.0))


def _ref_mean_satellite(log10_m: float) -> float:
    """Reference Zheng07 mean satellite occupation."""
    diff = 10.0**log10_m - 10.0**LOG_M0
    if diff <= 0.0:
        return 0.0
    return (diff / 10.0**LOG_M1) ** ALPHA


def test_means_match_reference() -> None:
    """The C means reproduce the Zheng07 formulas at the default parameters."""
    hod = Nc.GalaxyHODZheng07.new()
    for log10_m in (12.0, 12.72, 13.5, 14.5, 15.0):
        ln_m = log10_m * math.log(10.0)
        assert math.isclose(
            hod.mean_n_central(ln_m), _ref_mean_central(log10_m), rel_tol=1e-12
        )
        assert math.isclose(
            hod.mean_n_satellite(ln_m), _ref_mean_satellite(log10_m), rel_tol=1e-12
        )


def test_is_ncm_model() -> None:
    """Zheng07 is an NcGalaxyHOD and an NcmModel with the five parameters."""
    hod = Nc.GalaxyHODZheng07.new()
    assert isinstance(hod, Nc.GalaxyHOD)
    assert isinstance(hod, Ncm.Model)
    assert math.isclose(hod.param_get_by_name("logMmin"), LOG_MMIN)
    assert math.isclose(hod.param_get_by_name("alpha"), ALPHA)


def test_gen_satellites_require_central() -> None:
    """gen returns a 0/1 central and satellites only when a central is present."""
    hod = Nc.GalaxyHODZheng07.new()
    rng = Ncm.RNG.seeded_new(None, 42)
    ln_m = 14.0 * math.log(10.0)
    for _ in range(2000):
        n_cen, n_sat = hod.gen(ln_m, rng)
        assert n_cen in (0, 1)
        assert n_sat >= 0
        if n_cen == 0:
            assert n_sat == 0


def test_deterministic_central() -> None:
    """With stochastic-central off the central follows the 0.5 threshold."""
    hod = Nc.GalaxyHODZheng07.new()
    hod.set_stochastic_central(False)
    assert not hod.get_stochastic_central()
    rng = Ncm.RNG.seeded_new(None, 1)

    # Below logMmin the mean is < 0.5 -> no central; well above it is 1.
    low_cen, _ = hod.gen(11.0 * math.log(10.0), rng)
    high_cen, _ = hod.gen(14.0 * math.log(10.0), rng)
    assert low_cen == 0
    assert high_cen == 1


def test_stochastic_central_statistics() -> None:
    """The stochastic central rate matches the mean occupation."""
    hod = Nc.GalaxyHODZheng07.new()
    rng = Ncm.RNG.seeded_new(None, 7)
    ln_m = LOG_MMIN * math.log(10.0)  # mean central = 0.5 here
    n = 20000
    centrals = sum(hod.gen(ln_m, rng)[0] for _ in range(n))
    assert abs(centrals / n - 0.5) < 0.02


def test_serialization_roundtrip() -> None:
    """The model survives a serialization round-trip."""
    hod = Nc.GalaxyHODZheng07.new()
    hod.param_set_by_name("alpha", 1.0)
    hod.set_stochastic_central(False)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    dup = ser.from_string(ser.to_string(hod, True))

    assert isinstance(dup, Nc.GalaxyHODZheng07)
    assert math.isclose(dup.param_get_by_name("alpha"), 1.0)
    assert not dup.get_stochastic_central()
