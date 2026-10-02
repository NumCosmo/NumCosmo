#!/usr/bin/env python
#
# test_mass_projection.py
#
# Copyright  2026  Cinthia N. Lima
# <cinthia.n.lima@uel.br>
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

"""Tests for NcClusterMassProjection model."""

import math
import pytest
import numpy as np
from scipy.integrate import quad
from numcosmo_py import Nc, Ncm
from numcosmo_py.helper import duplicate_via_serialization

Ncm.cfg_init()

LN_RICHNESS_CUT = np.log(20.0)
LN_RICHNESS_MAX = np.log(600.0)
M0 = 3.0e14
Z0 = 0.6

LNM = math.log(2.0e14)
Z = 0.4
OBS_PARAMS = [0.0]


def _set_relation(model) -> None:
    """Set the same mass-richness relation on Ascaso and Projection."""
    model.param_set_by_name("mup0", 3.19)
    model.param_set_by_name("mup1", 2.0 / math.log(10.0))
    model.param_set_by_name("mup2", -0.7 / math.log(10.0))
    model.param_set_by_name("sigmap0", 0.33)
    model.param_set_by_name("sigmap1", -0.08 / math.log(10.0))
    model.param_set_by_name("sigmap2", 0.0)
    model.param_set_by_name("cut", LN_RICHNESS_CUT)


@pytest.fixture(name="cluster_m")
def fixture_cluster_m() -> Nc.ClusterMassProjection:
    """Create the mass-richness relation with projection effects."""
    cluster_m = Nc.ClusterMassProjection(
        lnRichness_min=np.log(5.0), lnRichness_max=LN_RICHNESS_MAX, M0=M0, z0=Z0
    )
    _set_relation(cluster_m)
    cluster_m.param_set_by_name("fprj", 0.25)
    cluster_m.param_set_by_name("tau", 0.08)
    return cluster_m


@pytest.fixture(name="cosmo")
def fixture_cosmo() -> Nc.HICosmo:
    """Create a cosmology; the richness model does not depend on it."""
    return Nc.HICosmoDEXcdm()


def test_defaults(cluster_m: Nc.ClusterMassProjection) -> None:
    """Test the model exposes its parameters and tolerance."""
    assert cluster_m.get_cut() == LN_RICHNESS_CUT
    assert cluster_m.f_prj() == pytest.approx(0.25)
    assert cluster_m.tau() == pytest.approx(0.08)
    assert cluster_m.get_reltol() > 0.0


def test_serialization(cluster_m: Nc.ClusterMassProjection) -> None:
    """Test the model survives a serialization round trip."""
    dup = duplicate_via_serialization(cluster_m)
    assert isinstance(dup, Nc.ClusterMassProjection)
    assert dup.f_prj() == pytest.approx(cluster_m.f_prj())
    assert dup.tau() == pytest.approx(cluster_m.tau())
    assert dup.get_cut() == pytest.approx(cluster_m.get_cut())


def test_no_projection_matches_ascaso(cosmo: Nc.HICosmo) -> None:
    """With f_prj = 0 the model must reduce to NcClusterMassAscaso."""
    ascaso = Nc.ClusterMassAscaso(
        lnRichness_min=np.log(5.0), lnRichness_max=LN_RICHNESS_MAX, M0=M0, z0=Z0
    )
    projection = Nc.ClusterMassProjection(
        lnRichness_min=np.log(5.0), lnRichness_max=LN_RICHNESS_MAX, M0=M0, z0=Z0
    )
    _set_relation(ascaso)
    _set_relation(projection)
    projection.param_set_by_name("fprj", 0.0)

    for ln_richness in np.linspace(LN_RICHNESS_CUT, LN_RICHNESS_MAX, 32):
        assert projection.p(cosmo, LNM, Z, [ln_richness], OBS_PARAMS) == pytest.approx(
            ascaso.p(cosmo, LNM, Z, [ln_richness], OBS_PARAMS), abs=1.0e-14
        )

    assert projection.intp(cosmo, LNM, Z) == pytest.approx(
        ascaso.intp(cosmo, LNM, Z), rel=1.0e-12
    )


def test_p_below_cut_vanishes(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """The density must vanish below the richness cut."""
    assert cluster_m.p(cosmo, LNM, Z, [LN_RICHNESS_CUT - 0.5], OBS_PARAMS) == 0.0


def test_projection_adds_richness(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """Projection moves probability to high richness, never to low."""
    no_prj = duplicate_via_serialization(cluster_m)
    no_prj.param_set_by_name("fprj", 0.0)

    high = math.log(120.0)
    assert cluster_m.p(cosmo, LNM, Z, [high], OBS_PARAMS) > no_prj.p(
        cosmo, LNM, Z, [high], OBS_PARAMS
    )


def test_intp_bin_matches_quadrature(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """The closed-form bin integral must match a quadrature of the density."""
    edges = np.log([20.0, 40.0, 80.0, 600.0])
    for lo, hi in zip(edges[:-1], edges[1:]):
        analytic = cluster_m.intp_bin(cosmo, LNM, Z, [lo], [hi], OBS_PARAMS)
        numeric = quad(
            lambda t: cluster_m.p(cosmo, LNM, Z, [t], OBS_PARAMS), lo, hi, limit=300
        )[0]
        assert analytic == pytest.approx(numeric, rel=1.0e-4)


def test_intp_is_sum_of_bins(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """intP over the full range must equal the sum of a partition of it."""
    edges = np.linspace(LN_RICHNESS_CUT, LN_RICHNESS_MAX, 9)
    total = sum(
        cluster_m.intp_bin(cosmo, LNM, Z, [lo], [hi], OBS_PARAMS)
        for lo, hi in zip(edges[:-1], edges[1:])
    )
    assert total == pytest.approx(cluster_m.intp(cosmo, LNM, Z), rel=1.0e-10)


def test_reltol_controls_accuracy(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """A looser tolerance must stay close to a tight one."""
    tight = duplicate_via_serialization(cluster_m)
    tight.set_reltol(1.0e-11)
    loose = duplicate_via_serialization(cluster_m)
    loose.set_reltol(1.0e-5)

    for ln_richness in np.linspace(LN_RICHNESS_CUT, LN_RICHNESS_MAX, 16):
        assert loose.p(cosmo, LNM, Z, [ln_richness], OBS_PARAMS) == pytest.approx(
            tight.p(cosmo, LNM, Z, [ln_richness], OBS_PARAMS), rel=1.0e-2, abs=1.0e-10
        )


def test_resample_matches_density(cluster_m: Nc.ClusterMassProjection, cosmo) -> None:
    """Sampling must reproduce the density: the mock-versus-truth check."""
    rng = Ncm.RNG.seeded_new(None, 1234)
    obs = Ncm.Vector.new(1)
    pars = Ncm.Vector.new_array(OBS_PARAMS)

    nsamples = 60000
    samples = np.empty(nsamples)
    for i in range(nsamples):
        cluster_m.resample_vec(cosmo, LNM, Z, obs, pars, rng)
        samples[i] = obs.get(0)

    inside = samples[(samples >= LN_RICHNESS_CUT) & (samples <= LN_RICHNESS_MAX)]

    # The fraction kept must match intP within Binomial errors.
    fraction = len(inside) / nsamples
    expected = cluster_m.intp(cosmo, LNM, Z)
    sigma = math.sqrt(expected * (1.0 - expected) / nsamples)
    assert abs(fraction - expected) < 4.0 * sigma

    # The shape must match the density, bin by bin, within Poisson errors.
    counts, edges = np.histogram(
        inside, bins=10, range=(LN_RICHNESS_CUT, LN_RICHNESS_MAX)
    )
    for count, lo, hi in zip(counts, edges[:-1], edges[1:]):
        if count < 25:
            continue
        empirical = count / nsamples / (hi - lo)
        predicted = cluster_m.intp_bin(cosmo, LNM, Z, [lo], [hi], OBS_PARAMS) / (hi - lo)
        error = math.sqrt(count) / nsamples / (hi - lo)
        assert abs(empirical - predicted) < 4.0 * error


def test_reltol_survives_serialization(cluster_m: Nc.ClusterMassProjection) -> None:
    """The reltol property must round trip, not pick up another property's value."""
    cluster_m.set_reltol(1.0e-9)
    dup = duplicate_via_serialization(cluster_m)
    assert dup.get_reltol() == pytest.approx(1.0e-9)
