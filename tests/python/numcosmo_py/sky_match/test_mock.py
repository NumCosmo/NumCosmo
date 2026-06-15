#!/usr/bin/env python
#
# test_mock.py
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

"""Tests for the mock catalog generator's cosmology/geometry helpers."""

import pytest
import numpy as np

pytest.importorskip("astropy")
# flake8: noqa: E402
# pylint: disable=wrong-import-position

from astropy.table import Table

from numcosmo_py import Nc, Ncm
from numcosmo_py.catalog import (
    MockGenerator,
    ConstantCompleteness,
    ConstantPurity,
    identity_scaling_relation,
)
from numcosmo_py.cosmology import Cosmology

Ncm.cfg_init()

LOG_M_MIN = 12.72  # HOD_model default logMmin


def _cosmo() -> Nc.HICosmo:
    """A simple flat XCDM cosmology for tests."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.props.H0 = 70.0
    cosmo.props.Omegac = 0.25
    cosmo.props.Omegab = 0.05
    cosmo.props.Omegax = 0.70
    return cosmo


def _mock_with_hmf(**kwargs) -> MockGenerator:
    """A MockGenerator wired with a real halo mass function and cluster mass.

    Small footprint and a narrow high-mass range keep the abundance sampling fast
    and the realization modest (~100 halos).
    """
    cosmology = Cosmology.default()
    mulf = Nc.MultiplicityFuncTinker.new()
    mulf.set_mdef(Nc.MultiplicityFuncMassDef.CRITICAL)
    hmf = Nc.HaloMassFunction.new(cosmology.dist, cosmology.psf_tophat, mulf)
    cluster_m = Nc.ClusterMassNodist(lnM_min=np.log(1e14), lnM_max=np.log(1e15))
    params = dict(
        ra_interval=(-2.0, 2.0),
        dec_interval=(-2.0, 2.0),
        z_interval=(0.2, 0.5),
        halo_mass_interval=(1e14, 1e15),
        seed=42,
    )
    params.update(kwargs)
    return MockGenerator(cosmo=cosmology.cosmo, hmf=hmf, cluster_m=cluster_m, **params)


def test_rho_crit_matches_numcosmo() -> None:
    """rho_crit(z) == rho_crit0 * E2(z), and grows with z (matter ~ (1+z)^3)."""
    cosmo = _cosmo()
    mg = MockGenerator(cosmo=cosmo)
    rho_crit0 = Ncm.C.crit_mass_density_h2_solar_mass_Mpc3() * cosmo.h2()

    z = np.array([0.0, 0.5, 1.0, 2.0])
    expected = rho_crit0 * np.array([cosmo.E2(zi) for zi in z])
    assert np.allclose(np.asarray(mg.rho_crit(z)), expected, rtol=1e-12)

    # Critical density must increase with redshift (the previous bug made it decrease).
    rho = np.asarray(mg.rho_crit(z))
    assert np.all(np.diff(rho) > 0.0)


def test_get_r200c_matches_spherical_overdensity_definition() -> None:
    """R200c == (3 M / (4 pi 200 rho_crit(z)))^(1/3)."""
    cosmo = _cosmo()
    mg = MockGenerator(cosmo=cosmo)

    mass = np.array([1e13, 1e14, 1e15])
    z = np.array([0.1, 0.3, 0.7])
    expected = ((3.0 * mass) / (4.0 * np.pi * 200.0 * np.asarray(mg.rho_crit(z)))) ** (
        1.0 / 3.0
    )
    assert np.allclose(np.asarray(mg.get_R200c(mass, z)), expected, rtol=1e-12)
    # More massive halos are larger at fixed redshift.
    r_fixed_z = np.asarray(mg.get_R200c(mass, np.full_like(mass, 0.3)))
    assert np.all(np.diff(r_fixed_z) > 0.0)


def test_generate_clusters_reproducible_and_within_footprint() -> None:
    """Same seed reproduces the catalog; positions/redshifts respect the bounds."""
    cosmo = _cosmo()
    box = dict(
        ra_interval=(-10.0, 10.0),
        dec_interval=(-5.0, 5.0),
        z_interval=(0.2, 0.5),
        cluster_mass_interval=(1e14, 1e15),
        cluster_set_size=200,
    )
    clusters_a = MockGenerator(cosmo=cosmo, seed=42, **box).generate_clusters()
    clusters_b = MockGenerator(cosmo=cosmo, seed=42, **box).generate_clusters()

    # Reproducible across instances with the same seed.
    for col in ("RA", "DEC", "z", "Mass"):
        assert np.array_equal(clusters_a[col], clusters_b[col])

    # Positions and redshifts within the requested footprint.
    assert np.all(clusters_a["RA"] >= -10.0) and np.all(clusters_a["RA"] <= 10.0)
    assert np.all(clusters_a["DEC"] >= -5.0) and np.all(clusters_a["DEC"] <= 5.0)
    assert np.all(clusters_a["z"] >= 0.2) and np.all(clusters_a["z"] <= 0.5)

    # R200c column is consistent with the helper.
    mg = MockGenerator(cosmo=cosmo, **box)
    expected_r200c = mg.get_R200c(
        10.0 ** np.asarray(clusters_a["Mass"]), clusters_a["z"]
    )
    assert np.allclose(clusters_a["R200c"], np.asarray(expected_r200c), rtol=1e-12)


def test_sky_area_is_rectangular_solid_angle() -> None:
    """sky_area is the solid angle (deg^2) of the RA/DEC rectangle."""
    mg = MockGenerator(
        cosmo=_cosmo(), ra_interval=(-10.0, 10.0), dec_interval=(-5.0, 5.0)
    )
    dra = np.radians(10.0 - (-10.0))
    dsin = np.sin(np.radians(5.0)) - np.sin(np.radians(-5.0))
    expected = dra * dsin * (180.0 / np.pi) ** 2
    assert mg.sky_area() == pytest.approx(expected, rel=1e-12)


def test_get_3d_coordinates_matches_astro_convention() -> None:
    """get_3D_coordinates agrees with NumCosmo's TriVec astro convention (degrees)."""
    mg = MockGenerator(cosmo=_cosmo())
    for ra, dec, r in [(30.0, 10.0, 100.0), (-45.0, -20.0, 50.0), (120.0, 60.0, 200.0)]:
        x1, x2, x3 = mg.get_3D_coordinates(ra, dec, r)
        ref = Ncm.TriVec.new_astro_ra_dec(r, ra, dec).c
        assert np.allclose([x1, x2, x3], list(ref), rtol=0, atol=1e-9)


def test_hod_model_deterministic_central_is_mass_dependent() -> None:
    """With stochastic_central=False, the central is present only above threshold."""
    mg = MockGenerator(cosmo=_cosmo(), stochastic_central=False)
    # Well below logMmin -> mean_n_cen ~ 0 -> no central; well above -> central.
    n_cen_low, _ = mg.HOD_model(10**11.0)
    n_cen_high, _ = mg.HOD_model(10**14.0)
    assert n_cen_low == 0
    assert n_cen_high == 1


def test_hod_model_stochastic_central_fraction() -> None:
    """With stochastic_central=True, the realized central fraction tracks <N_cen>."""
    mg = MockGenerator(cosmo=_cosmo(), stochastic_central=True)
    # At m = 10**logMmin, <N_cen> = 0.5 exactly.
    np.random.seed(1234)
    draws = np.array([mg.HOD_model(10**LOG_M_MIN)[0] for _ in range(5000)])
    assert set(np.unique(draws)).issubset({0, 1})
    assert draws.mean() == pytest.approx(0.5, abs=0.03)


def test_hod_model_per_call_override() -> None:
    """The per-call stochastic_central overrides the instance default."""
    mg = MockGenerator(cosmo=_cosmo(), stochastic_central=True)
    # Deterministic override at low mass -> no central, regardless of instance default.
    assert mg.HOD_model(10**11.0, stochastic_central=False)[0] == 0


def test_generate_galaxies_empty_returns_empty_table() -> None:
    """No galaxies generated returns an empty Table instead of raising."""
    mg = MockGenerator(cosmo=_cosmo(), stochastic_central=False)
    catalog = Table()
    catalog["Mass"] = [10.0, 10.5]  # log10 mass, far below logMmin -> no members
    catalog["RA"] = [0.0, 1.0]
    catalog["DEC"] = [0.0, 1.0]
    catalog["z"] = [0.3, 0.3]
    catalog["R200c"] = [0.5, 0.5]

    result = mg.generate_galaxies(catalog, "cluster")
    assert isinstance(result, Table)
    assert len(result) == 0


def test_generate_halos_from_hmf_invariants() -> None:
    """HMF-sampled halo catalog respects footprint, mass cut and derived columns."""
    mg = _mock_with_hmf()
    halos = mg.generate_halos_from_hmf()

    assert len(halos) > 0
    assert set(halos.colnames) == {
        "halo_id",
        "RA",
        "DEC",
        "z",
        "Mass",
        "Mass_obs",
        "R200c",
        "x1",
        "x2",
        "x3",
        "halo_r",
        "is_detected",
    }

    # Footprint and mass cut.
    assert np.all(halos["RA"] >= -2.0) and np.all(halos["RA"] <= 2.0)
    assert np.all(halos["DEC"] >= -2.0) and np.all(halos["DEC"] <= 2.0)
    assert np.all(halos["z"] >= 0.2) and np.all(halos["z"] <= 0.5)
    assert np.all(halos["Mass"] >= np.log(1e14) - 1e-9)
    assert np.all(halos["Mass"] <= np.log(1e15) + 1e-9)

    # No completeness model -> everything detected.
    assert np.all(np.asarray(halos["is_detected"]) == 1)

    # Derived columns consistent with the helpers.
    expected_r200c = mg.get_R200c(np.exp(np.asarray(halos["Mass"])), halos["z"])
    assert np.allclose(halos["R200c"], np.asarray(expected_r200c), rtol=1e-12)
    x1, x2, x3 = mg.get_3D_coordinates(halos["RA"], halos["DEC"], halos["halo_r"])
    assert np.allclose(halos["x1"], x1, rtol=1e-12)
    assert np.allclose(halos["x2"], x2, rtol=1e-12)
    assert np.allclose(halos["x3"], x3, rtol=1e-12)


def test_generate_halos_from_hmf_reproducible() -> None:
    """Same seed reproduces the HMF-sampled catalog across instances."""
    halos_a = _mock_with_hmf().generate_halos_from_hmf()
    halos_b = _mock_with_hmf().generate_halos_from_hmf()
    assert len(halos_a) == len(halos_b)
    for col in ("RA", "DEC", "z", "Mass", "Mass_obs"):
        assert np.array_equal(halos_a[col], halos_b[col])


def test_generate_halos_from_hmf_completeness() -> None:
    """A completeness model yields a binary, non-trivial detection flag."""
    halos = _mock_with_hmf().generate_halos_from_hmf(ConstantCompleteness(0.5))
    detected = np.asarray(halos["is_detected"])
    assert set(np.unique(detected)).issubset({0, 1})
    assert 0 < detected.mean() < 1


def test_generate_clusters_from_halos_real_only() -> None:
    """Without a purity model, every cluster is a real detection of a parent halo."""
    mg = _mock_with_hmf()
    halos = mg.generate_halos_from_hmf()
    n_detected = int(np.sum(np.asarray(halos["is_detected"]) == 1))

    clusters = mg.generate_clusters_from_halos(halos)

    assert set(clusters.colnames) == {
        "cluster_id",
        "RA",
        "DEC",
        "z",
        "Mass_obs",
        "Mass",
        "R200c",
        "x1",
        "x2",
        "x3",
        "cluster_r",
        "parent_id",
    }
    assert len(clusters) == n_detected
    assert mg.cluster_set_size == n_detected
    # All real: parent_id non-zero and drawn from the detected halo ids.
    assert np.all(np.asarray(clusters["parent_id"]) != 0)
    assert set(np.asarray(clusters["parent_id"])).issubset(
        set(np.asarray(halos["halo_id"]))
    )
    # R200c consistent with the helper (uses true mass).
    expected = mg.get_R200c(np.exp(np.asarray(clusters["Mass"])), clusters["z"])
    assert np.allclose(clusters["R200c"], np.asarray(expected), rtol=1e-12)


def test_generate_clusters_from_halos_purity_injects_fakes() -> None:
    """A purity model appends fake (parent_id == 0) clusters inside the footprint."""
    mg = _mock_with_hmf()
    halos = mg.generate_halos_from_hmf()
    n_detected = int(np.sum(np.asarray(halos["is_detected"]) == 1))

    clusters = mg.generate_clusters_from_halos(
        halos, 2.0, ConstantPurity(0.9), identity_scaling_relation
    )

    parent = np.asarray(clusters["parent_id"])
    n_real = int(np.sum(parent != 0))
    n_fake = int(np.sum(parent == 0))
    assert n_real == n_detected
    assert n_fake > 0
    assert len(clusters) == n_detected + n_fake
    assert mg.cluster_set_size == len(clusters)

    # Fakes are drawn uniformly within the footprint.
    fakes = clusters[parent == 0]
    assert np.all(fakes["RA"] >= -2.0) and np.all(fakes["RA"] <= 2.0)
    assert np.all(fakes["DEC"] >= -2.0) and np.all(fakes["DEC"] <= 2.0)
    assert np.all(fakes["z"] >= 0.2) and np.all(fakes["z"] <= 0.5)


def test_generate_clusters_from_halos_reproducible() -> None:
    """Same seed reproduces the cluster catalog (real + fakes) across instances."""

    def run():
        mg = _mock_with_hmf()
        halos = mg.generate_halos_from_hmf()
        return mg.generate_clusters_from_halos(
            halos, 2.0, ConstantPurity(0.9), identity_scaling_relation
        )

    a, b = run(), run()
    assert len(a) == len(b)
    for col in ("RA", "DEC", "z", "Mass", "Mass_obs", "parent_id"):
        assert np.array_equal(a[col], b[col])
