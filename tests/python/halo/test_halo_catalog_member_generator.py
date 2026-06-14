#!/usr/bin/env python
#
# test_halo_catalog_member_generator.py
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

"""Tests for NcHaloCatalogMemberGenerator."""

import math

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm, GLib
from numcosmo_py.catalog import catalog_to_table

Ncm.cfg_init()

HOST_COLS = ["ra", "dec", "z_true", "lnM_true", "r_Delta"]


def _cosmo_mset():
    """A minimal cosmology and mset for the member generator."""
    cosmo = Nc.HICosmoDEXcdm()
    cosmo.add_submodel(Nc.HIReionCamb())
    cosmo.add_submodel(Nc.HIPrimPowerLaw())
    for name, value in (
        ("H0", 70.0),
        ("Omegab", 0.05),
        ("Omegac", 0.25),
        ("Omegax", 0.70),
        ("Tgamma0", 2.72),
        ("w", -1.0),
    ):
        cosmo.param_set_by_name(name, value)
    return cosmo, Ncm.MSet.new_array([cosmo])


def _host_catalog(rows):
    """Build a HALO host catalog from a list of (ra, dec, z, lnM, r_Delta)."""
    host = Nc.HaloCatalog.new(
        Nc.HaloCatalogKind.HALO, None, None, len(rows), HOST_COLS, None
    )
    data = host.peek_data()
    for i, row in enumerate(rows):
        for j, val in enumerate(row):
            data.set(i, j, val)
    return host


def _angular_sep_deg(ra1, dec1, ra2, dec2):
    """Great-circle separation (degrees) between two sky points in degrees."""
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    cos_sep = math.sin(d1) * math.sin(d2) + math.cos(d1) * math.cos(d2) * math.cos(
        r1 - r2
    )
    return math.degrees(math.acos(min(1.0, max(-1.0, cos_sep))))


def test_generate_basic_schema() -> None:
    """The result is a MEMBER catalog with the expected columns and linkage."""
    cosmo, mset = _cosmo_mset()
    host = _host_catalog(
        [
            (100.0, 20.0, 0.3, math.log(1.0e15), 2.0),
            (150.0, -10.0, 0.5, math.log(5.0e14), 1.5),
        ]
    )

    hod = Nc.GalaxyHODZheng07.new()
    memgen = Nc.HaloCatalogMemberGenerator.new(hod)
    assert memgen.peek_hod() is hod

    rng = Ncm.RNG.seeded_new(None, 1)
    members = memgen.generate(host, mset, rng)

    assert isinstance(members, Nc.HaloCatalog)
    assert members.get_kind() == Nc.HaloCatalogKind.MEMBER
    assert members.peek_id_col() == "galaxy_id"
    assert members.peek_parent_id_col() == "parent_id"
    assert list(members.peek_columns()) == [
        "galaxy_id",
        "parent_id",
        "is_central",
        "ra",
        "dec",
        "z",
    ]
    assert members.len() > 0

    # galaxy ids are unique and sequential from 0.
    table = catalog_to_table(members)
    assert list(table["galaxy_id"]) == list(range(members.len()))
    # every parent id points at a real host row.
    assert set(np.unique(table["parent_id"])).issubset({0, 1})


def test_central_sits_at_host_and_satellites_inside() -> None:
    """Centrals coincide with the host; satellites stay within r_Delta in 3D."""
    cosmo, mset = _cosmo_mset()
    ra_c, dec_c, z_c, r_delta = 100.0, 20.0, 0.3, 2.0
    host = _host_catalog([(ra_c, dec_c, z_c, math.log(1.0e15), r_delta)])

    hod = Nc.GalaxyHODZheng07.new()
    hod.set_stochastic_central(False)  # massive host -> central guaranteed
    memgen = Nc.HaloCatalogMemberGenerator.new(hod)

    rng = Ncm.RNG.seeded_new(None, 7)
    members = memgen.generate(host, mset, rng)
    assert members is not None
    table = catalog_to_table(members)

    centrals = table[np.asarray(table["is_central"], dtype=bool)]
    assert len(centrals) == 1
    assert math.isclose(float(centrals["ra"][0]), ra_c, abs_tol=1e-12)
    assert math.isclose(float(centrals["dec"][0]), dec_c, abs_tol=1e-12)
    assert math.isclose(float(centrals["z"][0]), z_c, abs_tol=1e-12)

    # Bounds for satellites: angular sep <= r_Delta / dA, |dz| <= H(z) r_Delta / c.
    dist = Nc.Distance.new(2.0)
    dist.prepare(cosmo)
    d_a = dist.angular_diameter(cosmo, z_c) * cosmo.RH_Mpc()
    sep_max_deg = math.degrees(r_delta / d_a)
    dz_max = cosmo.H(z_c) / (Ncm.C.c() / 1.0e3) * r_delta

    sats = table[~np.asarray(table["is_central"], dtype=bool)]
    assert len(sats) > 0
    for row in sats:
        assert _angular_sep_deg(row["ra"], row["dec"], ra_c, dec_c) <= sep_max_deg * (
            1.0 + 1e-9
        )
        assert abs(float(row["z"]) - z_c) <= dz_max * (1.0 + 1e-9)


def test_find_children_links_to_host() -> None:
    """Members can be grouped back to their host through the parent linkage."""
    cosmo, mset = _cosmo_mset()
    host = _host_catalog(
        [
            (100.0, 20.0, 0.3, math.log(1.0e15), 2.0),
            (150.0, -10.0, 0.4, math.log(1.0e15), 2.0),
        ]
    )

    hod = Nc.GalaxyHODZheng07.new()
    hod.set_stochastic_central(False)
    memgen = Nc.HaloCatalogMemberGenerator.new(hod)

    rng = Ncm.RNG.seeded_new(None, 3)
    members = memgen.generate(host, mset, rng)
    assert members is not None

    total = sum(len(members.find_children(p)) for p in (0, 1))
    assert total == members.len()
    # each host with a guaranteed central has at least one member.
    assert len(members.find_children(0)) >= 1
    assert len(members.find_children(1)) >= 1


def test_missing_column_raises() -> None:
    """A host lacking a required column raises a GLib.Error."""
    cosmo, mset = _cosmo_mset()
    bad = Nc.HaloCatalog.new(
        Nc.HaloCatalogKind.HALO, None, None, 1, ["ra", "dec", "z_true"], None
    )

    hod = Nc.GalaxyHODZheng07.new()
    memgen = Nc.HaloCatalogMemberGenerator.new(hod)
    rng = Ncm.RNG.seeded_new(None, 1)

    with pytest.raises(GLib.Error):
        memgen.generate(bad, mset, rng)


def test_distance_property_roundtrip() -> None:
    """A provided distance object is held and reused."""
    hod = Nc.GalaxyHODZheng07.new()
    memgen = Nc.HaloCatalogMemberGenerator.new(hod)
    assert memgen.peek_distance() is None

    dist = Nc.Distance.new(3.0)
    memgen.set_distance(dist)
    assert memgen.peek_distance() is dist

    memgen.set_distance(None)
    assert memgen.peek_distance() is None
