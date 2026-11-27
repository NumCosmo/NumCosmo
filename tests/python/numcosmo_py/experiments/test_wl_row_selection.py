#
# test_wl_row_selection.py
#
# Fri Nov 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_wl_row_selection.py
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

"""Tests of the row selections of a weak-lensing catalog."""

from pathlib import Path

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.experiments.wl_row_selection import (
    AndRowSelector,
    ColumnRowSelector,
    NotRowSelector,
    ObjectIdRowSelector,
    OrRowSelector,
    RadiusRowSelector,
    RowSelectionContext,
    ZSupportRowSelector,
    load_row_selector,
    row_selection_mask,
)

pytestmark = pytest.mark.app

Ncm.cfg_init()

N_ROWS = 30
Z_CL = 0.3


@pytest.fixture(name="ctx")
def fixture_ctx() -> RowSelectionContext:
    """A catalog of N_ROWS galaxies around a cluster at (0, 0, Z_CL), each with a p(z)."""
    rng = np.random.default_rng(20261127)
    cols = ["ra", "dec", "z", "object_id"]
    types = [Ncm.CatalogColType.DOUBLE] * 3 + [Ncm.CatalogColType.INT]
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE,
        Nc.WLEllipticityFrame.CELESTIAL,
        N_ROWS,
        cols,
        types,
    )
    z = rng.uniform(0.1, 1.5, N_ROWS)
    for i in range(N_ROWS):
        obs.set("ra", i, float(rng.uniform(-0.2, 0.2)))
        obs.set("dec", i, float(rng.uniform(-0.2, 0.2)))
        obs.set("z", i, float(z[i]))
        obs.set_int("object_id", i, 2**60 + 7 * i)
        zs = np.linspace(max(z[i] - 0.3, 0.01), z[i] + 0.3, 16)
        obs.set_pz(
            i,
            Nc.GalaxyPzBSplineSqrt.new_from_table(
                Ncm.Vector.new_array(zs.tolist()),
                Ncm.Vector.new_array(np.exp(-0.5 * ((zs - z[i]) / 0.1) ** 2).tolist()),
            ),
        )

    cosmo = Nc.HICosmoDEXcdm()
    hp = Nc.HaloPosition.new(Nc.Distance.new(6.0))
    hp["ra"], hp["dec"], hp["z"] = 0.0, 0.0, Z_CL
    hp.prepare(cosmo)

    return RowSelectionContext(obs=obs, cosmo=cosmo, halo_position=hp)


def _radii(ctx: RowSelectionContext) -> np.ndarray:
    """The projected radius of each row, one row at a time."""
    obs = ctx.obs
    return np.array(
        [
            ctx.halo_position.projected_radius_from_ra_dec(
                ctx.cosmo, obs.get("ra", i), obs.get("dec", i)
            )
            for i in range(obs.len())
        ]
    )


def test_radius(ctx: RowSelectionContext) -> None:
    """The radius selector keeps the rows with r_min <= r <= r_max."""
    r = _radii(ctx)
    r_min, r_max = np.quantile(r, [0.25, 0.75])
    keep = RadiusRowSelector(r_min=r_min, r_max=r_max).keep(ctx)

    np.testing.assert_array_equal(keep, (r >= r_min) & (r <= r_max))
    np.testing.assert_array_equal(RadiusRowSelector(r_min=r_min).keep(ctx), r >= r_min)


def test_column(ctx: RowSelectionContext) -> None:
    """The column selector keeps the rows inside its bounds; a missing bound is open."""
    z = np.array([ctx.obs.get("z", i) for i in range(N_ROWS)])

    np.testing.assert_array_equal(
        ColumnRowSelector(column="z", min=0.5, max=1.0).keep(ctx),
        (z >= 0.5) & (z <= 1.0),
    )
    np.testing.assert_array_equal(
        ColumnRowSelector(column="z", max=1.0).keep(ctx), z <= 1.0
    )


def test_z_support(ctx: RowSelectionContext) -> None:
    """The z-support selector keeps the rows with P(z > z_cl + dz) >= p_min."""
    p_above = np.array(
        [1.0 - ctx.obs.peek_pz(i).cdf(Z_CL + 0.2) for i in range(N_ROWS)]
    )
    keep = ZSupportRowSelector(dz=0.2, p_min=0.9).keep(ctx)

    np.testing.assert_array_equal(keep, p_above >= 0.9)
    assert 0 < keep.sum() < N_ROWS


def test_object_id(ctx: RowSelectionContext) -> None:
    """The object-id selector keeps the listed identifiers, exact above 2^53."""
    ids = [2**60 + 7 * 3, 2**60 + 7 * 11, 12345]
    keep = ObjectIdRowSelector(ids=ids).keep(ctx)

    assert np.flatnonzero(keep).tolist() == [3, 11]


def test_combinators(ctx: RowSelectionContext) -> None:
    """The and, or and not selectors combine the rows of their selectors."""
    near = RadiusRowSelector(r_max=float(np.median(_radii(ctx))))
    low_z = ColumnRowSelector(column="z", max=0.8)
    a, b = near.keep(ctx), low_z.keep(ctx)

    np.testing.assert_array_equal(
        AndRowSelector(selectors=[near, low_z]).keep(ctx), a & b
    )
    np.testing.assert_array_equal(
        OrRowSelector(selectors=[near, low_z]).keep(ctx), a | b
    )
    np.testing.assert_array_equal(NotRowSelector(selector=near).keep(ctx), ~a)


def test_load_yaml(ctx: RowSelectionContext, tmp_path: Path) -> None:
    """A YAML row selection reads into the selectors it describes."""
    path = tmp_path / "selection.yaml"
    path.write_text("""
kind: and
selectors:
  - kind: radius
    r_min: 0.1
  - kind: column
    column: z
    max: 2.5
  - kind: not
    selector:
      kind: object-id
      ids: [1152921504606846997]
""")
    selector = load_row_selector(path)

    assert selector == AndRowSelector(
        selectors=[
            RadiusRowSelector(r_min=0.1),
            ColumnRowSelector(column="z", max=2.5),
            NotRowSelector(selector=ObjectIdRowSelector(ids=[2**60 + 21])),
        ]
    )
    assert len(row_selection_mask(selector, ctx)) == N_ROWS


def test_load_inline(ctx: RowSelectionContext, tmp_path: Path) -> None:
    """Inline YAML reads as a file with the same text does."""
    path = tmp_path / "selection.yaml"
    path.write_text("kind: column\ncolumn: z\nmax: 1.0\n")
    inline = load_row_selector("{kind: column, column: z, max: 1.0}")

    assert inline == load_row_selector(path)
    assert inline == ColumnRowSelector(column="z", max=1.0)


def test_load_neither(tmp_path: Path) -> None:
    """A value that is neither an existing file nor a YAML mapping raises ValueError."""
    with pytest.raises(ValueError, match="neither an existing file"):
        load_row_selector(tmp_path / "missing.yaml")


@pytest.mark.parametrize(
    "text, match",
    [
        ("kind: box\n", "box"),
        ("kind: radius\nr_max: 1.0\nextra: 2\n", "extra"),
        ("kind: radius\nr_min: 2.0\nr_max: 1.0\n", "is not below"),
        ("kind: column\ncolumn: z\n", "give min, max or both"),
        ("kind: z-support\ndz: 0.2\np_min: 1.5\n", "p_min"),
        ("kind: and\nselectors: []\n", "at least 1"),
    ],
)
def test_load_yaml_invalid(tmp_path: Path, text: str, match: str) -> None:
    """An invalid row selection raises ValueError naming the problem."""
    path = tmp_path / "selection.yaml"
    path.write_text(text)

    with pytest.raises(ValueError, match=match):
        load_row_selector(path)


def test_errors_at_evaluation(ctx: RowSelectionContext) -> None:
    """A missing column, a non-INT identifier column and an empty selection raise."""
    with pytest.raises(ValueError, match="no column w"):
        ColumnRowSelector(column="w", max=1.0).keep(ctx)
    with pytest.raises(ValueError, match="not an INT column"):
        ObjectIdRowSelector(ids=[1], column="z").keep(ctx)
    with pytest.raises(ValueError, match="keeps no row"):
        row_selection_mask(ColumnRowSelector(column="z", min=10.0), ctx)
