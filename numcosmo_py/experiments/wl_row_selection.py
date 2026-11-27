#
# wl_row_selection.py
#
# Fri Nov 27 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# wl_row_selection.py
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

"""Row selections of a weak-lensing catalog, read from YAML.

A cluster weak-lensing likelihood uses the rows of its catalog that its mask keeps
(``NcDataClusterWLFactor:mask``). A row selection describes that mask by criteria
instead of a list of booleans: leaf selectors test one property of each row, and the
combinators ``and``, ``or`` and ``not`` compose them. Each selector is a mapping with a
``kind`` key::

    kind: and
    selectors:
      - kind: radius
        r_min: 0.446
        r_max: 4.46
      - kind: column
        column: z
        max: 2.5
      - kind: z-support
        dz: 0.2
        p_min: 0.98

The selection is evaluated once, when the experiment is built, at the cluster center,
redshift and cosmology of the experiment; the likelihood then uses the resulting fixed
mask, so the selected rows do not follow a fitted center.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Annotated, Literal

import numpy as np
import yaml
from pydantic import BaseModel, ConfigDict, Field, TypeAdapter, model_validator

from numcosmo_py import Nc, Ncm
from numcosmo_py.helper import catalog_column


@dataclass(frozen=True)
class RowSelectionContext:
    """The catalog and the cluster at which a row selection is evaluated."""

    obs: Nc.GalaxyWLObs
    cosmo: Nc.HICosmo
    halo_position: Nc.HaloPosition


class _RowSelectorBase(BaseModel):
    """Common configuration of the row selectors."""

    model_config = ConfigDict(extra="forbid", frozen=True)


class RadiusRowSelector(_RowSelectorBase):
    """Rows whose projected radius from the cluster center is in [r_min, r_max] Mpc."""

    kind: Literal["radius"] = "radius"
    r_min: Annotated[float, Field(ge=0.0)] = 0.0
    r_max: Annotated[float | None, Field(gt=0.0)] = None

    @model_validator(mode="after")
    def _check_range(self) -> "RadiusRowSelector":
        if (self.r_max is not None) and not self.r_min < self.r_max:
            raise ValueError(
                f"radius: r_min = {self.r_min} is not below r_max = {self.r_max}"
            )
        return self

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        r = np.array(
            ctx.halo_position.projected_radius_array(
                ctx.cosmo,
                catalog_column(ctx.obs, "ra").tolist(),
                catalog_column(ctx.obs, "dec").tolist(),
            )
        )
        keep = r >= self.r_min
        if self.r_max is not None:
            keep &= r <= self.r_max
        return keep


class ColumnRowSelector(_RowSelectorBase):
    """Rows whose value of a numeric column is in [min, max]; a missing bound is open."""

    kind: Literal["column"] = "column"
    column: str
    min: float | None = None
    max: float | None = None

    @model_validator(mode="after")
    def _check_bounds(self) -> "ColumnRowSelector":
        if (self.min is None) and (self.max is None):
            raise ValueError(f"column {self.column}: give min, max or both")
        if (self.min is not None) and (self.max is not None) and self.min > self.max:
            raise ValueError(
                f"column {self.column}: min = {self.min} is above max = {self.max}"
            )
        return self

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        if not ctx.obs.has_column(self.column):
            raise ValueError(f"column: the catalog has no column {self.column}")
        values = catalog_column(ctx.obs, self.column)
        keep = np.ones(len(values), dtype=bool)
        if self.min is not None:
            keep &= values >= self.min
        if self.max is not None:
            keep &= values <= self.max
        return keep


class ZSupportRowSelector(_RowSelectorBase):
    """Rows whose redshift density puts at least p_min above z_cl + dz.

    The probability is P(z > z_cl + dz) = 1 - CDF(z_cl + dz) of the row's NcGalaxyPz,
    with z_cl the redshift of the cluster.
    """

    kind: Literal["z-support"] = "z-support"
    dz: float = 0.0
    p_min: Annotated[float, Field(ge=0.0, le=1.0)]

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        z_cl = ctx.halo_position["z"]
        cdf = np.array(ctx.obs.pz_cdf_array(z_cl + self.dz))
        return 1.0 - cdf >= self.p_min


class ObjectIdRowSelector(_RowSelectorBase):
    """Rows whose identifier, an INT column, is in a list."""

    kind: Literal["object-id"] = "object-id"
    ids: list[int]
    column: str = "object_id"

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        if not ctx.obs.has_column(self.column):
            raise ValueError(f"object-id: the catalog has no column {self.column}")
        if ctx.obs.get_col_type(self.column) != Ncm.CatalogColType.INT:
            raise ValueError(
                f"object-id: the column {self.column} is not an INT column"
            )
        return np.isin(
            catalog_column(ctx.obs, self.column), np.array(self.ids, dtype=np.int64)
        )


class AndRowSelector(_RowSelectorBase):
    """Rows kept by every selector of the list."""

    kind: Literal["and"] = "and"
    selectors: Annotated[list["RowSelector"], Field(min_length=1)]

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        keep = self.selectors[0].keep(ctx)
        for selector in self.selectors[1:]:
            keep &= selector.keep(ctx)
        return keep


class OrRowSelector(_RowSelectorBase):
    """Rows kept by at least one selector of the list."""

    kind: Literal["or"] = "or"
    selectors: Annotated[list["RowSelector"], Field(min_length=1)]

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        keep = self.selectors[0].keep(ctx)
        for selector in self.selectors[1:]:
            keep |= selector.keep(ctx)
        return keep


class NotRowSelector(_RowSelectorBase):
    """Rows not kept by the selector."""

    kind: Literal["not"] = "not"
    selector: "RowSelector"

    def keep(self, ctx: RowSelectionContext) -> np.ndarray:
        """Return whether each row is kept."""
        return ~self.selector.keep(ctx)


RowSelector = Annotated[
    RadiusRowSelector
    | ColumnRowSelector
    | ZSupportRowSelector
    | ObjectIdRowSelector
    | AndRowSelector
    | OrRowSelector
    | NotRowSelector,
    Field(discriminator="kind"),
]

for _model in (AndRowSelector, OrRowSelector, NotRowSelector):
    _model.model_rebuild()

_ROW_SELECTOR_ADAPTER: TypeAdapter[RowSelector] = TypeAdapter(RowSelector)


def load_row_selector(source: str | Path) -> RowSelector:
    """Read a row selection from a YAML file or from inline YAML.

    :param source: the name of an existing file, read as YAML, or the YAML text itself,
        such as ``"{kind: radius, r_max: 4.46}"``.
    :return: the row selector.
    :raises ValueError: the source is not a valid row selection; a value that is neither
        an existing file nor a YAML mapping is reported as such.
    """
    path = Path(source)
    text = path.read_text(encoding="utf-8") if path.is_file() else str(source)
    data = yaml.safe_load(text)
    if isinstance(data, dict):
        return _ROW_SELECTOR_ADAPTER.validate_python(data)
    raise ValueError(
        f"{source} is neither an existing file nor a YAML row selection (a mapping)"
    )


def row_selection_mask(selector: RowSelector, ctx: RowSelectionContext) -> list[bool]:
    """Evaluate a row selection: whether the likelihood uses each row of the catalog.

    :param selector: the row selector.
    :param ctx: the catalog and the cluster.
    :return: one boolean per row.
    :raises ValueError: the selection keeps no row.
    """
    keep = selector.keep(ctx)
    if not keep.any():
        raise ValueError("the row selection keeps no row of the catalog")
    return keep.tolist()
