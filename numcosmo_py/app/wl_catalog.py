#
# wl_catalog.py
#
# Fri Nov 13 2026
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
# with this program. If not, see <http://www.gnu.org/licenses/>.

"""NumCosmo APP commands that describe weak lensing catalogs and their caches."""

from __future__ import annotations

import dataclasses
import json
from pathlib import Path
from typing import Annotated

import typer
from rich.console import Console

from numcosmo_py import Nc, Ncm
from numcosmo_py.experiments.wl_interp_cache import (
    format_coverage,
    read_shape_interp_cache,
    recorded_coverage,
    resolve_shape_interp_cache,
    shape_interp_cache_info,
    wl_obs_fingerprint,
)


def _read_wl_obs_file(path: Path) -> tuple[Nc.GalaxyWLObs | None, str]:
    """Return the NcGalaxyWLObs of path and what the file is: catalog, cache or dataset."""
    obj = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED).from_binfile(
        path.absolute().as_posix()
    )
    if isinstance(obj, Ncm.Dataset):
        for i in range(obj.get_length()):
            data = obj.get_data(i)
            if isinstance(data, Nc.DataClusterWLFactor):
                return data.peek_obs(), "dataset"
        return None, "dataset"
    if isinstance(obj, Nc.GalaxyWLObs):
        meta = obj.peek_meta()
        if meta is not None and meta.get_string("format")[0]:
            return obj, "cache"
        return obj, "catalog"
    return None, type(obj).__name__


def _dataset_mask(path: Path) -> list[bool]:
    """Return the mask of the NcDataClusterWLFactor of the dataset in path, [] for none."""
    dset = Ncm.Serialize.new(Ncm.SerializeOpt.KEEP_SHARED).from_binfile(
        path.absolute().as_posix()
    )
    assert isinstance(dset, Ncm.Dataset)
    for i in range(dset.get_length()):
        data = dset.get_data(i)
        if isinstance(data, Nc.DataClusterWLFactor):
            return list(data.get_mask() or [])
    return []


def _interp_summary(obs: Nc.GalaxyWLObs) -> tuple[int, int]:
    """Return the number of rows of obs with a fixed interpolant and with modes."""
    lens = [m.len() for m in map(obs.peek_shape_interp, range(obs.len())) if m]
    n_fixed = sum(1 for n in lens if n == 1)
    n_modes = sum(1 for n in lens if n > 1)
    return n_fixed, n_modes


@dataclasses.dataclass(kw_only=True)
class InfoWLCatalog:
    """Describe a WL catalog, a cache of its shape interpolants, or an experiment dataset.

    For a catalog (or the catalog of a .dataset.gvar): its rows, columns, ellipticity
    convention, cluster metadata and the cache it names, whether that cache exists and
    matches.
    For a cache: the settings, population, coverage and order that built it.
    """

    file: Annotated[
        Path,
        typer.Argument(
            help="A catalog, a cache of interpolants or a .dataset.gvar.",
            exists=True,
            readable=True,
        ),
    ]

    def __post_init__(self) -> None:
        """Print the description."""
        Ncm.cfg_init()
        console = Console()
        obs, kind = _read_wl_obs_file(self.file)
        if obs is None:
            raise typer.BadParameter(f"{self.file} holds no WL catalog ({kind})")

        size = self.file.stat().st_size / 2**20
        console.print(f"{self.file}: {kind}, {size:.1f} MiB, {obs.len()} rows")
        if kind == "dataset":
            mask = _dataset_mask(self.file)
            if mask:
                console.print(f"  the likelihood uses {sum(mask)} of the rows (mask)")
        if kind == "cache":
            cache = read_shape_interp_cache(self.file)
            info = shape_interp_cache_info(cache.holder, cache.fingerprint)
            info["numcosmo_version_reading"] = info.pop("numcosmo_version")
            console.print_json(json.dumps(info))
            return

        console.print(
            f"  convention {obs.get_ellip_conv().value_nick}, frame "
            f"{obs.get_coord().value_nick}"
        )
        console.print(f"  columns: {', '.join(obs.peek_columns())}")
        meta = obs.peek_meta()
        if meta is not None:
            for key in ("cluster_ra", "cluster_dec", "cluster_z", "cluster_c"):
                ok, value = meta.get_double(key)
                if ok:
                    console.print(f"  {key} = {value}")

        path = resolve_shape_interp_cache(obs, self.file.absolute().parent)
        if path is None:
            console.print("  names no cache of interpolants")
            return
        coverage = recorded_coverage(obs)
        console.print(
            f"  names the cache {path}"
            + (f", coverage {format_coverage(coverage)}" if coverage else "")
        )
        if not path.exists():
            console.print("  the cache does not exist yet: running builds it")
            return
        cache = read_shape_interp_cache(path)
        ok = cache.fingerprint == wl_obs_fingerprint(obs)
        n_fixed, n_modes = _interp_summary(cache.holder)
        console.print(
            f"  the cache exists: {n_fixed} fixed, {n_modes} expansions, "
            + ("matches the rows" if ok else "built for OTHER rows (fingerprint)")
        )
