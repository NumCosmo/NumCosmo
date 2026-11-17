#
# wl_interp_cache.py
#
# Mon Nov 16 2026
# Copyright  2026  Caio Lima de Oliveira
# <caiolimadeoliveira@pm.me>
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

"""Caches of the shape interpolants of a weak-lensing catalog, kept apart from it.

The interpolants of ln P(eps_obs | g) that NcGalaxyShapeFactorSpectral builds for each
galaxy depend on the observed ellipticity and noise of the galaxy, on the population and on
the settings of the shape factor, but not on the cluster. They are kept in a cache file
beside the catalog, one per configuration, instead of inside it:

* the cache file is a NcGalaxyWLObs with the rows of the catalog and no columns, holding the
  interpolants of its rows, the record of what built them (shape-interp-source) and, in its
  metadata, the format and the fingerprint of the catalog;
* a JSON file beside it repeats what lookup needs (the fingerprint, the population and the
  settings), so that caches can be listed without reading them;
* :func:`attach_shape_interp_cache` copies the interpolants of a cache into a catalog after
  checking that the rows are those the cache was built for (the fingerprint) and that the
  record agrees with the fit (nc_galaxy_shape_factor_spectral_check_record(), which never
  aborts);
* :func:`update_shape_interp_cache` builds the interpolants a fit is missing, and those of
  a cache that do not cover the shears the galaxies reach under the priors of the fit, and
  writes the cache.
"""

from __future__ import annotations

import dataclasses
import hashlib
import json
import os
from pathlib import Path
from typing import Any, cast

import numpy as np

from numcosmo_py import Nc, Ncm

# Format of the cache file, kept in its metadata
SHAPE_INTERP_CACHE_FORMAT = "nc-shape-interp-cache-v1"

# Keys of the metadata of a catalog (ncm_catalog_peek_meta()) naming its cache; the first is
# also read by NcDataClusterWLFactor, which warns when it names a cache that is not attached.
CATALOG_META_INTERP_CACHE = "shape_interp_cache"
CATALOG_META_INTERP_FINGERPRINT = "shape_interp_cache_fingerprint"

# Keys of the metadata of a cache file
_CACHE_META_FORMAT = "format"
_CACHE_META_FINGERPRINT = "fingerprint"

# The columns of a catalog that the interpolant of a row depends on
_FINGERPRINT_COLUMNS = (
    "epsilon_obs_1",
    "epsilon_obs_2",
    "std_noise",
    "e_rms",
    "std_shape",
)


def wl_obs_fingerprint(obs: Nc.GalaxyWLObs) -> str:
    """Return the fingerprint of the rows of obs that its interpolants depend on.

    A SHA-256 of the number of rows, the ellipticity convention and frame, and the
    columns epsilon_obs_1, epsilon_obs_2, std_noise and, when present, e_rms (or its
    former name std_shape), as little-endian float64 in that order. It is the only link
    between row i of a catalog and interpolant i of a cache: the library does not compare
    them when it loads an interpolant.
    """
    h = hashlib.sha256()
    h.update(f"rows={obs.len()};".encode())
    h.update(f"conv={obs.get_ellip_conv().value_nick};".encode())
    h.update(f"coord={obs.get_coord().value_nick};".encode())
    columns = list(obs.peek_columns())
    data = obs.peek_data().to_numpy() if obs.len() > 0 and columns else None
    for name in _FINGERPRINT_COLUMNS:
        if name not in columns:
            continue
        assert data is not None
        h.update(f"{name}:".encode())
        h.update(
            np.ascontiguousarray(data[:, columns.index(name)], dtype="<f8").tobytes()
        )
    return h.hexdigest()


@dataclasses.dataclass
class ShapeInterpCache:
    """A cache of the shape interpolants of a catalog, as read from its file."""

    path: Path
    holder: Nc.GalaxyWLObs
    fingerprint: str

    @property
    def n_rows(self) -> int:
        """Return the number of rows of the catalog of the cache."""
        return self.holder.len()

    @property
    def source(self) -> Ncm.ObjDictStr | None:
        """Return the record of what built the interpolants."""
        return self.holder.peek_shape_interp_source()

    @property
    def n_modes(self) -> int:
        """Return the number of modes of an expansion over the population, 0 for none."""
        for i in range(self.holder.len()):
            modes = self.holder.peek_shape_interp_modes(i)
            if modes is not None:
                return modes.len()
        return 0


def _cache_info(holder: Nc.GalaxyWLObs, fingerprint: str) -> dict[str, Any]:
    """Return what lookup needs of a cache, for its JSON file."""
    source = holder.peek_shape_interp_source()
    assert source is not None
    gsfs = cast(Nc.GalaxyShapeFactorSpectral, source.peek("shape-factor"))
    pop = cast(Nc.GalaxyShapePop, source.peek("shape-pop"))
    params: dict[str, float] = {}
    coverage: dict[str, list[float]] = {}
    for i in range(pop.len()):
        name = pop.param_name(i)
        params[name] = pop.param_get(i)
        if pop.param_get_ftype(i) == Ncm.ParamType.FREE:
            coverage[name] = [
                pop.param_get_lower_bound(i),
                pop.param_get_upper_bound(i),
            ]
    n_interp = sum(
        1
        for i in range(holder.len())
        if holder.peek_shape_interp(i) is not None
        or holder.peek_shape_interp_modes(i) is not None
    )
    return {
        "format": SHAPE_INTERP_CACHE_FORMAT,
        "fingerprint": fingerprint,
        "rows": holder.len(),
        "rows_with_interpolants": n_interp,
        "ellip_conv": holder.get_ellip_conv().value_nick,
        "shape_factor": {
            "type": type(gsfs).__name__,
            "reltol": gsfs.get_reltol(),
            "interp_g_max": gsfs.get_interp_g_max(),
            "interp_abstol": gsfs.get_interp_abstol(),
            "interp_min_radial_panels": gsfs.get_interp_min_radial_panels(),
            "restrict_range": gsfs.get_restrict_range(),
            "pop_order": gsfs.get_pop_order(),
        },
        "population": {
            "type": type(pop).__name__,
            "params": params,
            "coverage": coverage,
        },
        "numcosmo_version": Ncm.cfg_get_version_string(),
    }


def _replace_atomically(path: Path, write) -> None:
    """Write path through write(part) on a temporary file, then rename it into place."""
    part = path.with_name(f".{path.name}.{os.getpid()}.part")
    try:
        write(part)
        os.replace(part, path)
    finally:
        if part.exists():
            part.unlink()


def write_shape_interp_cache(obs: Nc.GalaxyWLObs, path: Path) -> Path:
    """Write the interpolants of obs and their record to the cache file path.

    obs must carry current interpolants (a record). The JSON file of the cache is written
    beside it, with the suffix .json, after the cache file itself.

    Returns:
        The path of the JSON file.

    Raises:
        ValueError: obs carries no current interpolants.
    """
    source = obs.peek_shape_interp_source()
    if source is None:
        raise ValueError("the catalog carries no current interpolants to cache")

    holder = Nc.GalaxyWLObs.new(obs.get_ellip_conv(), obs.get_coord(), obs.len(), [])
    for i in range(obs.len()):
        modes = obs.peek_shape_interp_modes(i)
        if modes is not None:
            holder.set_shape_interp_modes(i, modes)
            continue
        spp = obs.peek_shape_interp(i)
        if spp is not None:
            holder.set_shape_interp(i, spp)
    holder.props.shape_interp_source = source

    fingerprint = wl_obs_fingerprint(obs)
    meta = Ncm.VarDict.new()
    meta.set_string(_CACHE_META_FORMAT, SHAPE_INTERP_CACHE_FORMAT)
    meta.set_string(_CACHE_META_FINGERPRINT, fingerprint)
    holder.set_meta(meta)

    path = Path(path)
    json_path = path.with_suffix(".json")
    info = _cache_info(holder, fingerprint)

    _replace_atomically(
        path,
        lambda part: Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(
            holder, part.as_posix()
        ),
    )
    _replace_atomically(
        json_path, lambda part: part.write_text(json.dumps(info, indent=2) + "\n")
    )
    return json_path


def read_shape_interp_cache(path: Path) -> ShapeInterpCache:
    """Read a cache file written by :func:`write_shape_interp_cache`.

    Raises:
        ValueError: The file is not a cache of shape interpolants of a known format.
    """
    holder = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).from_binfile(
        Path(path).as_posix()
    )
    if not isinstance(holder, Nc.GalaxyWLObs):
        raise ValueError(f"{path} is not a cache of shape interpolants")  # noqa: TRY004
    meta = holder.peek_meta()
    ok_format, fmt = (
        (False, None) if meta is None else meta.get_string(_CACHE_META_FORMAT)
    )
    ok_fp, fingerprint = (
        (False, None) if meta is None else meta.get_string(_CACHE_META_FINGERPRINT)
    )
    if not ok_format or fmt != SHAPE_INTERP_CACHE_FORMAT or not ok_fp:
        raise ValueError(
            f"{path} is not a cache of shape interpolants of format "
            f"{SHAPE_INTERP_CACHE_FORMAT} (format: {fmt})"
        )
    return ShapeInterpCache(
        path=Path(path), holder=holder, fingerprint=str(fingerprint)
    )


def read_shape_interp_cache_info(path: Path) -> dict[str, Any]:
    """Return the JSON description of the cache file path."""
    return json.loads(Path(path).with_suffix(".json").read_text())


def check_shape_interp_cache(
    cache: ShapeInterpCache,
    obs: Nc.GalaxyWLObs,
    gsfs: Nc.GalaxyShapeFactorSpectral,
    pop: Nc.GalaxyShapePop,
) -> str | None:
    """Return why cache cannot be attached to obs for gsfs and pop, or None when it can.

    The rows of obs must be those the cache was built for (the fingerprint), and the
    record of the cache must agree with gsfs and pop as loading the catalog would require
    (nc_galaxy_shape_factor_spectral_check_record()).
    """
    if cache.n_rows != obs.len():
        return f"the cache has {cache.n_rows} rows, the catalog {obs.len()}"
    if cache.fingerprint != wl_obs_fingerprint(obs):
        return "the cache was built for other rows (its fingerprint differs)"
    ok, diff = gsfs.check_record(pop, cache.source, cache.n_modes)
    if not ok:
        return f"the cache was built for other settings:{diff}"
    return None


def attach_shape_interp_cache(
    cache: ShapeInterpCache,
    obs: Nc.GalaxyWLObs,
    gsfs: Nc.GalaxyShapeFactorSpectral,
    pop: Nc.GalaxyShapePop,
) -> None:
    """Copy the interpolants of cache into obs, and name the cache in its metadata.

    Raises:
        ValueError: The cache cannot be attached (:func:`check_shape_interp_cache`).
    """
    reason = check_shape_interp_cache(cache, obs, gsfs, pop)
    if reason is not None:
        raise ValueError(f"cannot attach {cache.path}: {reason}")

    holder = cache.holder
    for i in range(holder.len()):
        modes = holder.peek_shape_interp_modes(i)
        if modes is not None:
            obs.set_shape_interp_modes(i, modes)
            continue
        spp = holder.peek_shape_interp(i)
        if spp is not None:
            obs.set_shape_interp(i, spp)
    obs.props.shape_interp_source = cache.source
    name_shape_interp_cache(obs, cache.path, cache.fingerprint)


def name_shape_interp_cache(obs: Nc.GalaxyWLObs, path: Path, fingerprint: str) -> None:
    """Name the cache of obs in its metadata (the path and the fingerprint)."""
    meta = obs.peek_meta() or Ncm.VarDict.new()
    meta.set_string(CATALOG_META_INTERP_CACHE, Path(path).as_posix())
    meta.set_string(CATALOG_META_INTERP_FINGERPRINT, fingerprint)
    obs.set_meta(meta)


def strip_shape_interp(obs: Nc.GalaxyWLObs) -> None:
    """Remove the interpolants and their record from obs, keeping its rows and metadata."""
    obs.props.shape_interp = Ncm.ObjDictInt.new()
    obs.props.shape_interp_modes = Ncm.ObjDictInt.new()
    obs.props.shape_interp_source = Ncm.ObjDictStr.new()


def find_shape_interp_cache(
    directory: Path,
    obs: Nc.GalaxyWLObs,
    gsfs: Nc.GalaxyShapeFactorSpectral,
    pop: Nc.GalaxyShapePop,
) -> ShapeInterpCache | None:
    """Return the first cache in directory that can be attached to obs, or None.

    The JSON files select the candidates (format, fingerprint, population type and
    pop-order); each candidate is then read and checked by
    :func:`check_shape_interp_cache`.
    """
    directory = Path(directory)
    if not directory.is_dir():
        return None
    fingerprint = wl_obs_fingerprint(obs)
    for json_path in sorted(directory.glob("*.json")):
        try:
            info = json.loads(json_path.read_text())
        except (OSError, ValueError):
            continue
        if (
            info.get("format") != SHAPE_INTERP_CACHE_FORMAT
            or info.get("fingerprint") != fingerprint
            or info.get("population", {}).get("type") != type(pop).__name__
            or info.get("shape_factor", {}).get("pop_order") != gsfs.get_pop_order()
        ):
            continue
        gvar = json_path.with_suffix(".gvar")
        if not gvar.exists():
            continue
        cache = read_shape_interp_cache(gvar)
        if check_shape_interp_cache(cache, obs, gsfs, pop) is None:
            return cache
    return None


@dataclasses.dataclass(frozen=True)
class ShapeInterpCacheUpdate:
    """What :func:`update_shape_interp_cache` did."""

    rows: int
    built: int
    reused: int
    path: Path


def _builder_mset(
    mset: Ncm.MSet, coverage: dict[str, tuple[float, float]] | None
) -> Ncm.MSet:
    """Return mset, or a copy whose population parameters have the bounds of coverage."""
    if not coverage:
        return mset
    builder = mset.dup(Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP))
    pop = builder.peek(Nc.GalaxyShapePop.id())
    assert isinstance(pop, Nc.GalaxyShapePop)
    for name, (lo, hi) in coverage.items():
        found, pid = pop.param_index_from_name(name)
        if not found:
            raise ValueError(f"the population has no parameter {name}")
        if not lo <= pop.param_get(pid) <= hi:
            pop.param_set(pid, 0.5 * (lo + hi))
        pop.param_set_lower_bound(pid, lo)
        pop.param_set_upper_bound(pid, hi)
        pop.param_set_ftype(pid, Ncm.ParamType.FREE)
    builder.prepare_fparam_map()
    return builder


def update_shape_interp_cache(
    data: Nc.DataClusterWLFactor,
    mset: Ncm.MSet,
    path: Path,
    *,
    coverage: dict[str, tuple[float, float]] | None = None,
) -> ShapeInterpCacheUpdate:
    """Bring the cache file path up to date for the fit of data and mset, and attach it.

    The interpolants a galaxy reaches under the priors of mset (restrict-range) are those
    the cache must hold. The cache at path, when it exists, can be attached and was built
    over coverage, is attached; then the galaxies without an interpolant, and those whose
    interpolant does not cover the shears they reach under the priors of mset, are built,
    in parallel, and the cache is written again. coverage gives the bounds of the
    population parameters the expansion spans (pop-order), when they are not the priors of
    mset. Afterwards obs carries the interpolants and names the cache, and data reads them
    again, so that the fit loads them.

    Raises:
        ValueError: The shape factor of data does not use interpolants.
    """
    obs = data.peek_obs()
    gsfs = data.props.shape_factor
    if not isinstance(gsfs, Nc.GalaxyShapeFactorSpectral) or not gsfs.get_use_interp():
        raise ValueError(
            "the shape factor must be NcGalaxyShapeFactorSpectral with use-interp"
        )

    builder = _builder_mset(mset, coverage)
    pop = builder.peek(Nc.GalaxyShapePop.id())
    assert isinstance(pop, Nc.GalaxyShapePop)

    strip_shape_interp(obs)
    path = Path(path)
    if path.exists():
        cache = read_shape_interp_cache(path)
        if check_shape_interp_cache(cache, obs, gsfs, pop) is None and _same_coverage(
            cache, pop
        ):
            attach_shape_interp_cache(cache, obs, gsfs, pop)
    data.set_obs(obs)

    n0 = gsfs.get_n_built()
    gsfs.set_rebuild_short(True)
    try:
        data.write_shape_interp(builder, obs)
    finally:
        gsfs.set_rebuild_short(False)
    built = gsfs.get_n_built() - n0

    write_shape_interp_cache(obs, path)
    name_shape_interp_cache(obs, path, wl_obs_fingerprint(obs))
    # The fit reads the interpolants again as loaded ones, for its own priors
    data.set_obs(obs)
    return ShapeInterpCacheUpdate(
        rows=obs.len(), built=built, reused=obs.len() - built, path=path
    )


def _same_coverage(cache: ShapeInterpCache, pop: Nc.GalaxyShapePop) -> bool:
    """Return whether the expansion of cache spans the free parameters of pop over their bounds."""
    if cache.n_modes == 0:
        return True
    source = cache.source
    assert source is not None
    src = cast(Nc.GalaxyShapePop, source.peek("shape-pop"))
    for i in range(pop.len()):
        free = pop.param_get_ftype(i) == Ncm.ParamType.FREE
        if free != (src.param_get_ftype(i) == Ncm.ParamType.FREE):
            return False
        if free and (
            pop.param_get_lower_bound(i) != src.param_get_lower_bound(i)
            or pop.param_get_upper_bound(i) != src.param_get_upper_bound(i)
        ):
            return False
    return True
