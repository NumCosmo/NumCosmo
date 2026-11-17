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
* its name hashes the fingerprint, the settings, the population and its coverage
  (:func:`shape_interp_cache_name`), so finding the cache of a configuration is checking that
  the file exists;
* :func:`attach_shape_interp_cache` copies the interpolants of a cache into a catalog after
  checking that the rows are those the cache was built for (the fingerprint) and that the
  record agrees with the fit (nc_galaxy_shape_factor_spectral_check_record(), which never
  aborts);
* :func:`update_shape_interp_cache` builds the interpolants a fit is missing, and those of
  a cache that do not cover the shears the galaxies reach under the priors of the fit, and
  writes the cache.
"""

from __future__ import annotations

import contextlib
import dataclasses
import fcntl
import hashlib
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
CATALOG_META_INTERP_COVERAGE = "shape_interp_cache_coverage"

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
        """Return the number of modes of each galaxy, 1 for a fixed population, 0 for none."""
        for i in range(self.holder.len()):
            modes = self.holder.peek_shape_interp(i)
            if modes is not None:
                return modes.len()
        return 0


def shape_interp_cache_info(holder: Nc.GalaxyWLObs, fingerprint: str) -> dict[str, Any]:
    """Return the description of a cache: its format, rows, settings and population."""
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
        1 for i in range(holder.len()) if holder.peek_shape_interp(i) is not None
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

    obs must carry current interpolants (a record).

    Returns:
        The path of the cache file.

    Raises:
        ValueError: obs carries no current interpolants.
    """
    source = obs.peek_shape_interp_source()
    if source is None:
        raise ValueError("the catalog carries no current interpolants to cache")

    holder = Nc.GalaxyWLObs.new(obs.get_ellip_conv(), obs.get_coord(), obs.len(), [])
    for i in range(obs.len()):
        modes = obs.peek_shape_interp(i)
        if modes is not None:
            holder.set_shape_interp(i, modes)
    holder.props.shape_interp_source = source

    fingerprint = wl_obs_fingerprint(obs)
    meta = Ncm.VarDict.new()
    meta.set_string(_CACHE_META_FORMAT, SHAPE_INTERP_CACHE_FORMAT)
    meta.set_string(_CACHE_META_FINGERPRINT, fingerprint)
    holder.set_meta(meta)

    path = Path(path)
    _replace_atomically(
        path,
        lambda part: Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(
            holder, part.as_posix()
        ),
    )
    return path


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
    source = cache.source
    if source is None:
        return "the cache holds no record of what built it"
    ok, diff = gsfs.check_record(pop, source, cache.n_modes)
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
        modes = holder.peek_shape_interp(i)
        if modes is not None:
            obs.set_shape_interp(i, modes)
    source = cache.source
    assert source is not None
    obs.props.shape_interp_source = source
    name_shape_interp_cache(obs, cache.path, cache.fingerprint)


def name_shape_interp_cache(
    obs: Nc.GalaxyWLObs,
    path: Path,
    fingerprint: str,
    *,
    base_dir: Path | None = None,
) -> None:
    """Name the cache of obs in its metadata: its path and the fingerprint of obs.

    With base_dir the path is kept relative to it (the directory of the file that will
    hold obs), so that the files can be moved together; :func:`resolve_shape_interp_cache`
    resolves it against the same directory.
    """
    path = Path(path)
    if base_dir is not None:
        path = Path(os.path.relpath(path.absolute(), Path(base_dir).absolute()))
    meta = obs.peek_meta() or Ncm.VarDict.new()
    meta.set_string(CATALOG_META_INTERP_CACHE, path.as_posix())
    meta.set_string(CATALOG_META_INTERP_FINGERPRINT, fingerprint)
    obs.set_meta(meta)


def resolve_shape_interp_cache(obs: Nc.GalaxyWLObs, base_dir: Path) -> Path | None:
    """Return the path of the cache obs names, resolved against base_dir, or None."""
    meta = obs.peek_meta()
    if meta is None:
        return None
    ok, name = meta.get_string(CATALOG_META_INTERP_CACHE)
    if not ok:
        return None
    path = Path(name)
    return path if path.is_absolute() else Path(base_dir) / path


def strip_shape_interp(obs: Nc.GalaxyWLObs) -> None:
    """Remove the interpolants and their record from obs, keeping its rows and metadata."""
    obs.props.shape_interp = Ncm.ObjDictInt.new()
    obs.props.shape_interp_source = Ncm.ObjDictStr.new()


@dataclasses.dataclass(frozen=True)
class ShapeInterpCacheUpdate:
    """What :func:`update_shape_interp_cache` did."""

    rows: int
    built: int
    reused: int
    unconverged: int
    path: Path


# Relative tolerance of the comparison of bounds, as the library compares a record
_COVERAGE_RELTOL = 1.0e-15


def _close(a: float, b: float) -> bool:
    """Return whether a and b agree within _COVERAGE_RELTOL relative."""
    return abs(a - b) <= _COVERAGE_RELTOL * max(abs(a), abs(b))


def _builder_mset(
    mset: Ncm.MSet, coverage: dict[str, tuple[float, float]] | None
) -> Ncm.MSet:
    """Return a copy of mset whose population parameters of coverage span it.

    The parameters of coverage are free over its bounds; their values are kept when
    inside, and set to the centre otherwise. The other parameters of mset, the priors of
    the fit, are unchanged: they decide the shears each galaxy reaches.
    """
    builder = mset.dup(Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP))
    pop = builder.peek(Nc.GalaxyShapePop.id())
    assert isinstance(pop, Nc.GalaxyShapePop)
    for pid in range(pop.len()):
        if pop.param_get_ftype(pid) == Ncm.ParamType.FREE and pop.param_name(
            pid
        ) not in (coverage or {}):
            pop.param_set_ftype(pid, Ncm.ParamType.FIXED)
    for name, (lo, hi) in (coverage or {}).items():
        found, pid = pop.param_index_from_name(name)
        if not found:
            raise ValueError(f"the population has no parameter {name}")
        if not lo < hi:
            raise ValueError(f"the coverage of {name} needs LO < HI: [{lo}, {hi}]")
        if not lo <= pop.param_get(pid) <= hi:
            pop.param_set(pid, 0.5 * (lo + hi))
        if lo >= pop.param_get_upper_bound(pid):
            pop.param_set_upper_bound(pid, hi)
            pop.param_set_lower_bound(pid, lo)
        else:
            pop.param_set_lower_bound(pid, lo)
            pop.param_set_upper_bound(pid, hi)
        pop.param_set_ftype(pid, Ncm.ParamType.FREE)
    builder.prepare_fparam_map()
    return builder


def shape_interp_cache_name(
    obs: Nc.GalaxyWLObs,
    gsfs: Nc.GalaxyShapeFactorSpectral,
    pop: Nc.GalaxyShapePop,
    coverage: dict[str, tuple[float, float]] | None,
) -> str:
    """Return the file name of the cache of obs for this configuration.

    A hash of everything the interpolants depend on: the rows of obs (their fingerprint),
    the settings of gsfs the record compares, restrict-range and pop-order, the type of the
    population, the values of its parameters outside coverage and the bounds of those in
    it. Two configurations never share a file.
    """
    h = hashlib.sha256()
    h.update(wl_obs_fingerprint(obs).encode())
    h.update(
        (
            f"conv={gsfs.get_ellip_conv().value_nick};reltol={gsfs.get_reltol()!r};"
            f"g_max={gsfs.get_interp_g_max()!r};abstol={gsfs.get_interp_abstol()!r};"
            f"panels={gsfs.get_interp_min_radial_panels()};"
            f"restrict={gsfs.get_restrict_range()};order={gsfs.get_pop_order()};"
            f"pop={type(pop).__name__};"
        ).encode()
    )
    cov = coverage or {}
    for pid in range(pop.len()):
        name = pop.param_name(pid)
        if name in cov:
            h.update(f"{name}=[{cov[name][0]!r},{cov[name][1]!r}];".encode())
        else:
            h.update(f"{name}={pop.param_get(pid)!r};".encode())
    nick = type(pop).__name__.removeprefix("GalaxyShapePop").lower() or "pop"
    kind = f"k{gsfs.get_pop_order()}" if cov else "fixed"
    return f"{nick}-{kind}-{h.hexdigest()[:16]}.gvar"


@contextlib.contextmanager
def _locked(path: Path):
    """Hold an exclusive lock on a file beside path, for jobs updating it at once."""
    lock_path = path.with_name(f".{path.name}.lock")
    lock_path.parent.mkdir(parents=True, exist_ok=True)
    with open(lock_path, "a", encoding="utf-8") as lock:
        fcntl.flock(lock.fileno(), fcntl.LOCK_EX)
        try:
            yield
        finally:
            fcntl.flock(lock.fileno(), fcntl.LOCK_UN)


def update_shape_interp_cache(
    data: Nc.DataClusterWLFactor,
    mset: Ncm.MSet,
    cache_dir: Path,
    *,
    coverage: dict[str, tuple[float, float]] | None = None,
    build: bool = True,
) -> ShapeInterpCacheUpdate:
    """Bring the cache of the fit of data and mset up to date in cache_dir, and attach it.

    The cache is the file :func:`shape_interp_cache_name` names in cache_dir. When it
    exists it is attached; then the galaxies without an interpolant, and those whose
    interpolant does not cover the shears they reach under the priors of mset
    (restrict-range), are built, in parallel, and the cache is written again: a cache built
    for wider priors is reused whole. coverage gives the population parameters the
    expansion spans (pop-order) and their bounds; the other population parameters are
    fixed. The update holds a lock beside the file, so that jobs started together build
    once. Afterwards obs carries the interpolants and names the cache, and data holds them
    as loaded ones.

    Raises:
        ValueError: The shape factor of data does not use interpolants; the file of this
            configuration exists but cannot be attached (it is never overwritten); or
            there is no cache and build is False.
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
    path = Path(cache_dir) / shape_interp_cache_name(obs, gsfs, pop, coverage)

    with _locked(path):
        strip_shape_interp(obs)
        if path.exists():
            cache = read_shape_interp_cache(path)
            reason = check_shape_interp_cache(cache, obs, gsfs, pop)
            if reason is None and not _same_coverage(cache, pop):
                reason = "its coverage of the population differs"
            if reason is not None:
                raise ValueError(f"refusing to overwrite {path}: {reason}")
            attach_shape_interp_cache(cache, obs, gsfs, pop)
        data.set_obs(obs)

        n0 = gsfs.get_n_built()
        u0 = gsfs.get_pop_n_unconverged()
        if not build:
            # Attached as it is: shears it does not cover are computed without it, counted
            if obs.peek_shape_interp_source() is None:
                raise ValueError(
                    f"there is no cache {path} for this fit, and building is off"
                )
        else:
            gsfs.set_rebuild_short(True)
            try:
                data.write_shape_interp(builder, obs)
            finally:
                gsfs.set_rebuild_short(False)
        built = gsfs.get_n_built() - n0
        unconverged = gsfs.get_pop_n_unconverged() - u0

        if built > 0 or not path.exists():
            write_shape_interp_cache(obs, path)
    name_shape_interp_cache(obs, path, wl_obs_fingerprint(obs))
    # Galaxies built here hold the interpolants of the coverage as built ones: obs is set
    # again so that every galaxy is a loaded one, as in any later load of the cache. With
    # nothing built they all are already, and setting obs would check every galaxy again
    # at the next prepare (10 s for the 12,460 expansions of s19a HWL16a-002).
    if built > 0:
        data.set_obs(obs)
    return ShapeInterpCacheUpdate(
        rows=obs.len(),
        built=built,
        reused=data.get_length() - built,
        unconverged=max(unconverged, 0),
        path=path,
    )


def _sub_catalog(obs: Nc.GalaxyWLObs, rows: list[int]) -> Nc.GalaxyWLObs:
    """Return the catalog of the rows of obs, with their p(z) and the metadata of obs."""
    cols = list(obs.peek_columns())
    sub = Nc.GalaxyWLObs.new(obs.get_ellip_conv(), obs.get_coord(), len(rows), cols)
    for k, i in enumerate(rows):
        for c in cols:
            sub.set(c, k, obs.get(c, i))
        sub.copy_pz(k, obs, i)
    if obs.peek_meta() is not None:
        sub.set_meta(obs.peek_meta())
    return sub


# Orders tried by choose_pop_order(), and the galaxies of the sample it builds
POP_ORDER_CANDIDATES = (7, 9, 11, 13, 15, 19, 23)
POP_ORDER_SAMPLE = 48


def choose_pop_order(
    data: Nc.DataClusterWLFactor,
    mset: Ncm.MSet,
    coverage: dict[str, tuple[float, float]],
    *,
    candidates: tuple[int, ...] = POP_ORDER_CANDIDATES,
    n_sample: int = POP_ORDER_SAMPLE,
) -> tuple[int, int]:
    """Return the lowest order of candidates whose expansions over coverage converge.

    The expansions of a sample of n_sample galaxies, spread in |eps_obs| from the smallest
    to the largest, are built at each order in turn, in parallel, until none fails to
    converge; their accuracy is then within interp-abstol. The width of the coverage, not
    the catalog, sets the order: on s19a HWL16a-002, Beta alpha and beta in [1, 4] need
    15, in [1.25, 2.75] 9, and a box of +-0.2 needs 7.

    Returns:
        The order and the number of galaxies of the sample whose expansion did not
        converge at it (non-zero only when no candidate converges).
    """
    obs = data.peek_obs()
    columns = list(obs.peek_columns())
    eps = np.hypot(
        obs.peek_data().to_numpy()[:, columns.index("epsilon_obs_1")],
        obs.peek_data().to_numpy()[:, columns.index("epsilon_obs_2")],
    )
    # The galaxies of the fit: the rows of the mask of data, or every row
    mask = data.get_mask()
    kept = np.flatnonzero(mask) if mask else np.arange(obs.len())
    order = kept[np.argsort(eps[kept])]
    n = min(n_sample, len(kept))
    rows = sorted({int(order[i]) for i in np.linspace(0, len(kept) - 1, n).astype(int)})
    sub = _sub_catalog(obs, rows)

    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    unconverged = n
    for k in candidates:
        probe = cast(Nc.DataClusterWLFactor, ser.dup_obj(data))
        ser.reset(True)
        gsfs = probe.props.shape_factor
        assert isinstance(gsfs, Nc.GalaxyShapeFactorSpectral)
        gsfs.set_pop_order(k)
        probe.set_mask(None)
        probe.set_obs(_sub_catalog(sub, list(range(sub.len()))))
        builder = _builder_mset(mset, coverage)
        probe.write_shape_interp(builder, probe.peek_obs())
        unconverged = gsfs.get_pop_n_unconverged()
        if unconverged == 0:
            return k, 0
    return candidates[-1], unconverged


def _same_coverage(cache: ShapeInterpCache, pop: Nc.GalaxyShapePop) -> bool:
    """Return whether the expansion of cache spans the free parameters of pop over their bounds."""
    if cache.n_modes <= 1:
        return True
    source = cache.source
    assert source is not None
    src = cast(Nc.GalaxyShapePop, source.peek("shape-pop"))
    for i in range(pop.len()):
        free = pop.param_get_ftype(i) == Ncm.ParamType.FREE
        if free != (src.param_get_ftype(i) == Ncm.ParamType.FREE):
            return False
        if free and not (
            _close(pop.param_get_lower_bound(i), src.param_get_lower_bound(i))
            and _close(pop.param_get_upper_bound(i), src.param_get_upper_bound(i))
        ):
            return False
    return True


def format_coverage(coverage: dict[str, tuple[float, float]] | None) -> str:
    """Return coverage as NAME=LO:HI;..., exact (repr) and sorted by name."""
    return ";".join(
        f"{k}={lo!r}:{hi!r}" for k, (lo, hi) in sorted((coverage or {}).items())
    )


def parse_coverage(text: str) -> dict[str, tuple[float, float]]:
    """Return the coverage of a string written by :func:`format_coverage`.

    Raises:
        ValueError: A malformed entry.
    """
    coverage: dict[str, tuple[float, float]] = {}
    for item in filter(None, text.split(";")):
        name, sep, bounds = item.partition("=")
        lo, sep2, hi = bounds.partition(":")
        if not (sep and sep2 and name):
            raise ValueError(f"malformed coverage {item!r}: expected NAME=LO:HI")
        coverage[name.strip()] = (float(lo), float(hi))
    return coverage


# Default population coverage, by population type: the box the expansion spans when
# --pop-coverage does not give one, intersected with the prior of each free parameter.
# Measured on 48 galaxies of s19a HWL16a-002 (2026-11-16): Beta [1.25, 2.75]^2 converges
# at pop-order 9 (the whole default prior [1, 4]^2 needs 15, ~2.4x the size), Gauss
# sigma in [0.25, 0.75] at 19, while its prior [0.1, 1] does not converge at 47. They cover
# the 5th-95th percentiles of the HSC cluster and strata fits.
DEFAULT_POP_COVERAGE: dict[str, dict[str, tuple[float, float]]] = {
    "GalaxyShapePopBeta": {"alpha": (1.25, 2.75), "beta": (1.25, 2.75)},
    "GalaxyShapePopGauss": {"sigma": (0.25, 0.75)},
}


def default_coverage(mset: Ncm.MSet) -> dict[str, tuple[float, float]] | None:
    """Return the default coverage of the free population parameters of mset, or None.

    Each free parameter of the population is covered over its prior, intersected with
    DEFAULT_POP_COVERAGE for its population type when that intersection is not empty: a
    prior narrower than the default box is covered whole, a wider one only over the box
    (outside it the galaxies fall back, with a warning). None when the population is
    fixed.
    """
    pop = mset.peek(Nc.GalaxyShapePop.id())
    if pop is None:
        return None
    box = DEFAULT_POP_COVERAGE.get(type(pop).__name__, {})
    coverage = {}
    for i in range(pop.len()):
        if pop.param_get_ftype(i) != Ncm.ParamType.FREE:
            continue
        name = pop.param_name(i)
        lo, hi = pop.param_get_lower_bound(i), pop.param_get_upper_bound(i)
        if name in box:
            blo, bhi = max(lo, box[name][0]), min(hi, box[name][1])
            if blo < bhi:
                lo, hi = blo, bhi
        coverage[name] = (lo, hi)
    return coverage or None


def coverage_short_of_prior(
    mset: Ncm.MSet, coverage: dict[str, tuple[float, float]] | None
) -> dict[str, tuple[float, float]]:
    """Return the priors of the free population parameters coverage does not span."""
    pop = mset.peek(Nc.GalaxyShapePop.id())
    short: dict[str, tuple[float, float]] = {}
    if pop is None:
        return short
    for i in range(pop.len()):
        if pop.param_get_ftype(i) != Ncm.ParamType.FREE:
            continue
        name = pop.param_name(i)
        lo, hi = pop.param_get_lower_bound(i), pop.param_get_upper_bound(i)
        clo, chi = (coverage or {}).get(name, (hi, lo))
        if clo > lo or chi < hi:
            short[name] = (lo, hi)
    return short


def record_shape_interp_cache(
    obs: Nc.GalaxyWLObs,
    path: Path,
    coverage: dict[str, tuple[float, float]] | None,
    base_dir: Path,
) -> None:
    """Name the cache of obs, relative to base_dir, and its coverage, then strip obs.

    What an experiment saves: the catalog without interpolants, naming the cache that
    :func:`load_shape_interp_cache` attaches (or builds) when the experiment is loaded.
    """
    name_shape_interp_cache(obs, path, wl_obs_fingerprint(obs), base_dir=base_dir)
    meta = obs.peek_meta()
    assert meta is not None
    meta.set_string(CATALOG_META_INTERP_COVERAGE, format_coverage(coverage))
    obs.set_meta(meta)
    strip_shape_interp(obs)


def recorded_coverage(obs: Nc.GalaxyWLObs) -> dict[str, tuple[float, float]] | None:
    """Return the coverage recorded by :func:`record_shape_interp_cache`, or None."""
    meta = obs.peek_meta()
    if meta is None:
        return None
    ok, text = meta.get_string(CATALOG_META_INTERP_COVERAGE)
    return (parse_coverage(text) or None) if ok else None


def load_shape_interp_cache(
    data: Nc.DataClusterWLFactor,
    mset: Ncm.MSet,
    base_dir: Path,
    *,
    build: bool = True,
) -> ShapeInterpCacheUpdate | None:
    """Attach the cache the catalog of data names, building what it misses.

    The cache is resolved against base_dir (the directory of the experiment) and brought up
    to date for the priors of mset by :func:`update_shape_interp_cache`, with the coverage
    recorded beside its name: a missing file, or galaxies the priors take beyond it, are
    built and saved, unless build is False, which makes a missing file an error. Returns
    None when the catalog names no cache, or carries its interpolants itself (a catalog
    written before the caches).

    Raises:
        ValueError: See :func:`update_shape_interp_cache`.
    """
    obs = data.peek_obs()
    path = resolve_shape_interp_cache(obs, base_dir)
    if path is None or obs.peek_shape_interp_source() is not None:
        return None
    gsfs = data.props.shape_factor
    if not isinstance(gsfs, Nc.GalaxyShapeFactorSpectral) or not gsfs.get_use_interp():
        return None
    return update_shape_interp_cache(
        data, mset, path.parent, coverage=recorded_coverage(obs), build=build
    )


# ---------------------------------------------------------------------------------------
# Caches of the redshift-node configuration (auto-nodes)
#
# The calibration of the redshift nodes of each galaxy (NcDataClusterWLFactor:auto-nodes)
# depends on the cluster (its redshift), the rows of the catalog with their p(z), the mask
# and the settings of the factors and of the nodes, unlike the interpolants, which do not
# depend on the cluster. It is kept in its own file beside the interpolant cache, named by
# a hash of all of these. The file holds NcDataClusterWLFactor:node-config, which the
# library checks again when it adopts it (its galaxy count, cluster redshift and node
# settings), recalibrating with a warning on any difference.
# ---------------------------------------------------------------------------------------

# Format of the node-config cache file, kept in its metadata
NODE_CONFIG_CACHE_FORMAT = "nc-node-config-cache-v1"
_NODE_META_KEY = "key"
# Prefix of the entries of the node-config in the metadata of the cache file
_NODE_META_PREFIX = "node-config:"


def wl_obs_full_fingerprint(obs: Nc.GalaxyWLObs) -> str:
    """Return a SHA-256 of every column and every p(z) table of obs.

    The calibration of the redshift nodes reads the positions, the shapes, the biases and
    the p(z) of each galaxy, so the fingerprint covers all of them (as little-endian
    float64), unlike :func:`wl_obs_fingerprint`.
    """
    h = hashlib.sha256()
    h.update(f"rows={obs.len()};".encode())
    h.update(f"conv={obs.get_ellip_conv().value_nick};".encode())
    h.update(f"coord={obs.get_coord().value_nick};".encode())
    columns = list(obs.peek_columns())
    if obs.len() > 0 and columns:
        data = obs.peek_data().to_numpy()
        for k, name in enumerate(columns):
            h.update(f"{name}:".encode())
            h.update(np.ascontiguousarray(data[:, k], dtype="<f8").tobytes())
    pz_table = obs.props.pz_table
    if pz_table is not None:
        for i in range(obs.len()):
            tab = pz_table.peek(i)
            h.update(f"pz{i}:".encode())
            if tab is not None:
                mat = cast(Ncm.Matrix, tab)
                h.update(np.ascontiguousarray(mat.to_numpy(), dtype="<f8").tobytes())
    return h.hexdigest()


def node_config_cache_name(data: Nc.DataClusterWLFactor, mset: Ncm.MSet) -> str:
    """Return the file name of the node-config cache of data at the cluster of mset.

    A hash of what the calibration depends on: every column and p(z) of the catalog
    (:func:`wl_obs_full_fingerprint`), the mask, the radial cuts, the node settings, the
    cluster redshift and the settings of the redshift, position and shape factors. The
    model point the calibration runs at is not part of it: any calibration reaching the
    tolerance is valid at every point.
    """
    obs = data.peek_obs()
    hp = mset.peek(Nc.HaloPosition.id())
    assert isinstance(hp, Nc.HaloPosition)
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    h = hashlib.sha256()
    h.update(wl_obs_full_fingerprint(obs).encode())
    mask = data.get_mask()
    h.update(
        f"mask={'' if mask is None else ''.join('1' if k else '0' for k in mask)};".encode()
    )
    h.update(
        (
            f"r=[{data.props.r_min!r},{data.props.r_max!r}];"
            f"method={data.get_integ_method().value_nick};auto={data.get_auto_nodes()};"
            f"reltol={data.get_node_reltol()!r};max={data.get_max_total_nodes()};"
            f"n={data.get_n_nodes()};rule={data.get_rule_n()};z_cl={hp['z']!r};"
        ).encode()
    )
    for factor in (
        data.props.redshift_factor,
        data.props.position_factor,
        data.props.shape_factor,
    ):
        h.update(ser.to_string(factor, True).encode())
        ser.reset(True)
    return f"nodes-{h.hexdigest()[:16]}.gvar"


def write_node_config_cache(
    obs: Nc.GalaxyWLObs, node_config: Ncm.VarDict, key: str, path: Path
) -> Path:
    """Write node_config to the cache file path, under key (its file name).

    The file is an NcGalaxyWLObs of one row and no columns (one with no rows cannot be
    read back) whose metadata holds the node-config.
    """
    holder = Nc.GalaxyWLObs.new(obs.get_ellip_conv(), obs.get_coord(), 1, [])
    meta = Ncm.VarDict.new()
    meta.set_string(_CACHE_META_FORMAT, NODE_CONFIG_CACHE_FORMAT)
    meta.set_string(_NODE_META_KEY, key)
    for name in node_config.keys():  # noqa: SIM118 (an NcmVarDict, not a dict)
        ok, value = node_config.get_variant(name)
        assert ok
        meta.set_variant(_NODE_META_PREFIX + name, value)
    holder.set_meta(meta)
    path = Path(path)
    _replace_atomically(
        path,
        lambda part: Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(
            holder, part.as_posix()
        ),
    )
    return path


def read_node_config_cache(path: Path) -> Ncm.VarDict:
    """Return the node-config of a file written by :func:`write_node_config_cache`.

    Raises:
        ValueError: The file is not a node-config cache of a known format, or was written
            under another name (another configuration).
    """
    path = Path(path)
    holder = Ncm.Serialize.new(Ncm.SerializeOpt.NONE).from_binfile(path.as_posix())
    meta = holder.peek_meta() if isinstance(holder, Nc.GalaxyWLObs) else None
    if meta is None:
        raise ValueError(f"{path} is not a node-config cache")
    ok_format, fmt = meta.get_string(_CACHE_META_FORMAT)
    ok_key, key = meta.get_string(_NODE_META_KEY)
    if not ok_format or fmt != NODE_CONFIG_CACHE_FORMAT:
        raise ValueError(
            f"{path} is not a node-config cache of format {NODE_CONFIG_CACHE_FORMAT} "
            f"(format: {fmt})"
        )
    if not ok_key or key != path.name:
        raise ValueError(f"{path} was written for another configuration ({key})")
    node_config = Ncm.VarDict.new()
    for name in meta.keys():  # noqa: SIM118 (an NcmVarDict, not a dict)
        if name.startswith(_NODE_META_PREFIX):
            ok, value = meta.get_variant(name)
            assert ok
            node_config.set_variant(name.removeprefix(_NODE_META_PREFIX), value)
    return node_config


@dataclasses.dataclass
class NodeConfigCacheUpdate:
    """What :func:`update_node_config_cache` did."""

    path: Path
    loaded: bool
    calibrated: int


def update_node_config_cache(
    data: Nc.DataClusterWLFactor,
    mset: Ncm.MSet,
    cache_dir: Path,
    *,
    prepare: bool = True,
) -> NodeConfigCacheUpdate | None:
    """Give data the node-config cached for this configuration, calibrating it if needed.

    The cache is the file :func:`node_config_cache_name` names in cache_dir. When it exists
    its node-config is given to data. With prepare, data is then prepared at mset
    (nc_data_cluster_wl_factor_data_prepare(), in parallel): that adopts the loaded
    configuration, or calibrates, in which case the result is written to the cache. A
    prepared dataset is also what a sampler shares with its copies, so preparing here
    avoids every copy of an MCMC calibrating. The update holds a lock beside the file, so
    that jobs started together calibrate once.

    Returns:
        None when data does not use auto-nodes (FIXED_NODES with auto-nodes), otherwise
        what was done.
    """
    if (
        data.get_integ_method() != Nc.DataClusterWLIntegMethod.FIXED_NODES
        or not data.get_auto_nodes()
    ):
        return None
    path = Path(cache_dir) / node_config_cache_name(data, mset)
    loaded = False
    with _locked(path):
        if path.exists():
            try:
                data.props.node_config = read_node_config_cache(path)
                loaded = True
            except ValueError as e:
                print(f"# Ignoring the node-config cache {path}: {e}")
        if not prepare:
            return NodeConfigCacheUpdate(path=path, loaded=loaded, calibrated=0)
        calib0 = data.get_calib_count()
        data.data_prepare(mset)
        calibrated = data.get_calib_count() - calib0
        if calibrated > 0:
            node_config = data.props.node_config
            assert node_config is not None
            write_node_config_cache(data.peek_obs(), node_config, path.name, path)
    return NodeConfigCacheUpdate(path=path, loaded=loaded, calibrated=calibrated)
