#
# test_wl_interp_cache.py
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

"""Tests of the caches of shape interpolants kept apart from the catalog."""

import json
from pathlib import Path
from typing import cast

import numpy as np
import pytest

from numcosmo_py import Nc, Ncm
from numcosmo_py.app.generate import parse_pop_dist, parse_shape_factor
from numcosmo_py.experiments.cluster_wl import (
    DEFAULT_INTEG_OPTIONS,
    apply_prior,
    centre_box_bounds,
    load_cluster_wl,
)
from numcosmo_py.experiments.wl_interp_cache import (
    CATALOG_META_INTERP_CACHE,
    SHAPE_INTERP_CACHE_FORMAT,
    attach_shape_interp_cache,
    check_shape_interp_cache,
    find_shape_interp_cache,
    read_shape_interp_cache,
    read_shape_interp_cache_info,
    strip_shape_interp,
    update_shape_interp_cache,
    wl_obs_fingerprint,
)

Ncm.cfg_init()

SPECTRAL = "spectral ellip_conv=trace ellip_coord=celestial use_interp=true restrict_range=true"


def _catalog(n: int = 12) -> Nc.GalaxyWLObs:
    """A small catalog around a cluster at (0, 0), z = 0.3, with p(z) splines."""
    rng = np.random.default_rng(20261116)
    cols = [
        "epsilon_int_1",
        "epsilon_int_2",
        "epsilon_obs_1",
        "epsilon_obs_2",
        "std_noise",
        "c1",
        "c2",
        "m",
        "ra",
        "dec",
        "z",
    ]
    obs = Nc.GalaxyWLObs.new(
        Nc.GalaxyWLObsEllipConv.TRACE, Nc.WLEllipticityFrame.CELESTIAL, n, cols
    )
    ra = rng.uniform(-0.05, 0.05, n)
    dec = rng.uniform(-0.05, 0.05, n)
    z = rng.uniform(0.6, 1.2, n)
    for i in range(n):
        obs.set("epsilon_int_1", i, 0.0)
        obs.set("epsilon_int_2", i, 0.0)
        obs.set("epsilon_obs_1", i, float(rng.normal(0.0, 0.2)))
        obs.set("epsilon_obs_2", i, float(rng.normal(0.0, 0.2)))
        obs.set("std_noise", i, float(rng.uniform(0.1, 0.2)))
        obs.set("c1", i, 0.0)
        obs.set("c2", i, 0.0)
        obs.set("m", i, 0.0)
        obs.set("ra", i, float(ra[i]))
        obs.set("dec", i, float(dec[i]))
        obs.set("z", i, float(z[i]))
        zs = np.linspace(max(z[i] - 0.3, 0.01), z[i] + 0.3, 16)
        spline = Ncm.SplineCubicNotaknot.new()
        spline.set(
            Ncm.Vector.new_array(zs.tolist()),
            Ncm.Vector.new_array(np.exp(-0.5 * ((zs - z[i]) / 0.05) ** 2).tolist()),
            True,
        )
        obs.set_pz(i, spline)
    meta = Ncm.VarDict.new()
    for key, value in (
        ("cluster_ra", 0.0),
        ("cluster_dec", 0.0),
        ("cluster_z", 0.3),
        ("cluster_c", 4.0),
    ):
        meta.set_double(key, value)
    obs.set_meta(meta)
    # Through a serialization, as a catalog read from a file
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.NONE)
    return cast(Nc.GalaxyWLObs, ser.from_variant(ser.to_variant(obs)))


def _experiment(
    obs: Nc.GalaxyWLObs, mass_max: float = 1.0e15, shape: str = SPECTRAL
) -> tuple[Nc.DataClusterWLFactor, Ncm.MSet, Ncm.Dataset]:
    """The likelihood of obs, mass and centre free (within 0.005 deg), with mass up to mass_max."""
    _sf, shape_gen = parse_shape_factor(shape)
    _pd, pop_gen = parse_pop_dist("gauss sigma=0.3")
    exp = load_cluster_wl(
        obs=obs,
        cluster_mass=1.0e14,
        cluster_mass_min=1.0e13,
        cluster_mass_max=mass_max,
        r_min=0.3 / 0.7,
        r_max=3.0 / 0.7,
        shape_gen=shape_gen,
        pop_gen=pop_gen,
        integ_options=DEFAULT_INTEG_OPTIONS,
        summary=False,
    )
    mset = cast(Ncm.MSet, exp.peek("model-set"))
    # A centre box smaller than the distance of the galaxies, so that their reach is bounded
    ra_lo, ra_hi, dec_lo, dec_hi = centre_box_bounds(mset, half_deg=0.005)
    apply_prior(mset, "NcHaloPosition:ra", ra_lo, ra_hi)
    apply_prior(mset, "NcHaloPosition:dec", dec_lo, dec_hi)
    for name in (
        "NcHaloMassSummary:log10MDelta",
        "NcHaloPosition:ra",
        "NcHaloPosition:dec",
    ):
        pindex = mset.param_get_by_full_name(name)
        mset.param_set_ftype(pindex.mid, pindex.pid, Ncm.ParamType.FREE)
    mset.prepare_fparam_map()
    dataset = cast(Ncm.Likelihood, exp.peek("likelihood")).peek_dataset()
    data = cast(Nc.DataClusterWLFactor, dataset.get_data(0))
    return data, mset, dataset


def _gsfs_pop(data, mset) -> tuple[Nc.GalaxyShapeFactorSpectral, Nc.GalaxyShapePop]:
    """The shape factor of data and the population of mset."""
    return data.props.shape_factor, mset.peek(Nc.GalaxyShapePop.id())


def test_update_and_attach(tmp_path: Path):
    """A cache built for a fit gives, attached to a fresh copy of the catalog, the same
    -2lnL bit for bit and builds nothing; a second update reuses every row."""
    path = tmp_path / "cache.gvar"
    data, mset, dataset = _experiment(_catalog())
    upd = update_shape_interp_cache(data, mset, path)
    assert upd.rows == 12
    assert upd.built == 12
    assert upd.reused == 0
    m2lnL = dataset.m2lnL_val(mset)

    info = read_shape_interp_cache_info(path)
    assert info["format"] == SHAPE_INTERP_CACHE_FORMAT
    assert info["rows_with_interpolants"] == 12
    assert info["shape_factor"]["restrict_range"] is True

    obs2 = _catalog()
    data2, mset2, dataset2 = _experiment(obs2)
    gsfs2, pop2 = _gsfs_pop(data2, mset2)
    cache = find_shape_interp_cache(tmp_path, obs2, gsfs2, pop2)
    assert cache is not None
    attach_shape_interp_cache(cache, obs2, gsfs2, pop2)
    data2.set_obs(obs2)
    assert dataset2.m2lnL_val(mset2) == m2lnL
    assert gsfs2.get_n_built() == 0
    ok, name = obs2.peek_meta().get_string(CATALOG_META_INTERP_CACHE)
    assert ok and name == path.as_posix()

    upd = update_shape_interp_cache(data2, mset2, path)
    assert (upd.built, upd.reused) == (0, 12)
    assert dataset2.m2lnL_val(mset2) == m2lnL


def test_update_wider_prior_rebuilds_short(tmp_path: Path):
    """With a wider mass prior the galaxies reach larger shears: the interpolants short
    of them, and only those, are built again, and the values are those of a cache built
    for the wider prior from scratch. A narrower prior then reuses the wider cache."""
    path = tmp_path / "cache.gvar"
    data, mset, dataset = _experiment(_catalog(), mass_max=1.0e14)
    update_shape_interp_cache(data, mset, path)
    m2lnL_narrow = dataset.m2lnL_val(mset)

    data2, mset2, dataset2 = _experiment(_catalog(), mass_max=1.0e16)
    gsfs2, _pop2 = _gsfs_pop(data2, mset2)
    upd = update_shape_interp_cache(data2, mset2, path)
    assert upd.built > 0
    assert upd.built == gsfs2.get_n_loaded_short()
    m2lnL = dataset2.m2lnL_val(mset2)

    data3, mset3, dataset3 = _experiment(_catalog(), mass_max=1.0e16)
    update_shape_interp_cache(data3, mset3, tmp_path / "fresh.gvar")
    assert dataset3.m2lnL_val(mset3) == m2lnL

    # Contained: nothing to build, and the values are those of the wider interpolants
    data4, mset4, dataset4 = _experiment(_catalog(), mass_max=1.0e14)
    upd = update_shape_interp_cache(data4, mset4, path)
    assert (upd.built, upd.reused) == (0, 12)
    assert dataset4.m2lnL_val(mset4) == pytest.approx(m2lnL_narrow, abs=1.0e-6)


def test_cache_refused(tmp_path: Path):
    """Other rows (the fingerprint) or other settings (the record) are refused,
    without aborting, and lookup finds nothing."""
    path = tmp_path / "cache.gvar"
    data, mset, _dataset = _experiment(_catalog())
    update_shape_interp_cache(data, mset, path)
    cache = read_shape_interp_cache(path)

    other = _catalog()
    other.set("epsilon_obs_1", 3, other.get("epsilon_obs_1", 3) + 1.0e-9)
    data_o, mset_o, _ = _experiment(other)
    gsfs_o, pop_o = _gsfs_pop(data_o, mset_o)
    assert wl_obs_fingerprint(other) != cache.fingerprint
    assert "fingerprint" in check_shape_interp_cache(cache, other, gsfs_o, pop_o)
    with pytest.raises(ValueError, match="fingerprint"):
        attach_shape_interp_cache(cache, other, gsfs_o, pop_o)
    assert find_shape_interp_cache(tmp_path, other, gsfs_o, pop_o) is None

    obs_s = _catalog()
    data_s, mset_s, _ = _experiment(obs_s, shape=SPECTRAL + " interp_abstol=1e-6")
    gsfs_s, pop_s = _gsfs_pop(data_s, mset_s)
    reason = check_shape_interp_cache(cache, obs_s, gsfs_s, pop_s)
    assert reason is not None and "interp-abstol" in reason
    assert find_shape_interp_cache(tmp_path, obs_s, gsfs_s, pop_s) is None


def test_check_record_without_record():
    """A record without the shape factor and the population is incompatible."""
    data, mset, _ = _experiment(_catalog())
    gsfs, pop = _gsfs_pop(data, mset)
    ok, diff = gsfs.check_record(pop, Ncm.ObjDictStr.new(), 0)
    assert not ok
    assert "lacks" in diff


def test_strip_and_bad_files(tmp_path: Path):
    """Stripping removes the interpolants and their record; a file that is not a cache,
    or the catalog itself, is refused."""
    path = tmp_path / "cache.gvar"
    data, mset, _ = _experiment(_catalog())
    update_shape_interp_cache(data, mset, path)
    obs = data.peek_obs()
    assert obs.peek_shape_interp_source() is not None
    strip_shape_interp(obs)
    assert obs.peek_shape_interp_source() is None
    assert all(obs.peek_shape_interp(i) is None for i in range(obs.len()))

    not_cache = tmp_path / "catalog.gvar"
    Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(
        _catalog(), not_cache.as_posix()
    )
    with pytest.raises(ValueError, match="format"):
        read_shape_interp_cache(not_cache)
    with pytest.raises(ValueError, match="no current interpolants"):
        update_shape_interp_cache.__globals__["write_shape_interp_cache"](
            _catalog(), tmp_path / "x.gvar"
        )
    assert json.loads(path.with_suffix(".json").read_text())["rows"] == 12
