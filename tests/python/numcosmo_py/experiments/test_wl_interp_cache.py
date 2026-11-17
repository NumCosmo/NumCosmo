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
    choose_pop_order,
    coverage_short_of_prior,
    default_coverage,
    name_shape_interp_cache,
    node_config_cache_name,
    read_node_config_cache,
    update_node_config_cache,
    read_shape_interp_cache,
    resolve_shape_interp_cache,
    shape_interp_cache_info,
    shape_interp_cache_name,
    strip_shape_interp,
    update_shape_interp_cache,
    wl_obs_fingerprint,
    write_shape_interp_cache,
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
    obs: Nc.GalaxyWLObs,
    mass_max: float = 1.0e15,
    shape: str = SPECTRAL,
    pop: str = "gauss sigma=0.3",
) -> tuple[Nc.DataClusterWLFactor, Ncm.MSet, Ncm.Dataset]:
    """The likelihood of obs, mass and centre free (within 0.005 deg), with mass up to mass_max."""
    _sf, shape_gen = parse_shape_factor(shape)
    _pd, pop_gen = parse_pop_dist(pop)
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
    data, mset, dataset = _experiment(_catalog())
    upd = update_shape_interp_cache(data, mset, tmp_path)
    assert (upd.rows, upd.built, upd.reused, upd.unconverged) == (12, 12, 0, 0)
    assert upd.path.parent == tmp_path
    assert upd.path.name.startswith("gauss-fixed-")
    m2lnL = dataset.m2lnL_val(mset)

    cache = read_shape_interp_cache(upd.path)
    info = shape_interp_cache_info(cache.holder, cache.fingerprint)
    assert info["format"] == SHAPE_INTERP_CACHE_FORMAT
    assert info["rows_with_interpolants"] == 12
    assert info["shape_factor"]["restrict_range"] is True

    obs2 = _catalog()
    data2, mset2, dataset2 = _experiment(obs2)
    gsfs2, pop2 = _gsfs_pop(data2, mset2)
    assert upd.path.name == shape_interp_cache_name(obs2, gsfs2, pop2, None)
    attach_shape_interp_cache(cache, obs2, gsfs2, pop2)
    data2.set_obs(obs2)
    assert dataset2.m2lnL_val(mset2) == m2lnL
    assert gsfs2.get_n_built() == 0
    ok, name = obs2.peek_meta().get_string(CATALOG_META_INTERP_CACHE)
    assert ok and name == upd.path.as_posix()

    upd2 = update_shape_interp_cache(data2, mset2, tmp_path)
    assert (upd2.built, upd2.reused, upd2.path) == (0, 12, upd.path)
    assert dataset2.m2lnL_val(mset2) == m2lnL

    # Without building, the cache is attached as it is; with no cache it is an error
    data3, mset3, dataset3 = _experiment(_catalog())
    update_shape_interp_cache(data3, mset3, tmp_path, build=False)
    assert dataset3.m2lnL_val(mset3) == m2lnL
    data4, mset4, _ = _experiment(_catalog())
    with pytest.raises(ValueError, match="building is off"):
        update_shape_interp_cache(data4, mset4, tmp_path / "empty", build=False)


def test_update_wider_prior_rebuilds_short(tmp_path: Path):
    """With a wider mass prior the galaxies reach larger shears: the interpolants short
    of them, and only those, are built again, and the values are those of a cache built
    for the wider prior from scratch. A narrower prior then reuses the wider cache."""
    data, mset, dataset = _experiment(_catalog(), mass_max=1.0e14)
    update_shape_interp_cache(data, mset, tmp_path)
    m2lnL_narrow = dataset.m2lnL_val(mset)

    data2, mset2, dataset2 = _experiment(_catalog(), mass_max=1.0e16)
    gsfs2, _pop2 = _gsfs_pop(data2, mset2)
    upd = update_shape_interp_cache(data2, mset2, tmp_path)
    assert upd.built > 0
    assert upd.built == gsfs2.get_n_loaded_short()
    m2lnL = dataset2.m2lnL_val(mset2)

    data3, mset3, dataset3 = _experiment(_catalog(), mass_max=1.0e16)
    update_shape_interp_cache(data3, mset3, tmp_path / "fresh")
    assert dataset3.m2lnL_val(mset3) == m2lnL

    # Contained: nothing to build, and the values are those of the wider interpolants
    data4, mset4, dataset4 = _experiment(_catalog(), mass_max=1.0e14)
    upd = update_shape_interp_cache(data4, mset4, tmp_path)
    assert (upd.built, upd.reused) == (0, 12)
    assert dataset4.m2lnL_val(mset4) == pytest.approx(m2lnL_narrow, abs=1.0e-6)


def test_expansion_rebuilds_short(tmp_path: Path):
    """An expansion over alpha and beta of Beta: a wider mass prior builds the short
    rows again, as expansions, without aborting, and gives the values of an expansion
    cache built for it from scratch, inside the coverage and at its edge."""
    beta = "beta alpha=1.58 beta=1.69"
    shape = SPECTRAL + " pop_order=5"
    cov = {"alpha": (1.48, 1.68), "beta": (1.59, 1.79)}
    data, mset, _ = _experiment(_catalog(), mass_max=1.0e14, shape=shape, pop=beta)
    upd = update_shape_interp_cache(data, mset, tmp_path, coverage=cov)
    assert upd.built == 12 and upd.unconverged == 0
    assert upd.path.name.startswith("beta-k5-")

    data2, mset2, dataset2 = _experiment(
        _catalog(), mass_max=1.0e16, shape=shape, pop=beta
    )
    gsfs2, _ = _gsfs_pop(data2, mset2)
    upd2 = update_shape_interp_cache(data2, mset2, tmp_path, coverage=cov)
    assert upd2.path == upd.path
    assert 0 < upd2.built == gsfs2.get_n_loaded_short()

    data3, mset3, dataset3 = _experiment(
        _catalog(), mass_max=1.0e16, shape=shape, pop=beta
    )
    update_shape_interp_cache(data3, mset3, tmp_path / "fresh", coverage=cov)
    for m, ds in ((mset2, dataset2), (mset3, dataset3)):
        for name in ("alpha", "beta"):
            pindex = m.param_get_by_full_name(f"NcGalaxyShapePop:{name}")
            m.param_set_ftype(pindex.mid, pindex.pid, Ncm.ParamType.FREE)
            model = m.peek(pindex.mid)
            model.param_set_lower_bound(pindex.pid, cov[name][0])
            model.param_set_upper_bound(pindex.pid, cov[name][1])
        m.prepare_fparam_map()
    for a, b in ((1.58, 1.69), (1.5, 1.75), (1.68, 1.59)):
        vals = []
        for m, ds in ((mset2, dataset2), (mset3, dataset3)):
            for name, x in (("alpha", a), ("beta", b)):
                pindex = m.param_get_by_full_name(f"NcGalaxyShapePop:{name}")
                m.param_set(pindex.mid, pindex.pid, x)
            vals.append(ds.m2lnL_val(m))
        assert vals[0] == vals[1]


def test_configurations_do_not_share_files(tmp_path: Path):
    """Caches of two populations coexist in one directory, and a file of the right name
    that cannot be attached is refused, never overwritten."""
    data_g, mset_g, _ = _experiment(_catalog())
    upd_g = update_shape_interp_cache(data_g, mset_g, tmp_path)
    data_b, mset_b, _ = _experiment(_catalog(), pop="beta alpha=1.58 beta=1.69")
    upd_b = update_shape_interp_cache(data_b, mset_b, tmp_path)
    assert upd_g.path != upd_b.path
    gauss_bytes = upd_g.path.read_bytes()

    # The Beta cache written under the name of the Gauss one
    upd_g.path.write_bytes(upd_b.path.read_bytes())
    data_g2, mset_g2, _ = _experiment(_catalog())
    with pytest.raises(ValueError, match="refusing to overwrite"):
        update_shape_interp_cache(data_g2, mset_g2, tmp_path)
    assert upd_g.path.read_bytes() == upd_b.path.read_bytes()
    upd_g.path.write_bytes(gauss_bytes)


def test_relative_name(tmp_path: Path):
    """The cache named relative to a directory is found from it after both move."""
    data, mset, _ = _experiment(_catalog())
    upd = update_shape_interp_cache(data, mset, tmp_path / "a" / "interp")
    obs = data.peek_obs()
    name_shape_interp_cache(
        obs, upd.path, wl_obs_fingerprint(obs), base_dir=tmp_path / "a"
    )
    ok, name = obs.peek_meta().get_string(CATALOG_META_INTERP_CACHE)
    assert ok and name == f"interp/{upd.path.name}"
    (tmp_path / "a").rename(tmp_path / "b")
    path = resolve_shape_interp_cache(obs, tmp_path / "b")
    assert path is not None and path.exists()


def test_choose_pop_order():
    """The order chosen is the lowest whose expansions of the sample converge: a narrow
    coverage takes the first candidate, the wide [1, 4] a higher one."""
    data, mset, _ = _experiment(
        _catalog(), shape=SPECTRAL + " pop_order=5", pop="beta alpha=1.58 beta=1.69"
    )
    k, unconv = choose_pop_order(
        data, mset, {"alpha": (1.57, 1.59)}, candidates=(5, 9), n_sample=6
    )
    assert (k, unconv) == (5, 0)
    k, unconv = choose_pop_order(
        data, mset, {"alpha": (1.0, 4.0)}, candidates=(5, 9, 15), n_sample=6
    )
    assert k > 5 and unconv == 0


def test_cache_refused(tmp_path: Path):
    """Other rows (the fingerprint) or other settings (the record) are refused,
    without aborting, and give other cache names."""
    data, mset, _dataset = _experiment(_catalog())
    upd = update_shape_interp_cache(data, mset, tmp_path)
    cache = read_shape_interp_cache(upd.path)

    other = _catalog()
    other.set("epsilon_obs_1", 3, other.get("epsilon_obs_1", 3) + 1.0e-9)
    data_o, mset_o, _ = _experiment(other)
    gsfs_o, pop_o = _gsfs_pop(data_o, mset_o)
    assert wl_obs_fingerprint(other) != cache.fingerprint
    assert "fingerprint" in check_shape_interp_cache(cache, other, gsfs_o, pop_o)
    with pytest.raises(ValueError, match="fingerprint"):
        attach_shape_interp_cache(cache, other, gsfs_o, pop_o)
    assert shape_interp_cache_name(other, gsfs_o, pop_o, None) != upd.path.name

    obs_s = _catalog()
    data_s, mset_s, _ = _experiment(obs_s, shape=SPECTRAL + " interp_abstol=1e-6")
    gsfs_s, pop_s = _gsfs_pop(data_s, mset_s)
    reason = check_shape_interp_cache(cache, obs_s, gsfs_s, pop_s)
    assert reason is not None and "interp-abstol" in reason
    assert shape_interp_cache_name(obs_s, gsfs_s, pop_s, None) != upd.path.name


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
    data, mset, _ = _experiment(_catalog())
    upd = update_shape_interp_cache(data, mset, tmp_path)
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
        write_shape_interp_cache(_catalog(), tmp_path / "x.gvar")
    assert read_shape_interp_cache(upd.path).n_rows == 12


def _free_pop(mset: Ncm.MSet, bounds: dict[str, tuple[float, float]]) -> None:
    """Free the population parameters of bounds over them."""
    pop = mset.peek(Nc.GalaxyShapePop.id())
    for name, (lo, hi) in bounds.items():
        _, pid = pop.param_index_from_name(name)
        pop.param_set_lower_bound(pid, lo)
        pop.param_set_upper_bound(pid, hi)
        pop.param_set_ftype(pid, Ncm.ParamType.FREE)


@pytest.mark.parametrize(
    "pop, prior, expected",
    [
        ("beta alpha=1.58 beta=1.69", None, None),
        (
            "beta alpha=1.58 beta=1.69",
            {"alpha": (1.0, 4.0), "beta": (1.0, 4.0)},
            {"alpha": (1.25, 2.75), "beta": (1.25, 2.75)},
        ),
        (
            "beta alpha=1.58 beta=1.69",
            {"alpha": (1.5, 1.7), "beta": (1.0, 2.0)},
            {"alpha": (1.5, 1.7), "beta": (1.25, 2.0)},
        ),
        ("gauss sigma=0.3", {"sigma": (0.1, 1.0)}, {"sigma": (0.25, 0.75)}),
        # A prior that misses the default box is covered whole
        ("gauss sigma=0.9", {"sigma": (0.8, 1.0)}, {"sigma": (0.8, 1.0)}),
    ],
)
def test_default_coverage(pop, prior, expected):
    """The default coverage is the prior of each free parameter within the default box of
    its population; coverage_short_of_prior() names the priors it does not span."""
    _, mset, _ = _experiment(_catalog(), pop=pop)
    if prior is not None:
        _free_pop(mset, prior)
    coverage = default_coverage(mset)
    assert coverage == expected
    short = coverage_short_of_prior(mset, coverage)
    assert set(short) == {
        k for k in (prior or {}) if (prior or {})[k] != (expected or {})[k]
    }


def _fixed_nodes(data: Nc.DataClusterWLFactor) -> Nc.DataClusterWLFactor:
    """data with FIXED_NODES redshift integrals, whose auto-nodes is on by default."""
    data.set_integ_method(Nc.DataClusterWLIntegMethod.FIXED_NODES)
    assert data.get_auto_nodes()
    return data


def test_node_config_cache(tmp_path: Path):
    """The redshift-node configuration calibrated for a fit is saved and, for a fresh copy
    of the experiment, adopted without calibrating, with the same -2lnL bit for bit."""
    data, mset, dataset = _experiment(_catalog())
    _fixed_nodes(data)
    upd = update_node_config_cache(data, mset, tmp_path)
    assert upd is not None
    assert (upd.loaded, upd.calibrated) == (False, data.get_length())
    assert upd.path.parent == tmp_path and upd.path.name.startswith("nodes-")
    assert upd.path.exists()
    assert data.is_data_prepared()
    m2lnL = dataset.m2lnL_val(mset)

    data2, mset2, dataset2 = _experiment(_catalog())
    _fixed_nodes(data2)
    upd2 = update_node_config_cache(data2, mset2, tmp_path)
    assert upd2 is not None
    assert (upd2.loaded, upd2.calibrated, upd2.path) == (True, 0, upd.path)
    assert dataset2.m2lnL_val(mset2) == m2lnL
    assert data2.get_calib_count() == 0

    # Without preparing, the configuration is only given to the data
    data3, mset3, dataset3 = _experiment(_catalog())
    _fixed_nodes(data3)
    upd3 = update_node_config_cache(data3, mset3, tmp_path, prepare=False)
    assert upd3 is not None and (upd3.loaded, upd3.calibrated) == (True, 0)
    assert dataset3.m2lnL_val(mset3) == m2lnL
    assert data3.get_calib_count() == 0

    # Not auto-nodes, or not FIXED_NODES: nothing to cache
    data4, mset4, _ = _experiment(_catalog())
    assert update_node_config_cache(data4, mset4, tmp_path) is None
    _fixed_nodes(data4).set_auto_nodes(False)
    assert update_node_config_cache(data4, mset4, tmp_path) is None


def test_node_config_cache_names(tmp_path: Path):
    """What the calibration depends on names another file: the node settings, the cluster
    redshift, the mask and the p(z) of a galaxy."""
    data, mset, _ = _experiment(_catalog())
    _fixed_nodes(data)
    name = node_config_cache_name(data, mset)
    assert name == node_config_cache_name(data, mset)

    reltol = data.get_node_reltol()
    data.set_node_reltol(reltol * 10.0)
    assert node_config_cache_name(data, mset) != name
    data.set_node_reltol(reltol)
    assert node_config_cache_name(data, mset) == name

    hp = mset.peek(Nc.HaloPosition.id())
    hp["z"] = 0.31
    assert node_config_cache_name(data, mset) != name
    hp["z"] = 0.3

    data.set_mask([i != 3 for i in range(data.peek_obs().len())])
    assert node_config_cache_name(data, mset) != name
    data.set_mask(None)
    assert node_config_cache_name(data, mset) == name

    obs = _catalog()
    zs = np.linspace(0.5, 1.1, 16)
    spline = Ncm.SplineCubicNotaknot.new()
    spline.set(
        Ncm.Vector.new_array(zs.tolist()),
        Ncm.Vector.new_array(np.exp(-0.5 * ((zs - 0.8) / 0.06) ** 2).tolist()),
        True,
    )
    obs.set_pz(5, spline)
    data2, mset2, _ = _experiment(obs)
    _fixed_nodes(data2)
    assert node_config_cache_name(data2, mset2) != name


def test_node_config_cache_refused(tmp_path: Path):
    """A file of another configuration under the name of this one, or a file that is not
    a cache, is ignored: the data calibrates and the file is written again."""
    data, mset, _ = _experiment(_catalog())
    _fixed_nodes(data)
    upd = update_node_config_cache(data, mset, tmp_path)
    assert upd is not None

    data2, mset2, _ = _experiment(_catalog())
    _fixed_nodes(data2).set_node_reltol(data2.get_node_reltol() * 10.0)
    other = tmp_path / node_config_cache_name(data2, mset2)
    other.write_bytes(upd.path.read_bytes())
    with pytest.raises(ValueError, match="another configuration"):
        read_node_config_cache(other)
    upd2 = update_node_config_cache(data2, mset2, tmp_path)
    assert upd2 is not None
    assert (upd2.path, upd2.loaded, upd2.calibrated) == (
        other,
        False,
        data2.get_length(),
    )
    assert read_node_config_cache(other) is not None

    bad = tmp_path / "nodes-bad.gvar"
    Ncm.Serialize.new(Ncm.SerializeOpt.NONE).to_binfile(_catalog(), bad.as_posix())
    with pytest.raises(ValueError, match="not a node-config cache"):
        read_node_config_cache(bad)
