"""Tests for numcosmo_py.datasets.dist_mu_obs."""

from pathlib import Path

import numpy as np
import pytest
from astropy.table import Table

from numcosmo_py import Nc, Ncm
from numcosmo_py.datasets.dist_mu_obs import (
    SPEED_OF_LIGHT_KM_S,
    DistMuObsCatalogID,
    SampleRule,
    dist_mu_obs_from_table,
    edd_cf4_snia_to_dist_mu_obs,
    samples_by_rule,
    survey_names,
)

Ncm.cfg_init()


def test_dist_mu_obs_from_table():
    """Table columns, types, NaN for masked cells and the metadata."""
    table = Table()
    table["ra"] = np.array([10.0, 20.0, 30.0])
    table["dec"] = np.array([-5.0, 0.0, 5.0])
    table["z"] = np.array([0.01, 0.02, 0.03])
    table["mu"] = np.array([33.0, 34.0, 35.0])
    table["sigma_mu"] = np.array([0.1, 0.2, 0.3])
    table["pgc"] = np.array([1, 2, 3], dtype=np.int64)
    table["ok"] = np.array([True, False, True])
    table["name"] = np.array(["a", "b", "c"])
    table["extra"] = np.ma.masked_array([1.0, 2.0, 3.0], mask=[False, True, False])

    obs = dist_mu_obs_from_table(table, {"source": "unit test"})

    assert isinstance(obs, Nc.DistMuObs)
    assert obs.len() == 3
    assert not obs.has_column("name")
    assert obs.get_col_type("pgc") == Ncm.CatalogColType.INT
    assert obs.get_col_type("ok") == Ncm.CatalogColType.BOOL
    assert obs.get_int("pgc", 2) == 3
    assert obs.get_bool("ok", 1) is False
    assert np.isnan(obs.get("extra", 1))
    assert obs.get("extra", 2) == 3.0
    meta = obs.peek_meta()
    assert meta is not None
    assert meta.get_string("source") == (True, "unit test")


def test_edd_cf4_snia_to_dist_mu_obs(edd_cf4_excerpt_file: Path):
    """The EDD excerpt gives the catalog columns and the metadata."""
    obs = edd_cf4_snia_to_dist_mu_obs(edd_cf4_excerpt_file, "2026-10-17")

    assert obs.len() == 6
    for col in ("ra", "dec", "z", "mu", "sigma_mu", "pgc", "1PGC", "Vcmb", "DM_av"):
        assert obs.has_column(col)
    for col in ("RA", "DE", "DM_zp", "eDM"):
        assert not obs.has_column(col)
    assert obs.get_col_type("pgc") == Ncm.CatalogColType.INT

    # First row of the export: pgc 250, RA 0.8958, DE 23.2008, Vcmb 6914, DM_zp 34.912, eDM 0.13
    assert obs.get_int("pgc", 0) == 250
    assert obs.get("ra", 0) == 0.8958
    assert obs.get("dec", 0) == 23.2008
    assert obs.get("z", 0) == pytest.approx(6914.0 / SPEED_OF_LIGHT_KM_S, rel=1.0e-15)
    assert obs.get("mu", 0) == 34.912
    assert obs.get("sigma_mu", 0) == 0.13

    meta = obs.peek_meta()
    assert meta is not None
    assert meta.get_string("downloaded") == (True, "2026-10-17")
    assert meta.get_string("source")[1].startswith("All CF4 SNIa Samples")

    n = obs.get_direction(0)
    assert n.norm() == pytest.approx(1.0, rel=1.0e-15)


def test_catalog_id_enum():
    """The catalog ids map to the C enum."""
    assert DistMuObsCatalogID.CF4_SNIA.genum == Nc.DistMuObsCatalogId.CF4_SNIA


def test_samples_by_rule(edd_cf4_excerpt_file: Path):
    """Sample rules on the EDD excerpt."""
    obs = edd_cf4_snia_to_dist_mu_obs(edd_cf4_excerpt_file, "2026-10-17")
    names = survey_names(obs)
    assert "pantheonplus" in names and "stahl" in names
    assert names == sorted(names)

    samples, labels = samples_by_rule(obs, SampleRule.SMALLEST_ERROR)
    assert labels == names
    assert len(samples) == obs.len()
    columns = list(obs.peek_columns())
    data = obs.peek_data().to_numpy()
    for i, s in enumerate(samples):
        errors = {n: data[i, columns.index(f"eDM_{n}")] for n in names}
        finite = {n: e for n, e in errors.items() if np.isfinite(e)}
        if finite:
            assert names[int(s)] == min(finite, key=finite.__getitem__)
        else:
            assert s == -1.0

    samples, labels = samples_by_rule(obs, SampleRule.HEMISPHERE, np.array([0, 1]))
    assert labels == ["north", "south"]
    assert samples == [0.0 if obs.get("dec", i) >= 0 else 1.0 for i in (0, 1)]

    samples, labels = samples_by_rule(obs, SampleRule.NONE)
    assert labels == [] and set(samples) == {-1.0}
