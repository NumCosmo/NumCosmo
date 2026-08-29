#!/usr/bin/env python
#
# test_planck_release.py
#
# Thu July 23 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_planck_release.py
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

"""Tests for the curated native Planck likelihood release (build + local load)."""

import os
import shutil

import numpy as np
import pytest

from python.fixtures_planck import make_commander_cldf, planck_mset
from numcosmo_py import Ncm, Nc
from numcosmo_py.experiments.planck_lite import find_baseline_file
from numcosmo_py.experiments.planck_commander import COMMANDER_RELPATH, build_commander
import numcosmo_py.experiments.planck_native_release as pnr
from numcosmo_py.experiments.planck_native_release import (
    PlanckReleaseId,
    build_release,
    load_planck_release,
    release_filename,
)

Ncm.cfg_init()

needs_data = pytest.mark.skipif(
    find_baseline_file(COMMANDER_RELPATH) is None,
    reason="Planck baseline (plc_3.0) clik data not found",
)

_EXPECTED_TYPE = {
    PlanckReleaseId.PR3_COMMANDER: "DataPlanckCommander",
    PlanckReleaseId.PR3_SIMALL_EE: "DataPlanckSimall",
    PlanckReleaseId.PR3_SIMALL_BB: "DataPlanckSimall",
    PlanckReleaseId.PR3_SIMALL_EEBB: "DataPlanckSimall",
    PlanckReleaseId.PR3_PLIK_TT: "DataPlanckSmica",
    PlanckReleaseId.PR3_PLIK_TTTEEE: "DataPlanckSmica",
    PlanckReleaseId.PR3_PLIK_LITE_TT: "DataPlanckPlikLite",
    PlanckReleaseId.PR3_PLIK_LITE_TTTEEE: "DataPlanckPlikLite",
    PlanckReleaseId.PR3_LENSING: "DataPlanckLensing",
    PlanckReleaseId.PR3_LENSING_MARGED: "DataPlanckLensing",
}


@pytest.mark.planck_data
@needs_data
def test_build_release_and_load(tmp_path):
    """build_release writes serialized objects that reload to the right types."""
    out = str(tmp_path)
    written = build_release(out_dir=out)
    assert written

    for rid in PlanckReleaseId:
        path = os.path.join(out, release_filename(rid))
        if not os.path.exists(path):
            continue  # source clik data for this id was missing
        data = load_planck_release(rid, cache_dir=out)
        assert data.__class__.__name__ == _EXPECTED_TYPE[rid]


@pytest.mark.planck_data
@needs_data
def test_release_block_evaluates(tmp_path):
    """A block loaded from the release self-configures its CBE and evaluates."""
    # pylint: disable=import-outside-toplevel
    from numcosmo_py.cosmology import create_cosmo, HIPrimModel
    from numcosmo_py.experiments.planck18 import mset_set_parameters, Planck18Types

    out = str(tmp_path)
    build_release(out_dir=out)

    cbe = Nc.HIPertBoltzmannCBE.new()  # bare: the loaded block configures it
    cosmo = create_cosmo(prim_model=HIPrimModel.POWER_LAW)
    planck = Nc.PlanckFICorTT()
    planck.params_set_default_ftype()
    mset = Ncm.MSet.new_array([planck, cosmo])
    mset_set_parameters(mset, Planck18Types.TT, HIPrimModel.POWER_LAW)
    mset.prepare_fparam_map()

    data = load_planck_release(PlanckReleaseId.PR3_COMMANDER, cbe, cache_dir=out)
    assert isinstance(data, Nc.DataPlanckCommander)
    data.prepare(mset)
    assert np.isfinite(data.m2lnL_val(mset))


# -----------------------------------------------------------------------------
# Synthetic-data tests: no plc_3.0 tree needed, so they also run where the real
# Planck data is absent (CI included). See tests/python/fixtures_planck.py.
# -----------------------------------------------------------------------------


def test_release_filenames_are_distinct():
    """Every curated id maps to its own release/cache filename."""
    names = {release_filename(rid) for rid in PlanckReleaseId}
    assert len(names) == len(list(PlanckReleaseId))
    assert all(n.startswith("planck_native_") and n.endswith(".gvar") for n in names)


def test_release_url_is_the_shared_datafile_release():
    """The artifacts are served from the shared NumCosmo data-asset release."""
    assert pnr.RELEASE_TAG == "datafile-release-v1.0.0"
    assert pnr.RELEASE_URL.endswith(pnr.RELEASE_TAG)


def _write_cached(tmp_path, rid, data):
    """Serialize @data where load_planck_release expects @rid's cache entry."""
    ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
    ser.to_binfile(data, os.path.join(str(tmp_path), release_filename(rid)))


def test_load_from_cache_without_download(tmp_path, monkeypatch):
    """A cached artifact is deserialized in place, with no network access."""
    clik = make_commander_cldf(tmp_path)
    _write_cached(tmp_path, PlanckReleaseId.PR3_COMMANDER, build_commander(clik))

    def _no_network(*args, **kwargs):
        raise AssertionError("must not download a cached artifact")

    monkeypatch.setattr(pnr.urllib.request, "urlretrieve", _no_network)

    cbe = Nc.HIPertBoltzmannCBE.new()
    mset, _ = planck_mset()
    data = load_planck_release(
        PlanckReleaseId.PR3_COMMANDER, cbe, cache_dir=str(tmp_path)
    )

    assert isinstance(data, Nc.DataPlanckCommander)
    assert data.peek_hipert_boltzmann() is not None

    data.prepare(mset)
    assert np.isfinite(data.m2lnL_val(mset))


def test_download_populates_the_cache(tmp_path, monkeypatch):
    """A missing artifact is fetched once and reused from the cache afterwards."""
    source = tmp_path / "source"
    cache = tmp_path / "cache"
    source.mkdir()
    clik = make_commander_cldf(source)
    _write_cached(source, PlanckReleaseId.PR3_COMMANDER, build_commander(clik))

    calls = []

    def _fake_download(url, path):
        calls.append(url)
        shutil.copyfile(
            os.path.join(str(source), release_filename(PlanckReleaseId.PR3_COMMANDER)),
            path,
        )

    monkeypatch.setattr(pnr.urllib.request, "urlretrieve", _fake_download)

    for _ in range(2):
        data = load_planck_release(PlanckReleaseId.PR3_COMMANDER, cache_dir=str(cache))
        assert isinstance(data, Nc.DataPlanckCommander)

    assert len(calls) == 1
    assert calls[0].startswith(pnr.RELEASE_URL)


def test_download_failure_is_reported_and_leaves_no_stub(tmp_path, monkeypatch):
    """A failed download raises with the URL and removes the partial file."""

    def _fail(url, path):
        with open(path, "wb") as fh:  # a partial file, as urlretrieve may leave
            fh.write(b"partial")
        raise OSError("network down")

    monkeypatch.setattr(pnr.urllib.request, "urlretrieve", _fail)

    with pytest.raises(RuntimeError, match="cannot download"):
        load_planck_release(PlanckReleaseId.PR3_LENSING, cache_dir=str(tmp_path))

    assert not os.path.exists(
        os.path.join(str(tmp_path), release_filename(PlanckReleaseId.PR3_LENSING))
    )


def test_build_release_skips_ids_without_source_data(tmp_path, monkeypatch):
    """build_release writes the ids it can build and silently skips the rest."""
    clik = make_commander_cldf(tmp_path)

    monkeypatch.setattr(
        pnr,
        "find_baseline_file",
        lambda relpath: clik if relpath == COMMANDER_RELPATH else None,
    )

    out = tmp_path / "release"
    written = build_release(out_dir=str(out))

    assert [os.path.basename(p) for p in written] == [
        release_filename(PlanckReleaseId.PR3_COMMANDER)
    ]
    assert (out / release_filename(PlanckReleaseId.PR3_COMMANDER)).exists()


def test_build_release_finds_nothing_without_data(tmp_path, monkeypatch):
    """With no clik tree at all, build_release writes nothing."""
    monkeypatch.setattr(pnr, "find_baseline_file", lambda relpath: None)

    assert build_release(out_dir=str(tmp_path / "empty")) == []


_REGISTRY_CHILD = """
import sys
from numcosmo_py import Ncm
from gi.repository import GObject

Ncm.cfg_init()

missing = []
for name in sys.argv[1:]:
    try:
        GObject.type_from_name(name)
    except RuntimeError:
        missing.append(name)
print(",".join(missing))
"""

_NATIVE_TYPE_NAMES = [
    "NcDataPlanckCommander",
    "NcDataPlanckLensing",
    "NcDataPlanckPlikLite",
    "NcDataPlanckSimall",
    "NcDataPlanckSmica",
]


def test_native_types_are_registered_in_a_fresh_process():
    """The native Planck types must resolve by name after cfg_init() alone.

    Loading a release artifact deserializes it in a process that never built one,
    and ncm_serialize_from_name_params() resolves the stored name with
    g_type_from_name(), which only succeeds once the GType has been realized --
    by ncm_cfg_register_objects() at startup, or incidentally by touching the
    class. Every other test here builds the object first, so they realize the
    type themselves and pass even when it is unregistered; that is exactly how
    load_planck_release() came to abort with "object `NcDataPlanckCommander' is
    not registered". The child must therefore construct nothing.
    """
    # pylint: disable=import-outside-toplevel
    import subprocess  # nosec B404 - fixed argv, no shell
    import sys

    proc = subprocess.run(  # nosec B603 - fixed argv, no shell
        [sys.executable, "-c", _REGISTRY_CHILD, *_NATIVE_TYPE_NAMES],
        capture_output=True,
        text=True,
        check=False,
        timeout=600,
    )

    assert proc.returncode == 0, proc.stderr
    missing = [n for n in proc.stdout.strip().split(",") if n]
    assert not missing, (
        f"not passed to ncm_cfg_register_obj(): {missing} -- a release artifact "
        f"of these types cannot be deserialized in a fresh process"
    )
