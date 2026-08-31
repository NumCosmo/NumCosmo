#!/usr/bin/env python
#
# test_data_download.py
#
# Sun August 30 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# test_data_download.py
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

"""Tests for the shared data-file downloader.

The catalogs NumCosmo fetches are tens of megabytes, which is why this exercises
``jla_snls3_sdss_sys_stat.fits.sig`` instead: 310 bytes in the same release,
reached by the same code. That keeps the download path covered on every run --
including once a CI cache stops the big files from ever being fetched, which is
how these defects went unnoticed for months.

Each test runs in a subprocess with its own HOME, so the base directory really
is empty and a failure (fatal by design) can be observed rather than taking the
test session with it.
"""

import os
import subprocess
import sys
import time

import pytest

# 310 bytes, same release and same code path as the multi-megabyte catalogs.
TINY_ASSET = "jla_snls3_sdss_sys_stat.fits.sig"

_REPAIR = """
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()
Nc.DataPlanckLKL.download_baseline(Ncm.cfg_get_fullpath_base())
"""

_FETCH = """
import sys
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()
sys.stdout.write(Nc.DataSNIACov.get_fits(sys.argv[1], False))
"""

_FETCH_WL = """
import sys
from numcosmo_py import Nc, Ncm

Ncm.cfg_init()
sys.stdout.write(
    Nc.galaxy_wl_obs_catalog_id_get_filename(Nc.GalaxyWLObsCatalogId.HSC_PDR1_HWL16A_094)
)
"""

# Smallest curated weak-lensing catalog, 3.1 MB. Neither test below fetches it:
# each arranges for the file to be there, and asserts it was left alone.
WL_ASSET = "wl_obs_HWL16a-094.gvar"


def run_isolated(script: str, home, *args) -> subprocess.CompletedProcess:
    """Run @script with an empty HOME, so the data directory starts bare."""
    env = dict(os.environ, HOME=str(home))

    return subprocess.run(
        [sys.executable, "-c", script, *args],
        capture_output=True,
        text=True,
        check=False,
        env=env,
        timeout=900,
    )


def test_concurrent_download_is_safe(tmp_path):
    """Several processes fetching the same file must not corrupt it.

    This is the failure that took CI down: the destination is shared by every
    NumCosmo process, so without a lock two transfers wrote one path while a
    third read it. The reader saw a truncated file and the process aborted.
    """
    env = dict(os.environ, HOME=str(tmp_path))
    procs = [
        subprocess.Popen(  # pylint: disable=consider-using-with
            [sys.executable, "-c", _FETCH, TINY_ASSET],
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            env=env,
        )
        for _ in range(4)
    ]
    outs = [p.communicate(timeout=900) for p in procs]
    codes = [p.returncode for p in procs]

    assert codes == [0, 0, 0, 0], f"a worker died: {codes}\n{outs}"

    landed = list(tmp_path.rglob(TINY_ASSET))

    assert len(landed) == 1, f"expected exactly one copy, got {landed}"
    assert landed[0].stat().st_size == 310

    # A transfer in flight is invisible outside its own process, and nothing is
    # left behind if one dies mid-download.
    assert not list(tmp_path.rglob("*.part"))
    assert not list(tmp_path.rglob("*.lock"))


def test_a_failed_download_says_so(tmp_path):
    """wget's exit status must be read.

    It was discarded, so a 404 or a refused connection looked like success and
    the first symptom was a corrupt file read much later, with nothing pointing
    back at the download.
    """
    result = run_isolated(_FETCH, tmp_path, "this-asset-does-not-exist.fits")

    assert result.returncode != 0
    assert "wget failed" in result.stderr
    assert "this-asset-does-not-exist.fits" in result.stderr


def test_a_failed_download_leaves_nothing_behind(tmp_path):
    """A failure must not poison the cache for every later run.

    The transfer used to go straight to the final name, so whatever wget left
    there -- an error page, a truncated body -- was indistinguishable from real
    data on the next run.
    """
    run_isolated(_FETCH, tmp_path, "this-asset-does-not-exist.fits")

    assert not list(tmp_path.rglob("this-asset-does-not-exist.fits"))
    assert not list(tmp_path.rglob("*.part"))


def test_a_partial_tree_is_replaced(tmp_path):
    """An interrupted extraction must not look finished.

    A tarball unpacked half-way leaves its first files present and its last
    missing. Keying "already done" on any one of them -- or on the tree
    existing at all -- makes every later run skip the repair and fail on
    whichever file did not make it, which is how a run died on a plik
    component while commander sat right beside it. Only a marker written after
    the rename can mean complete.
    """
    commander = (
        tmp_path
        / ".numcosmo/baseline/plc_3.0/low_l/commander/commander_dx12_v3_2_29.clik"
    )
    commander.parent.mkdir(parents=True)
    commander.touch()

    result = run_isolated(_REPAIR, tmp_path)

    assert result.returncode == 0, result.stderr

    missing = (
        tmp_path / ".numcosmo/baseline/plc_3.0/hi_l/plik"
        "/plik_rd12_HM_v22b_TTTEEE.clik/clik/lkl_0/component_1/_mdb"
    )

    assert missing.exists(), "the partial tree was left unrepaired"
    assert (tmp_path / ".numcosmo/baseline/.numcosmo-complete").exists()
    assert not list(tmp_path.glob(".numcosmo/baseline.*"))


def test_a_waiter_uses_what_the_holder_produced(tmp_path):
    """The waiting branch must return the other process's file, not fetch again.

    The four-way race above asserts the outcome -- one copy, no debris -- but on
    a 310-byte asset the losers find the file already there and never reach the
    lock, so the waiting path never runs and the test would pass even if it were
    wrong. Here the lock is held from the outside and the data appears while the
    caller waits, which is the sequence that actually happens when one worker is
    mid-download and the others queue behind it.
    """
    base = tmp_path / ".numcosmo"
    base.mkdir()
    target = base / TINY_ASSET

    # Held by nobody, which is what a live download looks like from outside.
    (base / f"{TINY_ASSET}.lock").mkdir()

    proc = subprocess.Popen(  # pylint: disable=consider-using-with
        [sys.executable, "-c", _FETCH, TINY_ASSET],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=dict(os.environ, HOME=str(tmp_path)),
    )

    try:
        # Long enough that the child is past its own existence check and sitting
        # in the lock's wait loop -- at two seconds it was still starting up, saw
        # the file appear before it ever reached the lock, and returned through a
        # path this test was not written to exercise.
        time.sleep(6.0)
        started = time.monotonic()
        # The holder finishes: distinguishable from a real download by content.
        target.write_bytes(b"SENTINEL")
        out, err = proc.communicate(timeout=120)
        waited_after = time.monotonic() - started
    finally:
        if proc.poll() is None:  # pragma: no cover - only on a hang
            proc.kill()

    assert proc.returncode == 0, err
    assert out.strip() == str(target)

    # It noticed within a poll of the file appearing, rather than having
    # returned long before through its own existence check.
    assert waited_after < 30.0

    # 8, not 310: the waiter took the holder's file instead of downloading over
    # it -- which is the whole point of waiting.
    assert target.read_bytes() == b"SENTINEL"
    assert not list(tmp_path.rglob("*.part"))


def test_wl_catalog_waits_for_the_holder(tmp_path):
    """The WL catalogs must take the lock like every other download.

    They did not: the fetch ran `wget` straight to the final name, ignoring a
    held lock, so a concurrent reader could open the file half-written.

    Holding the lock from outside and letting the data appear separates the
    two: the old code writes the real catalog over the sentinel, the shared
    path waits and returns what the holder produced. The lock and rename
    machinery is raced four ways over the tiny asset above; what this checks is
    that this caller goes through it.
    """
    base = tmp_path / ".numcosmo"
    base.mkdir()
    target = base / WL_ASSET

    # Held by nobody, which is what a live download looks like from outside.
    (base / f"{WL_ASSET}.lock").mkdir()

    proc = subprocess.Popen(  # pylint: disable=consider-using-with
        [sys.executable, "-c", _FETCH_WL],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=dict(os.environ, HOME=str(tmp_path)),
    )

    try:
        # Same six seconds as the SNIa waiter above: long enough that the child
        # is past its own existence check and sitting in the lock's wait loop.
        time.sleep(6.0)
        # The holder finishes: distinguishable from a real fetch by content.
        target.write_bytes(b"SENTINEL")
        out, err = proc.communicate(timeout=120)
    finally:
        if proc.poll() is None:  # pragma: no cover - only on a hang
            proc.kill()

    assert proc.returncode == 0, err

    # The discriminator, asserted first so a regression names itself: 8 bytes,
    # not 3.1 MB, and no transfer announced. The old code ignored the lock and
    # wrote the real catalog over the holder's file.
    assert "Downloading" not in out, "fetched the catalog despite the held lock"
    assert target.read_bytes() == b"SENTINEL"

    # Last line, not the whole of stdout: a fetch would have printed first, and
    # that is the assertion above's job to report, not this one's.
    assert out.strip().splitlines()[-1] == str(target)
    assert not list(tmp_path.rglob("*.part"))


def test_wl_catalog_already_there_is_not_refetched(tmp_path):
    """A catalog already in the data directory is returned untouched."""
    base = tmp_path / ".numcosmo"
    base.mkdir()
    target = base / WL_ASSET
    target.write_bytes(b"ALREADY")

    result = run_isolated(_FETCH_WL, tmp_path)

    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == str(target)
    assert target.read_bytes() == b"ALREADY"
    assert not list(tmp_path.rglob("*.part"))
    assert not list(tmp_path.rglob("*.lock"))
