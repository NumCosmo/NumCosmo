#!/usr/bin/env python3
#
# test_slicer.py
#
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

"""Duration-aware test slicing for the coverage job's C-tier shards.

meson's own `meson test --slice K/N` balances by test *count*, not wall
time, which lets one or two heavy tests dominate a shard's wall-clock
floor. This script replaces it: given historical per-test duration data
(cached from previous CI runs) and today's actual test list, it greedy
bin-packs the tests into N roughly-equal-wall-time slices.

Two subcommands:
  plan    -- print slice K's test names (one per line), computed from
             today's `meson test --list` output plus 0-or-more merged
             historical duration files. Tests with no historical duration
             get the median of known durations as an estimate; if there
             is no historical data at all, every test gets an equal
             weight, degrading gracefully to a count-balanced split.
  extract -- read a meson testlog.json and print a flat {name: duration}
             JSON object, for caching as tomorrow's historical duration
             data (see the `plan` subcommand's --durations).

See the `build-miniforge-coverage` job in
`.github/workflows/build_check.yml` for how these are wired into CI.
"""

from __future__ import annotations

import argparse
import json
import statistics
import subprocess
import sys
from pathlib import Path


def _load_durations(paths: list[Path]) -> dict[str, float]:
    """Merge {name: duration} dicts from possibly-missing JSON files."""
    merged: dict[str, float] = {}

    for path in paths:
        if not path.exists():
            continue
        with path.open(encoding="utf-8") as f:
            data = json.load(f)
        for name, duration in data.items():
            merged[name] = float(duration)

    return merged


def _current_tests(build_dir: Path, suites: list[str]) -> list[str]:
    """Return today's actual test names for the given suites, via meson."""
    cmd = ["meson", "test", "--list", "-C", str(build_dir)]
    for suite in suites:
        cmd += ["--suite", suite]
    result = subprocess.run(cmd, capture_output=True, text=True, check=True)

    return [line.strip() for line in result.stdout.splitlines() if line.strip()]


def _selector(name: str) -> str:
    """Convert a `--list`/testlog.json display name to a meson test selector.

    `meson test --list` and testlog.json's "name" field print e.g.
    "c-statistical+omp - numcosmo:ncm_stats_dist" -- the suite-tag prefix
    is for human reading only. The actual positional argument meson test
    accepts to select that specific test is "numcosmo:ncm_stats_dist"
    (the "subprojname:testname" part after the first " - ").
    """
    return name.split(" - ", 1)[-1]


def _lpt_partition(
    tests: list[str], weights: dict[str, float], num_slices: int
) -> list[list[str]]:
    """Greedy longest-processing-time-first bin packing.

    Sorts tests by weight descending, then drops each into whichever bin
    currently has the smallest running total -- within 4/3 of the optimal
    makespan in the worst case, and exact-balanced (as here) when weights
    are roughly uniform.
    """
    ordered = sorted(tests, key=lambda name: weights[name], reverse=True)
    bins: list[list[str]] = [[] for _ in range(num_slices)]
    totals = [0.0] * num_slices

    for name in ordered:
        i = min(range(num_slices), key=lambda i: totals[i])
        bins[i].append(name)
        totals[i] += weights[name]

    return bins


def cmd_plan(args: argparse.Namespace) -> None:
    """Print slice `args.slice`'s test names, one per line."""
    tests = _current_tests(Path(args.build_dir), args.suite)
    if not tests:
        raise SystemExit("no tests found for the given suite(s)")

    durations = _load_durations([Path(p) for p in args.durations])
    known = [durations[name] for name in tests if name in durations]
    fallback = statistics.median(known) if known else 1.0
    weights = {name: durations.get(name, fallback) for name in tests}

    bins = _lpt_partition(tests, weights, args.num_slices)

    # Safety net: the partition must be exhaustive and non-overlapping --
    # a bug here would silently drop or double-run tests in CI.
    flat = [name for b in bins for name in b]
    assert sorted(flat) == sorted(tests), "slice partition lost or duplicated tests"

    if not 1 <= args.slice <= args.num_slices:
        raise SystemExit(f"--slice must be between 1 and {args.num_slices}")

    for name in bins[args.slice - 1]:
        print(_selector(name))


def cmd_extract(args: argparse.Namespace) -> None:
    """Print a flat {name: duration} JSON object from a testlog.json."""
    testlog = Path(args.testlog)
    durations: dict[str, float] = {}

    if testlog.exists():
        with testlog.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                entry = json.loads(line)
                durations[entry["name"]] = float(entry.get("duration", 0.0))

    json.dump(durations, sys.stdout)
    print()


def main() -> None:
    """Parse arguments and dispatch to the requested subcommand."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    plan = sub.add_parser("plan", help="compute a duration-balanced test slice")
    plan.add_argument("--build-dir", default="build")
    plan.add_argument("--suite", action="append", required=True)
    plan.add_argument("--durations", action="append", default=[])
    plan.add_argument("--num-slices", type=int, required=True)
    plan.add_argument("--slice", type=int, required=True)
    plan.set_defaults(func=cmd_plan)

    extract = sub.add_parser(
        "extract", help="extract {name: duration} from a testlog.json"
    )
    extract.add_argument("--testlog", required=True)
    extract.set_defaults(func=cmd_extract)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
