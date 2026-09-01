#!/usr/bin/env python3
#
# test_summary.py
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

"""Per-shard test timing report for the CI job summary.

The coverage job splits the C tests with meson's own `--slice K/N`, which
balances by test count. This report is how the consequence of that becomes
visible: if one test grows big enough to unbalance its shard, it shows up at
the top of the list. The remedy is to split that test, which is worth doing
for isolation and for per-test timeouts anyway -- not to reweight the split.

See the `build-miniforge-coverage` job in `.github/workflows/build_check.yml`.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path


def _read_testlog(testlog: Path) -> list[tuple[str, float]]:
    """Read a meson testlog.json (JSON-lines) into (name, duration) pairs.

    Returns an empty list if the file doesn't exist -- meson never writes it if
    the build itself failed before any test ran.
    """
    entries: list[tuple[str, float]] = []

    if testlog.exists():
        with testlog.open(encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                entry = json.loads(line)
                entries.append((entry["name"], float(entry.get("duration", 0.0))))

    return entries


def cmd_summary(args: argparse.Namespace) -> None:
    """Print a Markdown top-N slowest-durations report from a testlog.json."""
    rows = sorted(
        _read_testlog(Path(args.testlog)), key=lambda row: row[1], reverse=True
    )
    total = sum(duration for _, duration in rows)

    print(
        f"### Test durations — {args.title} "
        f"(top {args.top} of {len(rows)}, total {total:.1f}s)"
    )
    print("```")
    for name, duration in rows[: args.top]:
        print(f"{duration:8.2f}s  {name}")
    print("```")


def main() -> None:
    """Parse arguments and dispatch to the requested subcommand."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest="command", required=True)

    summary = sub.add_parser(
        "summary", help="print a Markdown slowest-durations report"
    )
    summary.add_argument("--testlog", required=True)
    summary.add_argument("--title", required=True)
    summary.add_argument("--top", type=int, default=25)
    summary.set_defaults(func=cmd_summary)

    args = parser.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
