#!/usr/bin/env python3
"""Generate and validate the CI conda lock files (see CONTRIBUTING.md).

`generate` re-solves `environment.yml` with conda-lock, once per CI target,
and writes an explicit lock file (plain URL list, no solving needed to
install it) to .github/locks/.  The sha256 of the `environment.yml` it was
solved from is stored in a comment inside each lock file.

`check` compares that comment against the current `environment.yml` and fails
when they differ, i.e. when a dependency was added without regenerating the
locks.
"""

import argparse
import hashlib
import re
import shutil
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENV_FILE = ROOT / "environment.yml"
LOCK_DIR = ROOT / ".github" / "locks"
SHA_COMMENT = "# environment.yml sha256: "

# One entry per conda environment built by CI. Keep in sync with the
# build-miniforge/build-miniforge-coverage matrices of build_check.yml -- a
# missing combination is not fatal, the job just falls back to solving.
TARGETS = [
    ("linux-64", "3.13", "openmpi"),
    ("linux-64", "3.13", "mpich"),
    ("osx-arm64", "3.13", "mpich"),
]

# Installed next to the MPI package by both MPI legs of build_check.yml.
EXTRA_PACKAGES = ["libfabric-devel"]

PYTHON_RE = re.compile(r"^\s*-\s*python\s*[=<>!]", re.IGNORECASE)
MPI_RE = re.compile(r"^(?P<indent>\s*-\s*)(openmpi|mpich)\s*$", re.IGNORECASE)


def env_sha256() -> str:
    """Hex sha256 of environment.yml."""
    return hashlib.sha256(ENV_FILE.read_bytes()).hexdigest()


def lock_path(platform: str, python: str, mpi: str) -> Path:
    """Lock file of a target. The CI action builds this name too, in bash."""
    return LOCK_DIR / f"conda-{platform}-py{python}-{mpi}.lock"


def render_env_file(python: str, mpi: str) -> str:
    """environment.yml with python pinned, MPI swapped and extras appended."""
    out, found_python, found_mpi = [], False, False
    for line in ENV_FILE.read_text(encoding="utf-8").splitlines():
        mpi_match = MPI_RE.match(line)
        if PYTHON_RE.match(line):
            out.append(f"  - python ={python}")
            found_python = True
        elif mpi_match:
            out.append(f"{mpi_match.group('indent')}{mpi}")
            found_mpi = True
        else:
            out.append(line)
    if not (found_python and found_mpi):
        raise SystemExit(f"{ENV_FILE}: python and/or openmpi entry not found")
    return "\n".join(out + [f"  - {pkg}" for pkg in EXTRA_PACKAGES]) + "\n"


def solve_target(conda_lock: str, target: tuple[str, str, str]) -> Path:
    """Run conda-lock for one target and stamp the environment.yml hash in."""
    platform, python, mpi = target
    path = lock_path(platform, python, mpi)
    with tempfile.TemporaryDirectory() as tmp:
        spec = Path(tmp) / "environment.yml"
        spec.write_text(render_env_file(python, mpi), encoding="utf-8")
        # --filename-template must contain {platform}; we solve one at a time.
        template = str(path).replace(platform, "{platform}", 1)
        subprocess.run(
            [
                conda_lock,
                "lock",
                "--file",
                str(spec),
                "--platform",
                platform,
                "--kind",
                "explicit",
                "--without-cuda",
                "--filename-template",
                template,
            ],
            check=True,
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
        )
    # Record what was solved: `check` (and CI) compare this against
    # environment.yml. Comments are ignored by `conda create --file`.
    lines = path.read_text(encoding="utf-8").splitlines()
    lines.insert(1, f"{SHA_COMMENT}{env_sha256()}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def cmd_generate(args: argparse.Namespace) -> int:
    """Re-solve every selected target and rewrite its lock file."""
    conda_lock = shutil.which(args.conda_lock)
    if conda_lock is None:
        raise SystemExit(f"{args.conda_lock} not found (pip install conda-lock)")

    targets = [
        t
        for t in TARGETS
        if (not args.platform or t[0] in args.platform)
        and (not args.mpi or t[2] in args.mpi)
    ]
    if not targets:
        raise SystemExit("no target matches the given --platform/--mpi")

    LOCK_DIR.mkdir(parents=True, exist_ok=True)
    print(f"solving {len(targets)} target(s): " + ", ".join(t[0] for t in targets))
    # The solver itself is single-threaded (libsolv), so the only parallelism
    # available is one process per target.
    with ThreadPoolExecutor(max_workers=args.jobs or len(targets)) as pool:
        futures = {pool.submit(solve_target, conda_lock, t): t for t in targets}
        failed = 0
        for future in as_completed(futures):
            target = futures[future]
            try:
                print(f"wrote {future.result().relative_to(ROOT)}", flush=True)
            except subprocess.CalledProcessError as exc:
                print(f"FAILED {lock_path(*target).name}: {exc}", file=sys.stderr)
                failed += 1
    return 1 if failed else 0


def cmd_check(_args: argparse.Namespace) -> int:
    """Fail when a lock file is missing or was solved from a different env file."""
    expected = f"{SHA_COMMENT}{env_sha256()}"
    stale = []
    for target in TARGETS:
        path = lock_path(*target)
        if not path.is_file():
            stale.append(f"{path.relative_to(ROOT)}: missing")
        elif expected not in path.read_text(encoding="utf-8").splitlines():
            stale.append(
                f"{path.relative_to(ROOT)}: solved from a different environment.yml"
            )

    if stale:
        print("\n".join(f"error: {item}" for item in stale), file=sys.stderr)
        print(
            "\nRegenerate with `python .github/scripts/conda_locks.py generate`\n"
            "(needs conda-lock) or run the 'Update conda lock files' workflow,\n"
            "then commit the result. See CONTRIBUTING.md.",
            file=sys.stderr,
        )
        return 1
    print(f"conda lock files up to date ({len(TARGETS)} targets)")
    return 0


def main() -> int:
    """Parse arguments and run the requested subcommand."""
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(required=True)

    gen = sub.add_parser("generate", help="re-solve and rewrite the lock files")
    gen.add_argument("--platform", action="append", help="restrict to this platform")
    gen.add_argument("--mpi", action="append", help="restrict to this MPI")
    gen.add_argument("--jobs", type=int, default=0, help="concurrent solves")
    gen.add_argument("--conda-lock", default="conda-lock", help="conda-lock executable")
    gen.set_defaults(func=cmd_generate)

    chk = sub.add_parser("check", help="verify the locks match environment.yml")
    chk.set_defaults(func=cmd_check)

    args = parser.parse_args()
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
