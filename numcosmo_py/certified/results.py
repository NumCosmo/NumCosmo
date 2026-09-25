#
# results.py
#
# Thu Sep 24 2026
# Copyright  2026  Sandro Dias Pinto Vitenti
# <vitenti@uel.br>
#
# results.py
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

"""Results measured against certified tables, one file per campaign.

A campaign is one run of one script against one certified table: a code at a stated
version and configuration, on one machine. Its file has the layout of a truth table
(:mod:`numcosmo_py.certified.tables`): the campaign's metadata as entries, and one
:class:`~gi.repository.NumCosmoMath.Catalog` of measurements per window or pair under
``cells/<name>``, so :func:`numcosmo_py.certified.tables.read` and
:func:`~numcosmo_py.certified.tables.groups` read it.

The metadata says what was compared with what: the reference table and its SHA-256, the
code and its version, the configuration, and the environment captured at run time by
:func:`environment`. Every row carries the reference value and radius, the measured
value, and the ``scale`` its deviation is measured against, which the caller states (the
multipole's peak for a radial integral, the pair's scale for C_ell); ``dev_scaled`` and
``dev_rel`` are computed here from those, never passed in.
"""

import ctypes
import datetime
import hashlib
import json
import os
import pathlib
import platform
import re
import site
import socket
import subprocess
import sys
import typing

import numpy as np
from astropy.table import Table

from numcosmo_py import Ncm, dict_to_var_dict
from numcosmo_py.catalog.table import catalog_from_table
from numcosmo_py.certified import tables

THREAD_VARIABLES: typing.Final[tuple[str, ...]] = (
    "OMP_NUM_THREADS",
    "OMP_THREAD_LIMIT",
    "OMP_MAX_ACTIVE_LEVELS",
    "OPENBLAS_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "MKL_NUM_THREADS",
)

_LIBRARIES = re.compile(
    r"(/\S*(?:numcosmo|gsl|blas|lapack|flint|gomp|libomp|fftw|pylevin|ccl)\S*\.so\S*)"
)

REQUIRED_COLUMNS: typing.Final[tuple[str, ...]] = (
    "reference",
    "radius",
    "measured",
    "scale",
)


def _cpu_model() -> str:
    try:
        for line in (
            pathlib.Path("/proc/cpuinfo").read_text(encoding="utf-8").splitlines()
        ):
            if line.startswith("model name"):
                return line.split(":", 1)[1].strip()
    except OSError:
        pass

    return platform.processor()


def _mapped_libraries() -> list[str]:
    """Return the shared libraries this process has mapped, not what ldd predicts."""
    try:
        maps = pathlib.Path("/proc/self/maps").read_text(encoding="utf-8")
    except OSError:
        return []

    return sorted(set(_LIBRARIES.findall(maps)))


def _openblas_core(libraries: list[str]) -> str:
    """Return the core the loaded OpenBLAS selected, or an empty string."""
    for lib in libraries:
        if "openblas" in lib:
            try:
                handle = ctypes.CDLL(lib)
                handle.openblas_get_corename.restype = ctypes.c_char_p
                return handle.openblas_get_corename().decode()
            except (OSError, AttributeError):
                return ""

    return ""


def _git(tree: pathlib.Path, *args: str) -> str:
    try:
        return subprocess.run(
            ["git", "-C", str(tree), *args],
            capture_output=True,
            text=True,
            check=True,
            timeout=10,
        ).stdout.strip()
    except (OSError, subprocess.SubprocessError):
        return ""


def _library_build(libraries: list[str], source: pathlib.Path) -> dict[str, typing.Any]:
    """Describe the libnumcosmo this process loaded, and how stale it may be.

    The source tree's commit says nothing about a binary built earlier; the library's
    own modification time, its build directory's compiler and options, and the date of
    the tree's last change to the C sources do.
    """
    out: dict[str, typing.Any] = {}
    lib = next((path for path in libraries if "libnumcosmo" in path), "")

    if not lib:
        return out

    path = pathlib.Path(lib)
    out["env_libnumcosmo"] = lib
    out["env_libnumcosmo_mtime"] = datetime.datetime.fromtimestamp(
        path.stat().st_mtime, datetime.timezone.utc
    ).isoformat()
    out["env_numcosmo_last_c_commit_date"] = _git(
        source, "log", "-1", "--format=%cI", "--", "numcosmo"
    )

    for build in path.parents:
        info = build / "meson-info"

        if info.is_dir():
            options = {
                o["name"]: o["value"]
                for o in json.loads(
                    (info / "intro-buildoptions.json").read_text(encoding="utf-8")
                )
            }
            compilers = json.loads(
                (info / "intro-compilers.json").read_text(encoding="utf-8")
            )
            c_comp = compilers.get("host", {}).get("c", {})
            out["env_libnumcosmo_build_dir"] = str(build)
            out["env_libnumcosmo_buildtype"] = str(options.get("buildtype", ""))
            out["env_libnumcosmo_c_args"] = " ".join(options.get("c_args", []))
            out["env_libnumcosmo_compiler"] = (
                f"{c_comp.get('id', '')} {c_comp.get('version', '')}".strip()
            )
            break

    return out


PYTHON_MODULES: typing.Final[tuple[str, ...]] = (
    "numpy",
    "scipy",
    "astropy",
    "gi",
    "pylevin",
    "pyccl",
)


def _python_modules() -> dict[str, typing.Any]:
    """Record where the modules a campaign depends on were imported from.

    A per-user site directory (``~/.local/lib/pythonX.Y``) is searched ahead of an
    environment's own site-packages and can shadow it with an editable install; the
    file of each module says which one was used.
    """
    out: dict[str, typing.Any] = {
        "env_user_site_enabled": bool(site.ENABLE_USER_SITE),
        "env_user_site": site.getusersitepackages() if site.ENABLE_USER_SITE else "",
    }

    for name in PYTHON_MODULES:
        module = sys.modules.get(name)

        if module is not None:
            out[f"env_module_{name}"] = str(getattr(module, "__file__", ""))
            out[f"env_module_{name}_version"] = str(getattr(module, "__version__", ""))

    return out


def environment() -> dict[str, typing.Any]:
    """Capture the environment of this process, for a campaign's metadata.

    The NumCosmo commit comes from the source tree ``numcosmo_py`` was imported from,
    with whether it had uncommitted changes, or from the provenance a conda env built
    by ``setup_numcosmo_env.sh`` records. The libraries are those mapped at run time.

    :return: a flat dictionary of strings, numbers and string lists.
    """
    libraries = _mapped_libraries()
    source = pathlib.Path(__file__).resolve().parents[2]
    env: dict[str, typing.Any] = {
        "env_date": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "env_host": socket.gethostname(),
        "env_cpu": _cpu_model(),
        "env_python": sys.version.split()[0],
        "env_byte_order": sys.byteorder,
        "env_openblas_core": _openblas_core(libraries),
        "env_libraries": ";".join(libraries),
        "env_numcosmo_py": str(pathlib.Path(__file__).resolve().parents[1]),
        "env_numcosmo_commit": _git(source, "rev-parse", "HEAD"),
        "env_numcosmo_dirty": bool(
            _git(source, "status", "--porcelain", "--untracked-files=no")
        ),
    }

    env.update(_library_build(libraries, source))
    env.update(_python_modules())

    for name in THREAD_VARIABLES:
        env[f"env_{name}"] = os.environ.get(name, "")

    prefix = os.environ.get("CONDA_PREFIX", "")
    summary = pathlib.Path(prefix, "provenance", "summary.txt")

    if prefix and summary.is_file():
        env["env_conda_prefix"] = prefix
        env["env_conda_provenance"] = summary.read_text(encoding="utf-8")

    script = pathlib.Path(sys.argv[0]) if sys.argv and sys.argv[0] else None

    if script is not None and script.is_file():
        env["env_script"] = str(script.resolve())
        env["env_script_sha256"] = hashlib.sha256(script.read_bytes()).hexdigest()

    return env


class Campaign:
    """The measurements of one code against one certified table.

    :param reference: path of the certified table (a binfile of
        :mod:`numcosmo_py.certified.tables`) the measurements are against.
    :param code: what was measured, e.g. ``"numcosmo"`` or ``"pylevin"``.
    :param code_version: its version or commit.
    :param config: settings shared by every row, as plain values.
    """

    def __init__(
        self,
        reference: pathlib.Path,
        code: str,
        code_version: str,
        config: dict[str, typing.Any] | None = None,
    ) -> None:
        """Record the reference table's identity and capture the environment."""
        ref_meta = tables.metadata(tables.read(reference))
        self.meta: dict[str, typing.Any] = {
            "kind": "results",
            "reference_table": reference.name,
            "reference_sha256": hashlib.sha256(reference.read_bytes()).hexdigest(),
            "reference_kind": ref_meta["kind"],
            "reference_source_sha256": ref_meta["source_sha256"],
            "code": code,
            "code_version": code_version,
        }

        for key, value in (config or {}).items():
            self.meta[f"config_{key}"] = value

        self.meta.update(environment())
        self.groups: dict[str, Table] = {}

    def add(self, group: str, columns: dict[str, typing.Any]) -> None:
        """Add the measurements of one window or pair.

        :param group: the group's name in the reference table.
        :param columns: equal-length columns; ``reference``, ``radius``, ``measured``
            and ``scale`` are required, anything else (``ell``, ``k``, per-row
            settings) is stored as given.
        :raises ValueError: if a required column is missing or lengths differ.
        """
        missing = [c for c in REQUIRED_COLUMNS if c not in columns]

        if missing:
            raise ValueError(f"{group}: missing columns {missing}")

        table = Table({name: np.asarray(values) for name, values in columns.items()})
        diff = np.abs(table["measured"] - table["reference"])
        table["dev_scaled"] = diff / table["scale"]

        with np.errstate(divide="ignore", invalid="ignore"):
            table["dev_rel"] = np.where(
                table["reference"] != 0.0, diff / np.abs(table["reference"]), np.inf
            )

        self.groups[group] = table

    def finish(self) -> Ncm.VarDict:
        """Return the campaign as a VarDict."""
        ser = Ncm.Serialize.new(Ncm.SerializeOpt.CLEAN_DUP)
        out = dict_to_var_dict(self.meta)

        for name, table in self.groups.items():
            out.set_object(
                f"{tables.CELLS_PREFIX}{name}", ser, catalog_from_table(table)
            )

        return out

    def write(self, path: pathlib.Path) -> None:
        """Write the campaign as a NumCosmo binfile."""
        tables.write(self.finish(), path)
