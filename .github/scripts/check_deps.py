#!/usr/bin/env python3
"""Check that every third-party Python module NumCosmo imports is declared.

Three lists declare the Python dependencies, and this script compares each one
with the import statements in the tracked sources:

* `pyproject.toml`: the core `dependencies` must cover the third-party modules
  that `import numcosmo_py` loads, and each extra in ENTRY_POINTS must cover the
  modules its entry package loads (the `app` extra covers the `numcosmo` CLI).
  Every other unguarded module-level import in `numcosmo_py` must be in some
  extra.
* `environment.yml`: every third-party module imported anywhere in the sources
  read below (package, tests, docs, examples, tools, CI scripts), and every
  `pyproject.toml` entry except those of the `firecrown` extra.
* the conda-forge recipe (`feedstock` subcommand, needs network): its `run`
  requirements must cover the core, the extras in ENTRY_POINTS and the imports
  of `examples/`, which meson installs.

Imports are read statically with `ast`. An import inside a function or a `try`
block, or requested through `pytest.importorskip`, counts as used by the
environment.yml check but not as loaded by an entry package. `notebooks/` is
not read: no build, test or documentation page uses it.

Usage:
    python .github/scripts/check_deps.py repo
    python .github/scripts/check_deps.py feedstock
"""

import argparse
import ast
import functools
import re
import subprocess
import sys
import tomllib
import urllib.request
import warnings
from collections.abc import Callable, Iterator
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
ENV_FILE = ROOT / "environment.yml"
PYPROJECT = ROOT / "pyproject.toml"
FEEDSTOCK_RECIPE = (
    "https://raw.githubusercontent.com/conda-forge/numcosmo-feedstock/main/"
    "recipe/meta.yaml"
)

# Entry package of numcosmo_py -> the pyproject extra that must cover the
# third-party modules it loads, on top of the core dependencies.
ENTRY_POINTS = {
    "numcosmo_py.app": "app",
    "numcosmo_py.catalog": "catalog",
    "numcosmo_py.plotting": "plot",
}

# Extras installed only on request; the environment.yml check skips them.
OPTIONAL_EXTRAS = {"firecrown"}

# Top-level directories whose Python sources are not checked.
SKIPPED_DIRS = {"notebooks"}

# Import name -> distribution name, where the two differ.
IMPORT_TO_DIST = {
    "IPython": "ipython",
    "gi": "pygobject",
    "mpl_toolkits": "matplotlib",
    # pydantic requires one exact pydantic-core version.
    "pydantic_core": "pydantic",
    "pytest_cov": "pytest-cov",
    "pytest_lazy_fixtures": "pytest-lazy-fixtures",
    "pytest_timeout": "pytest-timeout",
    "xdist": "pytest-xdist",
    "yaml": "pyyaml",
}

# PyPI name -> conda-forge name, where the two differ.
PIP_TO_CONDA = {
    "matplotlib": "matplotlib-base",
}

# Distribution -> conda-forge packages that also provide it in the recipe. The
# recipe uses astropy-base: astropy >=7 pulls h5py, which lags hdf5 migrations.
RECIPE_ALTERNATIVES = {
    "astropy": {"astropy-base"},
}

STDLIB = set(sys.stdlib_module_names) | {"__future__"}
REQUIREMENT_NAME = re.compile(r"^\s*([A-Za-z0-9_.\-]+)")
ENV_ENTRY = re.compile(r"^\s+-\s+([A-Za-z0-9_.\-]+)")


def main() -> None:
    """Run the selected check and exit non-zero when it finds a problem."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("check", choices=["repo", "feedstock"])
    args = parser.parse_args()

    errors = check_repo() if args.check == "repo" else check_feedstock()
    for error in errors:
        print(f"error: {error}")
    if errors:
        sys.exit(1)
    print(f"python dependencies declared ({args.check})")


def check_repo() -> list[str]:
    """Compare pyproject.toml and environment.yml with the imports."""
    groups = pyproject_groups()
    core = groups["core"]
    sources = read_sources()
    package = {path: tree for path, tree in sources.items() if is_package(path)}
    errors = []

    reach = package_reach(package)
    declared_anywhere = set().union(*groups.values())
    for path, tree in package.items():
        for name in top_level_imports(tree) - declared_anywhere:
            errors.append(f"{path} imports {name}, which no pyproject list declares")

    for entry, extra in [("numcosmo_py", None), *ENTRY_POINTS.items()]:
        allowed = core | (groups[extra] if extra else set())
        where = f"the {extra} extra" if extra else "the core dependencies"
        for name, path in sorted(reach(entry).items()):
            if name not in allowed:
                errors.append(f"{entry} loads {name} ({path}), missing from {where}")

    env = environment_packages()
    optional = set().union(*(groups[g] for g in OPTIONAL_EXTRAS))
    required = set().union(*(groups[g] for g in groups if g not in OPTIONAL_EXTRAS))
    for name in sorted(required):
        if conda_name(name) not in env:
            errors.append(f"pyproject declares {name}, missing from environment.yml")
    for path, tree in sorted(sources.items()):
        for name in sorted(all_imports(tree) - optional):
            if conda_name(name) not in env:
                errors.append(f"{path} imports {name}, missing from environment.yml")
    return errors


def check_feedstock() -> list[str]:
    """Compare the conda-forge recipe's run requirements with pyproject.toml."""
    groups = pyproject_groups()
    with urllib.request.urlopen(FEEDSTOCK_RECIPE, timeout=30) as response:
        run = recipe_run_packages(response.read().decode())
    needed = {name: "pyproject core" for name in groups["core"]}
    for extra in ENTRY_POINTS.values():
        needed.update({name: f"pyproject {extra} extra" for name in groups[extra]})
    for path, tree in read_sources().items():
        if path.startswith("examples/"):
            needed.update({name: path for name in all_imports(tree)})
    return [
        f"recipe run requirements lack {conda_name(name)} (needed by {source})"
        for name, source in sorted(needed.items())
        if not ({conda_name(name)} | RECIPE_ALTERNATIVES.get(name, set())) & run
    ]


def pyproject_groups() -> dict[str, set[str]]:
    """Return the core dependencies and each extra as sets of distribution names.

    An extra listing `numcosmo_py[...]` is expanded to the extras it names.
    """
    project = tomllib.loads(PYPROJECT.read_text())["project"]
    raw = {"core": project["dependencies"], **project["optional-dependencies"]}
    groups: dict[str, set[str]] = {}

    def expand(group: str) -> set[str]:
        if group not in groups:
            names = set()
            for req in raw[group]:
                if req.startswith("numcosmo_py["):
                    for sub in req[len("numcosmo_py[") : -1].split(","):
                        names |= expand(sub.strip())
                elif match := REQUIREMENT_NAME.match(req):
                    names.add(normalize(match.group(1)))
            groups[group] = names
        return groups[group]

    for group in raw:
        expand(group)
    return groups


def environment_packages() -> set[str]:
    """Return the conda package names listed in environment.yml."""
    names = set()
    in_deps = False
    for line in ENV_FILE.read_text().splitlines():
        if line.startswith("dependencies:"):
            in_deps = True
        elif in_deps and (match := ENV_ENTRY.match(line)):
            names.add(normalize(match.group(1)))
    return names


def recipe_run_packages(text: str) -> set[str]:
    """Return the package names of the `run:` section of a meta.yaml recipe."""
    names = set()
    in_run = False
    for line in text.splitlines():
        if re.match(r"^  run:", line):
            in_run = True
        elif in_run and re.match(r"^\s{0,2}\S", line):
            break
        elif in_run and (match := ENV_ENTRY.match(line)):
            names.add(normalize(match.group(1)))
    return names


def read_sources() -> dict[str, ast.Module]:
    """Parse every tracked Python source, Python script and Quarto Python cell."""
    files = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    sources = {}
    for path in files:
        is_script = "." not in Path(path).name
        if path.split("/")[0] in SKIPPED_DIRS or not (
            path.endswith((".py", ".qmd")) or is_script
        ):
            continue
        text = (ROOT / path).read_text(encoding="utf-8", errors="replace")
        if path.endswith(".qmd"):
            cells = re.findall(r"```\{python[^}]*\}\n(.*?)```", text, re.DOTALL)
            text = "\n".join(cells)
        elif is_script and not re.match(r"#!.*python", text):
            continue
        with warnings.catch_warnings(action="ignore", category=SyntaxWarning):
            sources[path] = ast.parse(text, filename=path)
    return sources


def is_package(path: str) -> bool:
    """Return whether path is a module of the installed numcosmo_py package."""
    return path.startswith("numcosmo_py/") and path.endswith(".py")


def all_imports(tree: ast.Module) -> set[str]:
    """Return the third-party distributions a source imports anywhere."""
    names = set()
    for node in ast.walk(tree):
        if (
            isinstance(node, ast.Call)
            and getattr(node.func, "attr", "") == "importorskip"
            and node.args
            and isinstance(node.args[0], ast.Constant)
            and isinstance(node.args[0].value, str)
        ):
            names.add(node.args[0].value)
        else:
            names.update(imported_modules(node))
    return {d for d in map(distribution, names) if d}


def top_level_imports(tree: ast.Module) -> set[str]:
    """Return the third-party distributions imported unguarded at module level."""
    names = set()
    for node in module_level(tree):
        names.update(imported_modules(node))
    return {d for d in map(distribution, names) if d}


def module_level(tree: ast.Module) -> Iterator[ast.stmt]:
    """Yield the statements a module runs on import, skipping `try` blocks."""
    stack = list(tree.body)
    while stack:
        node = stack.pop()
        yield node
        if isinstance(node, ast.If):
            if not is_type_checking(node.test):
                stack.extend(node.body)
            stack.extend(node.orelse)
        elif isinstance(node, ast.With):
            stack.extend(node.body)


def is_type_checking(test: ast.expr) -> bool:
    """Return whether an if-test is `TYPE_CHECKING` or `typing.TYPE_CHECKING`."""
    name = getattr(test, "id", None) or getattr(test, "attr", None)
    return name == "TYPE_CHECKING"


def imported_modules(node: ast.AST) -> list[str]:
    """Return the absolute module names an import statement names."""
    if isinstance(node, ast.Import):
        return [alias.name for alias in node.names]
    if isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
        return [node.module]
    return []


def distribution(module: str) -> str | None:
    """Return the distribution providing module, or None for stdlib and local."""
    top = module.split(".")[0]
    if top in STDLIB or top in local_modules():
        return None
    return normalize(IMPORT_TO_DIST.get(top, top))


def package_reach(
    package: dict[str, ast.Module],
) -> Callable[[str], dict[str, str]]:
    """Return a function mapping an entry module to the distributions it loads.

    The function follows unguarded module-level imports between numcosmo_py
    modules, including the parent packages of each module, and returns each
    third-party distribution reached with one module that imports it.
    """
    modules = {module_name(path): (path, tree) for path, tree in package.items()}

    def edges(name: str) -> set[str]:
        path, tree = modules[name]
        base = name if path.endswith("__init__.py") else name.rpartition(".")[0]
        out = set()
        for node in module_level(tree):
            if isinstance(node, ast.Import):
                targets = [alias.name for alias in node.names]
            elif isinstance(node, ast.ImportFrom):
                if node.level:
                    parent = base.rsplit(".", node.level - 1)[0]
                    target = f"{parent}.{node.module}" if node.module else parent
                else:
                    target = node.module or ""
                targets = [target] + [f"{target}.{a.name}" for a in node.names]
            else:
                continue
            for target in targets:
                if target.split(".")[0] == "numcosmo_py":
                    out.update(with_parents(target))
        return {m for m in out if m in modules}

    def reach(entry: str) -> dict[str, str]:
        seen: set[str] = set()
        stack = with_parents(entry)
        found: dict[str, str] = {}
        while stack:
            name = stack.pop()
            if name in seen or name not in modules:
                continue
            seen.add(name)
            path, tree = modules[name]
            for dist in sorted(top_level_imports(tree)):
                found.setdefault(dist, path)
            stack.extend(edges(name))
        return found

    return reach


def module_name(path: str) -> str:
    """Return the dotted module name of a numcosmo_py source path."""
    name = path.removesuffix(".py").replace("/", ".")
    return name.removesuffix(".__init__")


def with_parents(name: str) -> list[str]:
    """Return name and every parent package of it."""
    parts = name.split(".")
    return [".".join(parts[: i + 1]) for i in range(len(parts))]


def conda_name(name: str) -> str:
    """Return the conda-forge package name of a distribution."""
    return PIP_TO_CONDA.get(name, name)


def normalize(name: str) -> str:
    """Return a distribution name in lower case with `-` separators."""
    return re.sub(r"[-_.]+", "-", name).lower()


@functools.cache
def local_modules() -> set[str]:
    """Return the names a source can import from the repository itself.

    These are the stems of tracked Python files and the names of tracked
    directories, which test and documentation code import through sys.path.
    """
    files = subprocess.run(
        ["git", "ls-files"], cwd=ROOT, capture_output=True, text=True, check=True
    ).stdout.split()
    names = {"numcosmo_py", "numcosmo"}
    for path in files:
        parts = path.split("/")
        names.update(parts[:-1])
        if path.endswith(".py"):
            names.add(parts[-1].removesuffix(".py"))
    return names


if __name__ == "__main__":
    main()
