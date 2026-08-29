# Contributing to NumCosmo

Thank you for your interest in contributing to NumCosmo. This guide covers the
mechanics of adding code, tests, and documentation so that a new contribution
integrates cleanly with the build, the test suite, and the generated API
reference.

For build and dependency setup, see the
[installation guide](https://numcosmo.readthedocs.io/en/latest/install.html).
NumCosmo uses the [meson](https://mesonbuild.com/) build system.

## Source layout

The library is split into two GObject namespaces, mirrored by the directory tree:

- `numcosmo/ncm/` — **NumCosmoMath**: foundation code with no cosmology
  (algebra, splines, integration, statistics, the model/MSet system, fitting).
- `numcosmo/nc/` — **NumCosmo**: cosmology models and likelihoods
  (background, perturbations, large-scale structure, clusters, CMB, …).

Each subdirectory groups a coherent family; an abstract base lives together with
its concrete subclasses. Vendored third-party code is quarantined under
`numcosmo/external/`.

## Adding a new class

1. **Place the files** in the subdirectory that matches the family, under
   `numcosmo/ncm/<area>/` or `numcosmo/nc/<area>/`.

2. **Register the sources** in `numcosmo/meson.build`. Add the `.c` to the
   `ncm_sources` (or `nc_sources`) `files()` list and the `.h` to `ncm_headers`
   (or `nc_headers`), using the subdirectory-qualified path, e.g.
   `'ncm/algebra/ncm_my_object.c'`.

3. **Export through the umbrella header.** Add the header to `numcosmo/numcosmo.h`
   (cosmology) or `numcosmo/numcosmo-math.h` (math). These two umbrellas are the
   only supported public include surface — never document reaching into a
   subdirectory header directly.

4. **Register the type** for serialization. Include the header in
   `numcosmo/ncm/core/ncm_cfg.c`, next to the other headers of its family, and
   register the object in `ncm_cfg_register_objects ()` with
   `ncm_cfg_register_obj (NCM_TYPE_MY_OBJECT);`. (`ncm_cfg.h` includes only what its
   own public API needs — do not add registration includes there.)

## Unit testing

Tests follow their sources as a hard invariant: a class in
`numcosmo/<ns>/<area>/` is tested in the mirrored `tests/c/<ns>/<area>/` and/or
`tests/python/<ns>/<area>/`. All new code must have at least minimal testing.

1. Add the test under the mirrored path, prefixed `test_`: `tests/c/<ns>/<area>/`
   for C, `tests/python/<ns>/<area>/` for Python.
2. **C tests must be registered** in `tests/c/meson.build`: add an entry giving the
   test `name` and its `sources`, since each C test is built as its own executable.
   Coverage of a class's ref/free/clear methods goes in the existing
   `tests/c/ncm/core/test_ncm_generic.c` rather than in a new file.
3. **Python tests are collected automatically** — the suites run `pytest` over
   `tests/python/`, so a new file needs no `meson.build` entry. Selection between
   lanes is by *marker*, not by path: a file with no marker joins the default lane,
   and one belonging to a dedicated shard (`xcor`, `powspec`, `app`, `omp`,
   `sphere_map`, `acceptance`) must declare a module-level
   `pytestmark = pytest.mark.<name>`.
4. Run the affected suite with `meson test -v -C build`, or a single file directly
   with `python -m pytest tests/python/<path> -q`.

## Documentation

Documentation is split by purpose, and the split is intentional:

- **API reference** (gi-docgen, generated from the GObject doc comments in the
  `.c` files) — brief and operational. State what the class does, the defining
  relation or signature, and the key methods. Keep it short.
- **Theoretical background** (the Quarto project under `docs/`, rendered to the
  website) — the physics, derivations, and full equations live here, in
  `docs/theory/<area>/<topic>.qmd`.

When a class involves non-trivial math, **do not put the derivation in the C doc
comment.** Put it on a theory page and link to it. The pattern is established by
`NcmCSQ1D`: see the short doc comment in
`numcosmo/ncm/dynamics/ncm_csq1d.c` and the corresponding
`docs/theory/csq1d.qmd`.

- From the C doc comment, link to the theory page with a plain anchor, e.g.
  `<a href="../../theory/csq1d.html">CSQ1D Formalism</a>`.
- From the theory page, link back to the API with wiki-style symbol references:
  `[[numcosmo-math|NcmCSQ1D]]`, `[[numcosmo|NcDistance]]`,
  `[[numcosmo-math|ncm_csq1d_prepare]]`. Unresolved references degrade to plain
  text, so verify the symbol name.

A single inline `$x$` in a doc comment is fine; a multi-line `\begin{align}`
derivation belongs on a theory page.

References use the bibliography only on the Quarto side: theory pages cite with
pandoc `[@Key]` against `docs/references.bib`; API doc comments (gi-docgen has no
bibliography support) use inline links instead, e.g.
`[Author (year)](https://arxiv.org/abs/...)`.

### Bibliography (`docs/references.bib`)

`bibtex-tidy` is the canonical formatter for `docs/references.bib`; the file is no
longer formatted by hand or by a reference manager. After adding or editing an
entry, run:

```bash
docs/tidy_references.sh          # re-format in place
docs/tidy_references.sh --check  # verify (this is what CI runs)
```

CI rejects a bibliography that is not tidy (`.github/workflows/bib_lint.yml`), so
run the formatter before committing. Install the pinned version with
`npm install -g bibtex-tidy@1.14.0` (different versions format differently). A
`pre-commit` hook is available — `pip install pre-commit && pre-commit install` —
to run this automatically. If you use a reference manager such as JabRef, treat
its output as a draft and re-tidy before committing; the local `file`, `owner`,
`timestamp`, and `__markedentry` fields are stripped from the repo copy.

### Adding a page

A new `.qmd` has to be registered in **four** places. Miss one and it either
fails to build, builds but is unreachable, or is reachable only from the
sidebar:

1. **`docs/<area>/meson.build`** — add the filename to the `files(...)` list, or
   it is never copied into the build tree. A brand-new subdirectory also needs
   its own `meson.build` and a `subdir('<name>')` line in the parent.
2. **`docs/_quarto.yml.in`, the render list** — the top-level list of pages to
   render (near `- tutorials/index.qmd`).
3. **`docs/_quarto.yml.in`, the sidebar** — the `contents:` entry under the
   right `section:`, with an `href:` and a short `text:` label.
4. **The index page for its area** — `docs/tutorials/index.qmd`,
   `docs/examples/index.qmd`, or `docs/theory/index.qmd`. The sidebar and the
   index are separate; adding a page to the sidebar alone leaves it off the
   landing page a reader actually browses.

Keep the section names in 2-4 consistent: if a page needs a new section, add it
to both the sidebar and the index.

Tutorial or example? A **tutorial** walks through a complete workflow and
explains every step, and is read start to finish; it goes in `docs/tutorials/`.
An **example** is a short self-contained script showing one part of the API in
use, read as reference; it goes in `docs/examples/`.

Match an existing page's front matter (`format: html` plus `ipynb`, and
`{{< include /_functions.qmd >}}`) so the notebook download link is generated.

### Building the documentation

Build the site with a documentation-enabled build directory. `BDocs` is the
conventional name, and the one CI uses:

```bash
meson setup BDocs -Ddocumentation=true
meson compile -C BDocs numcosmo-site
```

The rendered site lands in `BDocs/docs/numcosmo-site/`, with each executable
page as `.html`, `.ipynb` and `.py`. Build this before pushing a `.qmd`: the
page's code chunks run during the build, so a broken chunk is a build failure.

Do not invoke `quarto render` on a page directly. `_quarto.yml` is generated
from `_quarto.yml.in` by meson and is absent from a clean checkout, so the
`{{< include /_functions.qmd >}}` directive fails to resolve, and a hand-run
render also appends to `docs/.gitignore`.

Note: when you change **only** a `.c`/`.h` doc comment (no `.qmd`/`_quarto`
change), the site bundle may not re-render because the `numcosmo-site` target
does not detect changes in the generated API reference. Force it with:

```bash
ninja -C <builddir> -t clean docs/numcosmo-site
meson compile -C <builddir> numcosmo-site
```

This affects only local incremental builds; CI and Read the Docs build from
clean, so the published site is unaffected.

### The Read the Docs PR check

Read the Docs does not build the documentation — building it there was too slow.
The `docs-artifacts` GitHub Actions job builds the site and uploads it as
`numcosmo-site-<SHA>`, then restarts the RTD build, and RTD merely downloads
that artifact (`docs/download_docs_site.py`).

So `docs/readthedocs.org:numcosmo` **fails on every PR at first**: RTD's own
webhook fires on the push and looks for an artifact that does not exist yet. It
goes green once `docs-artifacts` finishes and restarts it, and RTD then serves a
per-PR preview. Judge documentation health from `docs-artifacts`, and do not
compare against merged PRs — those are green because their artifact landed long
ago.

## Code formatting

All C files, including headers, must be formatted with uncrustify using the
provided configuration `numcosmo_uncrustify.cfg`. Formatting is checked in CI.

Python code is checked with the configured `flake8`, `pylint`, and `mypy`
settings (see `.flake8`, `.pylintrc`, `.mypy.ini`).

## CI conda environment

The conda jobs in `.github/workflows/build_check.yml` do not solve
`environment.yml`. They install a committed lock file — an explicit list of
package URLs, already solved — so no solver runs on the runner:

```
.github/locks/conda-<subdir>-py<version>-<mpi>.lock
```

One file per CI target (`linux-64`/`osx-arm64`, python version, `openmpi`/`mpich`);
the target list lives in `TARGETS` in `.github/scripts/conda_locks.py` and must
match the `build-miniforge` matrices. Each lock is `environment.yml` with python
pinned, the MPI package swapped, and `libfabric-devel` added — the same
environment the jobs used to build by solving.

The `setup-miniforge` action picks the lock matching the runner and the job's
`python-version`/`mpi` inputs, and uses it only if its recorded
`environment.yml` sha256 still matches. Otherwise it warns and solves, so a
missing or stale lock is slow, never broken.

The environment is no longer cached between runs: installing a lock is a
download and link of a fixed package list, comparable to restoring a
multi-gigabyte cache entry, and dropping it leaves the repository's cache quota
to the test-duration caches.

### Changing dependencies

Editing `environment.yml` invalidates every lock, and the `check-conda-locks`
job fails until they are regenerated:

```bash
pip install conda-lock          # or conda install -c conda-forge conda-lock
python .github/scripts/conda_locks.py generate
python .github/scripts/conda_locks.py check
git add environment.yml .github/locks
```

`generate` solves every target from the local machine — conda-lock solves for
other platforms too, macOS included — one process per target, since the solver
itself (libsolv) is single-threaded. All three take well under a minute.
`--platform`, `--mpi` and `--jobs` restrict or serialize the run. Never
hand-edit a lock file: the sha256 comment on its second line records the
`environment.yml` it came from and is what `check` compares.

If a solve suddenly takes minutes rather than seconds, suspect a version
conflict the solver is backtracking around rather than a slow machine. The
`c-compiler`/`fortran-compiler` pins in `environment.yml` exist for exactly
that reason: `pyccl` pulls in `camb`, which needs the gcc 14 toolchain, and
letting the compilers float to 15 turned a 20-second solve into an
hours-long one. `mamba create --dry-run` on the same specs prints the
conflict.

### Keeping them current

`.github/workflows/conda_lock_update.yml` regenerates all locks weekly (and on
`workflow_dispatch`), smoke-installs the linux ones, and opens a pull request
with the diff. Because pull requests opened by `GITHUB_TOKEN` do not start
workflows, close and reopen that PR to run the build matrix against the new
environments before merging. Merging it is the only thing that moves the CI
environment forward — locks never drift on their own.

Local development environments are unaffected: `conda env create -f
environment.yml` still solves normally (see [docs/install.qmd](docs/install.qmd)).

## Submitting contributions

1. Fork the repository (skip if you are on the development team).
2. Branch from the latest `master`.
3. Make your changes, keeping each commit independently building so the series is
   bisectable.
4. Ensure the file is registered, tested, formatted, and (if it carries physics)
   documented per the sections above.
5. Open a pull request with a clear description of the change.
6. Address review feedback from the maintainers.

## Code of Conduct

By contributing to NumCosmo you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).
