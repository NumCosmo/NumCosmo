# Contributing to NumCosmo

Thank you for your interest in contributing to NumCosmo. This guide covers the
mechanics of adding code, tests, and documentation so that a new contribution
integrates cleanly with the build, the test suite, and the generated API
reference. [REVIEWING.md](REVIEWING.md) describes how contributions are reviewed, and
[CI.md](CI.md) how the CI jobs are set up.

For build and dependency setup, see the
[installation guide](https://numcosmo.readthedocs.io/en/latest/install.html).
NumCosmo uses the [meson](https://mesonbuild.com/) build system.

## Source layout

The library is split into two GObject namespaces, mirrored by the directory tree:

- `numcosmo/ncm/`, NumCosmoMath: foundation code with no cosmology
  (algebra, splines, integration, statistics, the model/MSet system, fitting).
- `numcosmo/nc/`, NumCosmo: cosmology models and likelihoods
  (background, perturbations, large-scale structure, clusters, CMB, …).

Each subdirectory groups a coherent family; an abstract base lives together with
its concrete subclasses. Vendored third-party code is quarantined under
`numcosmo/external/`.

## Adding a new class

1. Place the files in the subdirectory that matches the family, under
   `numcosmo/ncm/<area>/` or `numcosmo/nc/<area>/`.

2. Register the sources in `numcosmo/meson.build`. Add the `.c` to the
   `ncm_sources` (or `nc_sources`) `files()` list and the `.h` to `ncm_headers`
   (or `nc_headers`), using the subdirectory-qualified path, e.g.
   `'ncm/algebra/ncm_my_object.c'`.

3. Export through the umbrella header. Add the header to `numcosmo/numcosmo.h`
   (cosmology) or `numcosmo/numcosmo-math.h` (math). These two umbrellas are the
   only supported public include surface; never document reaching into a
   subdirectory header directly.

4. Register the type for serialization. Include the header in
   `numcosmo/ncm/core/ncm_cfg.c`, next to the other headers of its family, and
   register the object in `ncm_cfg_register_objects ()` with
   `ncm_cfg_register_obj (NCM_TYPE_MY_OBJECT);`. (`ncm_cfg.h` includes only what its
   own public API needs; do not add registration includes there.)

## Unit testing

Tests follow their sources as a hard invariant: a class in
`numcosmo/<ns>/<area>/` is tested in the mirrored `tests/c/<ns>/<area>/` and/or
`tests/python/<ns>/<area>/`. All new code must have at least minimal testing.

1. Add the test under the mirrored path, prefixed `test_`: `tests/c/<ns>/<area>/`
   for C, `tests/python/<ns>/<area>/` for Python.
2. C tests must be registered in `tests/c/meson.build`: add an entry giving the
   test `name` and its `sources`, since each C test is built as its own executable.
   Coverage of a class's ref/free/clear methods goes in the existing
   `tests/c/ncm/core/test_ncm_generic.c` rather than in a new file.
3. A C test file follows one layout, top to bottom: the license header and
   includes; `#define`s, types and static data tables; the prototypes of every
   fixture and test function; `main`, which only registers the tests; the helper
   functions; the fixture creators and destructors (the setup/teardown passed to
   `g_test_add`); then the test implementations. `tests/c/ncm/algebra/test_ncm_vector.c`
   is the model.
4. Python tests are collected automatically: the suites run `pytest` over
   `tests/python/`, so a new file needs no `meson.build` entry. Selection between
   lanes is by *marker*, not by path: a file with no marker joins the default lane,
   and one belonging to a dedicated shard (`xcor`, `powspec`, `app`, `omp`,
   `sphere_map`, `acceptance`) must declare a module-level
   `pytestmark = pytest.mark.<name>`.
5. Run the affected suite with `meson test -v -C build`, or a single file directly
   with `python -m pytest tests/python/<path> -q`.

A test is usually easier to develop in Python: the bindings, numpy and an
interactive session make it quick to probe a case and settle its tolerance. Once
the test is settled, move it to C when it exercises only the C API: C tests run
faster, are measured by the C coverage report, and do not depend on the bindings.
A Python test stays in Python when its subject is Python itself (the bindings,
`numcosmo_py`) or another code it compares against (CCL, CLASS). Remove the
Python test only after its C counterpart has every assertion it had.

## Documentation

Documentation is split by purpose, and the split is intentional:

- API reference (gi-docgen, generated from the GObject doc comments in the
  `.c` files): brief and operational. State what the class does, the defining
  relation or signature, and the key methods. Keep it short.
- Theoretical background (the Quarto project under `docs/`, rendered to the
  website): the physics, derivations, and full equations live here, in
  `docs/theory/<ns>/<area>/<topic>.qmd`, mirroring the `numcosmo/<ns>/<area>/`
  of the class it documents.

When a class involves non-trivial math, do not put the derivation in the C doc
comment. Put it on a theory page and link to it. The pattern is established by
`NcmCSQ1D`: see the short doc comment in
`numcosmo/ncm/dynamics/ncm_csq1d.c` and the corresponding
`docs/theory/ncm/dynamics/csq1d.qmd`.

- From the C doc comment, link to the theory page with a plain anchor, e.g.
  `<a href="../../theory/ncm/dynamics/csq1d.html">CSQ1D Formalism</a>`.
- From the theory page, link back to the API with wiki-style symbol references:
  `[[numcosmo-math|NcmCSQ1D]]`, `[[numcosmo|NcDistance]]`,
  `[[numcosmo-math|ncm_csq1d_prepare]]`. Unresolved references degrade to plain
  text, so verify the symbol name.

A single inline `$x$` in a doc comment is fine; a multi-line `\begin{align}`
derivation belongs on a theory page.

Only theory pages cite references: they cite with pandoc `[@Key]` against
`docs/references.bib`. API doc comments carry no references (gi-docgen has no
bibliography support); a doc comment that needs one links to the theory page that
cites it, and a theory page is written or extended when none covers it.

### Writing style for comments and documentation

Comments and docs are read by people holding a large amount of the library in
their head at once. Every phrase that has to be decoded before it can be used
costs them. Write so a sentence can be read once, at speed.

State facts, not impressions. Do not use evaluative or promotional wording
for the project's own work, and do not label something as self-evident. These
are rejected by the style check:

> shipped, blow up / blows up, blowup, comfortably, blindly, genuinely,
> decisive, safe zone, sweet spot, nets out, load-bearing, at length,
> worth knowing, needlessly, wins big, obvious / obviously, hot loop,
> of course, magic, trap / traps, story, honest / honestly,
> clever / cleverly, elegant, it turns out

Name the object instead of gesturing at it: not "each piece" but "the
`chi_I` series, the Jacobian series and their product".

State the current behavior, not how it got there. A comment says what the
code does and what it guarantees. The sequence of attempts that led to it, and
which of them failed, goes in the pull request description or an area history
document. These are rejected by the style check:

> historically, in the past, originally, nowadays, these days, at the time,
> the old code / path / version / implementation / way,
> used to abort / crash / hang / fail / refuse / stall / deadlock / leak /
> be broken / be wrong

Write the invariant instead: not "the pair that used to abort", but "a pair
with a panel-backed closure on either side". A past *interface* is sometimes
the point: a migration note has to say which property an old serialized file
set and what it now maps to, and "used to" and "no longer" on their own are
not flagged, so those notes are unaffected. A measured number keeps its
measurement: "measured at 975161 knots with the second differences never
coming down" is a fact about the code's behavior, not a history lesson.

Placement. A function or class doc states what the thing does and its error
modes, and stops. Derivations go on a theory page (see above). Measurement
tables, rejected approaches, and how a bug was found go in an area history
document such as `docs/theory/nc/lss/galaxy/wl_shape_factor_history.md`, with a
pointer left in the code.

Do not delete these when shortening a comment. Move them if they belong
elsewhere, but they must survive somewhere findable:

- error modes: what returns NaN, what aborts, what is clamped and why
- concurrency contracts: safe or unsafe to call in parallel, and what makes it so
- usage limits: the regime where the code gives a wrong answer, and what to use there
- cross-file invariants: "must be kept equal to `X` in `y.py`"
- measured accuracy envelopes, and the guidance drawn from them
- links: `<a href="../../theory/....html">` from C, `[[numcosmo|Sym]]` from `.qmd`

Checking. `.github/scripts/check_doc_style.sh` flags both word lists above in
added lines and reports theory-page links removed by the current branch. CI runs
it on pull requests. It is advisory for the link report and blocking for the two
word lists; to run it locally:

```bash
.github/scripts/check_doc_style.sh              # check the diff against origin/master
.github/scripts/check_doc_style.sh --all        # check the whole tree
```

### Bibliography (`docs/references.bib`)

`bibtex-tidy` formats `docs/references.bib`; do not format it by hand or with a
reference manager. After adding or editing an entry, run:

```bash
docs/tidy_references.sh          # re-format in place
docs/tidy_references.sh --check  # verify (this is what CI runs)
```

CI rejects a bibliography that is not tidy (`.github/workflows/bib_lint.yml`), so
run the formatter before committing. Install the pinned version with
`npm install -g bibtex-tidy@1.14.0`, since different versions format differently. A
`pre-commit` hook runs it automatically after
`pip install pre-commit && pre-commit install`. Output from a reference manager such
as JabRef is a draft to re-tidy; the local `file`, `owner`, `timestamp`, and
`__markedentry` fields are stripped from the repository copy.

### Adding a page

A new `.qmd` has to be registered in four places. Miss one and it either
fails to build, builds but is unreachable, or is reachable only from the
sidebar:

1. `docs/<area>/meson.build`: add the filename to the `files(...)` list, or
   it is never copied into the build tree. A new subdirectory also needs
   its own `meson.build` and a `subdir('<name>')` line in the parent.
2. `docs/_quarto.yml.in`, the render list: the top-level list of pages to
   render (near `- tutorials/index.qmd`).
3. `docs/_quarto.yml.in`, the sidebar: the `contents:` entry under the
   right `section:`, with an `href:` and a short `text:` label.
4. The index page for its area: `docs/tutorials/index.qmd`,
   `docs/examples/index.qmd`, or `docs/theory/index.qmd`. The sidebar and the
   index are separate; adding a page to the sidebar alone leaves it off the
   landing page a reader browses.

Keep the section names in 2-4 consistent: if a page needs a new section, add it
to both the sidebar and the index.

Tutorial or example? A tutorial walks through a complete workflow and
explains every step, and is read start to finish; it goes in `docs/tutorials/`.
An example is a short self-contained script showing one part of the API in
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

When you change only a `.c`/`.h` doc comment (no `.qmd`/`_quarto` change), the
site bundle may not re-render, because the `numcosmo-site` target does not detect
changes in the generated API reference. Force it with:

```bash
ninja -C <builddir> -t clean docs/numcosmo-site
meson compile -C <builddir> numcosmo-site
```

This affects only local incremental builds; CI and Read the Docs build from
clean, so the published site is unaffected. The Read the Docs check on a pull
request fails until the `docs-artifacts` job finishes; see [CI.md](CI.md).

## Code formatting

All C files, including headers, must be formatted with uncrustify using the
provided configuration `numcosmo_uncrustify.cfg`. Formatting is checked in CI.

Python code is formatted with `black` and checked with `ruff` and `mypy`, both
configured in `pyproject.toml`. CI also regenerates the stubs `numcosmo_py/nc.pyi`
and `numcosmo_py/ncm.pyi` and fails if they differ from the committed ones; after
changing the C API, build and run `numcosmo_py/update_pyi.sh`, then commit the stubs.

## Adding a dependency

A third-party Python module is declared in up to three lists, and
`.github/scripts/check_deps.py` compares each one with the import statements
in the tracked sources (`notebooks/` excluded):

| list | must contain | checked by |
|---|---|---|
| `pyproject.toml` | core: what `import numcosmo_py` loads; `app`, `catalog`, `plot`: what `numcosmo_py.app` (the `numcosmo` CLI), `numcosmo_py.catalog` and `numcosmo_py.plotting` load | `check_deps.py repo`, `lint` job |
| `environment.yml` | every module imported by the package, tests, docs, examples and tools, and every `pyproject.toml` entry except the `firecrown` extra | `check_deps.py repo`, `lint` job |
| conda-forge recipe `run:` | the core, the `app`, `catalog` and `plot` extras, and the imports of `examples/` | `check_deps.py feedstock`, weekly `feedstock_deps.yml` |

"Load" means an import at module level outside a `try` block, followed through
the package's own modules. An import inside a function, inside a `try` block or
through `pytest.importorskip` needs only the `environment.yml` entry. When a
module's import name differs from its distribution name (`gi` and
`pygobject`), or its PyPI name from its conda-forge name (`matplotlib` and
`matplotlib-base`), add it to `IMPORT_TO_DIST` or `PIP_TO_CONDA` in the script.
When the recipe uses another package that provides a distribution (`astropy-base`
for `astropy`), add it to `RECIPE_ALTERNATIVES`.

Editing `environment.yml` invalidates the CI lock files, and the `lint` job fails
until they are regenerated:

```bash
pip install conda-lock          # or conda install -c conda-forge conda-lock
python .github/scripts/conda_locks.py generate
python .github/scripts/conda_locks.py check
git add environment.yml .github/locks
```

[CI.md](CI.md) describes the lock files. Run
`python .github/scripts/check_deps.py feedstock` before a release: the recipe is
updated in its own repository, so the `lint` job cannot see it.

## Submitting contributions

1. Fork the repository (skip if you are on the development team) and branch from the
   latest `master`.
2. Keep one concern per commit, and keep each commit building so the series is
   bisectable. A change to the library, its Python helpers and its documentation is
   three commits, in that order: the C change with its C test and regenerated stubs,
   the Python code with its tests, then the pages.
3. Ensure each file is registered, tested, formatted, and (if it carries physics)
   documented per the sections above, and run the checks: uncrustify, `black`,
   `ruff`, `mypy`, `check_doc_style.sh`, `tidy_references.sh --check` when the
   bibliography changed, and the affected tests.
4. Write commit messages as described below, and the pull request description in the
   same form.
5. Open the pull request. Open separate pull requests only for pieces that are useful
   and reviewable alone; otherwise separate commits in one pull request are enough.
6. Address review feedback from the maintainers; [REVIEWING.md](REVIEWING.md) lists
   what a review checks.

### Commit messages

A commit message is a one-line title, a blank line, and a flat `*` list:

```
Add nc_halo_mass_function_f_lnR

* Add nc_halo_mass_function_f_lnR, the f(sigma) that dn_dlnR uses
* Results change: Bocquet 200c/500c dn_dlnR differs at rounding level
* Add test_nc_halo_mass_function.c
* Regenerate nc.pyi
```

- The title is imperative, in sentence case, without a `type:` prefix or a final period.
- Each bullet is one change, on one line, under about 100 characters. Name the affected
  symbol, file or option. A bullet that needs a clause of measurements or bounds
  becomes two bullets.
- Do not wrap bullets by hand: `generate-changelog.sh` builds `ChangeLog.md` from the
  messages and keeps manual line breaks.
- A change of behavior is the first bullet, "Behavior change: ...", with its size and
  how to restore the previous behavior. State measured effects and limitations as
  "Results change: ..." and "Known limitation: ...".
- No prose paragraphs, headings or tables, and no development history.
- Authors and co-authors are the people who contributed to the change, control its
  content, and answer for it. Credit a co-author with a `Co-authored-by: Name <email>`
  trailer. Tools, LLMs included, are never authors or co-authors.

### Numbers and claims

Every number in a commit message, comment, docstring or page must be reproducible:
give the test, the page cell or the command that computes it. A tolerance states the
measured value it was set from. A stated cause ("the difference is large because ...")
is checked by a measurement that varies that cause.

### Responsibility for the content

You answer for every line and every claim you submit, whatever tools you used to produce
them, as with a compiler, an editor or a spell checker. Read and run what you submit,
and check each number and stated cause yourself. Tools used are not disclosed.

## Code of Conduct

By contributing to NumCosmo you agree to abide by the
[Code of Conduct](CODE_OF_CONDUCT.md).
