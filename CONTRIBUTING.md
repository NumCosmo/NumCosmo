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
   `numcosmo/ncm/core/ncm_cfg.h` and register the object in `ncm_cfg.c` with
   `ncm_cfg_register_obj (NCM_TYPE_MY_OBJECT);`.

## Unit testing

Tests follow their sources as a hard invariant: a class in
`numcosmo/<ns>/<area>/` is tested in the mirrored `tests/c/<ns>/<area>/` and/or
`tests/python/<ns>/<area>/`. All new code must have at least minimal testing.

1. Add the test under the mirrored path (`tests/c/...` for C, prefixed
   `test_`; `tests/python/...` for Python, prefixed `test_py_`).
2. Register it in the corresponding `tests/c/meson.build` or
   `tests/python/meson.build`.
3. Run the affected suite with `meson test -v -C build`.

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

### Building the documentation

Build the site with a documentation-enabled build directory:

```bash
meson compile -C <builddir> numcosmo-site
```

Note: when you change **only** a `.c`/`.h` doc comment (no `.qmd`/`_quarto`
change), the site bundle may not re-render because the `numcosmo-site` target
does not detect changes in the generated API reference. Force it with:

```bash
ninja -C <builddir> -t clean docs/numcosmo-site
meson compile -C <builddir> numcosmo-site
```

This affects only local incremental builds; CI and Read the Docs build from
clean, so the published site is unaffected.

## Code formatting

All C files, including headers, must be formatted with uncrustify using the
provided configuration `numcosmo_uncrustify.cfg`. Formatting is checked in CI.

Python code is checked with the configured `flake8`, `pylint`, and `mypy`
settings (see `.flake8`, `.pylintrc`, `.mypy.ini`).

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
