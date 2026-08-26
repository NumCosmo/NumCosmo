# Truth tables

Reference values the test suite compares against, grouped by the subsystem that
owns them. They exist so a test can assert an *absolute* answer instead of only
checking a result against itself — a self-consistency test passes just as
happily when every path shares the same error.

Reached from a test through `Ncm.cfg_get_data_filename("truth_tables/<dir>/<file>", True)`,
never by a relative path: the file has to be found in the build tree and in an
installed prefix alike. The whole directory is installed by
`install_subdir('data', ...)` in the top-level `meson.build`, which recurses, so
adding a subdirectory here needs no build-system change.

## Layout

| directory | contents | consumed by |
|---|---|---|
| `sbessel/` | `int f(x) j_l(k x) dx` for a Gaussian and a rational `f`, over a grid in `l` and `k` | `tests/python/ncm/specfunc/test_sbessel_integrator_{levin,fftl}.py` |
| `halo/` | Castro multiplicity function and halo bias, against the CCToolkit implementation | `tests/python/nc/lss/halo/test_{bias_castro,multiplicity_func_castro_cctoolkit}.py` |
| `cluster/` | a seeded cluster-count resample, as a golden file | `tests/python/nc/data/test_ncount_resample_golden.py`, `tests/python/nc/lss/halo/test_halo_catalog_generator.py` |
| `sphere/` | a seeded rectangular sky footprint | `tests/python/ncm/sphere/test_sky_footprint.py` |
| `wl/` | generation/integration parity for the galaxy `Pop`/`Obs`/`Factor` hierarchy | `tests/python/nc/lss/galaxy/test_galaxy_*.py` |

## Formats

**`.bin`** — `NcmSerialize` binary dumps of the objects a test rebuilds and
compares. Written and read through the library, so they carry no separate
format description.

**`.json.gz`** — gzipped JSON, self-describing:

```
center, std, lower-bound, upper-bound   the shape and the integration range
kvals[i]                                 the k grid
lvals[j]                                 the multipole grid
table[j][i]                              the value at (lvals[j], kvals[i])
convention                               what table[j][i] means, exactly
```

The `convention` key is not decoration. These tables held `k` times the
documented integral for years, and nothing noticed: every test that used them
compared a result against another result at a different tolerance, where the
factor cancels exactly. It surfaced only when the values were first checked
against an outside reference. **State the convention, including the measure, in
any table added here.**

## Provenance, and how to regenerate

The `sbessel/` tables were produced in Mathematica; that notebook is not in the
repository, which makes them the least reproducible thing in this directory.
They have since been verified against certified Arb ball arithmetic (agreeing to
8 digits across `l` in [0, 500] and `k` in [0.1, 5623], down to values of
1e-17), and rescaled by `1/k` to match the documented `dx` measure.

Regenerating them from `tests/tools/nc_xcor_kernel_analytic_arb.c` would replace
that provenance with a program in-tree and attach a proved error radius to every
entry. The cost is not prohibitive — 501 x 101 entries at roughly 0.1 s each is
about 1.4 h on one core, a few minutes across twelve — and it is worth doing the
next time these tables are touched.

`gauss_jl_100.json.gz` currently has **no consumer**. It is kept only because
deleting data is not free to undo; delete it once you are sure.

## Adding a table

- Put it in the subdirectory of the subsystem that owns it, or add one.
- Prefer a self-describing format, and state the convention in the file.
- Say where the values came from, here, in enough detail to reproduce them.
- Reference it through `Ncm.cfg_get_data_filename`, not a relative path.
