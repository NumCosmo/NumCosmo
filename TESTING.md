# NumCosmo test organization policy

This document defines **where a test file lives** and **how it is labeled**. It is
normative: new tests must follow it, and existing tests are being migrated to it.

## Guiding principle: one concept, one mechanism

A test is classified along three *independent* axes. Each axis has exactly **one**
mechanism, so the axes never get conflated:

| Axis | Question it answers | Mechanism |
|------|---------------------|-----------|
| **Module** | *What part of NumCosmo is under test?* | **Directory** (mirrors the source tree) |
| **Tier** | *How fast/deterministic is it?* | **Marker** (Python) / **suite** (C) |
| **Capability** | *Does it need an optional dependency or runtime?* | **Marker + `--run-*` opt-in** |

CI shard balancing is **not** an axis: it is derived mechanically (`meson test --slice`),
never encoded as a hand-written label. Do not invent per-shard suites.

---

## 1. Directory layout (the Module axis)

Directories mirror the NumCosmo source tree. The same module names are used for **both**
C and Python tests, so a reader finds the same place in both trees.

> **The layout tracks the source tree.** The module list below reflects the *current*
> `numcosmo/` and `numcosmo_py/` organization. When the sources/packages are reorganized
> (a planned future PR), the tests move with them in that same PR — the mapping is expected
> to change, the *principle* (tests live where their source lives) does not. So the rule
> when adding a test today is simply: put it next to where its source lives *now*.

```
tests/
  c/
    <module>/test_<thing>.c
  python/
    <module>/test_<thing>.py          # tests of the C library (Ncm*/Nc* via GI)
    numcosmo_py/<subpkg>/test_*.py     # tests of the pure-Python numcosmo_py package
```

### C-library modules (mirror `numcosmo/<dir>/`)

| Dir | Covers | Source |
|-----|--------|--------|
| `math` | `Ncm*` core: vectors, splines, stats, fit, MCMC, FFT, special functions | `numcosmo/math` |
| `model` | model framework, `NcHICosmo`, recombination, distance | `numcosmo/model` |
| `data` | likelihood data objects (`NcmData*`, BAO, SNIa, Hubble, cluster counts) | `numcosmo/data` |
| `lss` | large-scale structure: halo mass function, bias, cluster abundance, transfer, power spectrum | `numcosmo/lss` |
| `galaxy` | galaxy sample distributions and weak lensing | `numcosmo/galaxy` |
| `perturbations` | `NcHIPert*` perturbation theory | `numcosmo/perturbations` |
| `xcor` | cross-correlations | `numcosmo/xcor` |
| `misc` | misc utilities | `numcosmo/misc` |

Vendored subtrees (`class`, `levmar`, `libcuba`, `lintegrate`, `plc`, `sundials`,
`toeplitz`) are **not** test modules.

### Pure-Python modules (mirror `numcosmo_py/<subpkg>/`, Python tree only)

`numcosmo_py/app`, `numcosmo_py/catalog`, `numcosmo_py/sky_match`,
`numcosmo_py/analysis`, `numcosmo_py/plotting`, … — one dir per `numcosmo_py` submodule.
These test Python code, not the GI bindings, and live under `tests/python/numcosmo_py/`.

### Python import mode

`tests/python/pytest.ini` sets **`--import-mode=importlib`**. This stops pytest from
inserting test directories onto `sys.path`, which is what lets a module dir be named
`math` (or `data`, etc.) without shadowing a stdlib/third-party module. Requirement:
**test file basenames are unique across the whole tree** (they are today; keep it so).
`__init__.py` files in test dirs are optional under importlib and may be removed.

---

## 2. Tiers (the Tier axis)

Every test belongs to exactly one tier. The default lane (plain `pytest`, plain
`meson test`) runs **unit only**; heavier tiers are opt-in or run in dedicated CI jobs.

| Tier | Meaning | Rules |
|------|---------|-------|
| **unit** | Pinpoints one function/behavior: *does the thing that does X actually do X.* | Fast (well under a second typically), deterministic. **Default** — no marker needed. Runs on every PR. |
| **statistical** | Exercises numerics whose result is random: sampling, resampling, estimators, convergence-free invariants. | Must use a **seeded** RNG and assert with tolerances, not exact values. Off the fast lane. |
| **acceptance** | Full integration / posterior recovery (e.g. ESMCMC convergence, end-to-end pipelines). | Slowest; minimal config on PRs, full config in the weekly job. Prefer replacing "posterior converged" with cheap invariants where possible. |

### How to mark

- **Python:** `@pytest.mark.statistical` / `@pytest.mark.acceptance`. Unmarked = unit.
- **C:** meson `suite` field: `['c']` (unit, default), `['c-statistical']`, `['c-acceptance']`.

### Fuzz exception (C only)

C tests using GLib `g_test_rand_*` pick a fresh seed each run; their randomness is a
feature (rare bugs surface over many CI runs). They are **not** made deterministic and are
**not** "statistical" in the sense above — they stay in the unit lane and additionally
carry a `fuzz` suite tag so the weekly job can `--repeat` them. See the repeat strategy in
the project test-overhaul plan.

---

## 3. Capabilities (the Capability axis)

For tests that need an optional dependency or a special runtime. There are **two** gating
mechanisms, by intent:

**(a) Opt-in capabilities** — runnable in most environments but skipped by default because
they are heavy or special. Gated by a `--run-*` flag (the skip logic lives in
`tests/python/conftest.py`), selected with `-m`:

| Marker | Needs | Enable with |
|--------|-------|-------------|
| `mpi` | MPI runtime (`mpiexec`) | `--run-mpi` |
| `app` | CLI dependencies / heavy app flows | `--run-app` |
| `powspec` | power-spectrum extras | `--run-powspec` |
| `xcor` | cross-correlation extras | `--run-xcor` |
| `sphere_map` | sphere-map extras | `--run-sphere-map` |
| `omp` | `OMP_NUM_THREADS>1` to exercise the OpenMP-parallel branch | dedicated lane (see below) |

**(b) Optional-dependency tests** — should run *whenever the dependency is installed* and
skip silently otherwise. Gated at module top by `pytest.importorskip("<dep>")`, **not** a
`--run-*` flag:

| Dependency | Marker | Gate |
|-----------|--------|------|
| `pyccl` | `ccl` (selection label) | `pytest.importorskip("pyccl")` |
| `getdist`, `astropy`, `healpy`, … | — | `pytest.importorskip(...)` |

Capability is orthogonal to tier: a test can be `unit` *and* `ccl`. **All markers are
declared once in `tests/python/pytest.ini`** (the single source); `conftest.py` only wires
the `--run-*` options and their skip logic — it must not re-declare markers.

### Parallelism and the `OMP_NUM_THREADS` pin

NumCosmo has three distinct parallel mechanisms, which must not be confused:

- **OpenMP threads** — `#pragma omp parallel [for]`, governed by `OMP_NUM_THREADS`. Used by
  `ncm_fit_esmcmc`, `ncm_fit_mc`, `nc_data_cluster_wl` (when `enable_parallel`),
  `ncm_stats_dist`, `ncm_stats_dist_vkde`. Note: ESMCMC's `set_nthreads(N)` only toggles the
  `if (nthreads > 1)` clause; the *actual* thread count is `OMP_NUM_THREADS`. There is **no**
  separate NumCosmo thread pool for these — it is OpenMP.
- **MPI** — separate processes via `mpiexec`; some objects in `perturbations/` branch on both
  OMP and MPI.
- **OpenMP SIMD** — `#pragma omp simd` (e.g. `ncm_sbessel_ode_solver`) is *vectorization*, not
  threading; **unaffected** by `OMP_NUM_THREADS`.

The two lanes are inverses of each other; in both, total live threads ≈ runner cores:

| Lane | Process parallelism | `OMP_NUM_THREADS` |
|------|--------------------|-------------------|
| **default** | many: `--num-processes=$(getconf _NPROCESSORS_ONLN)` / pytest `-n auto` | `1` |
| **omp** | **off**: `--num-processes=1` / pytest **without** xdist | `$(getconf _NPROCESSORS_ONLN)` (all cores) |

CI pins `OMP_NUM_THREADS=1` (and BLAS threads) on the **default lane** to avoid cores²
oversubscription and the mixed-OpenBLAS deadlock. Consequence: the OpenMP-parallel *branches*
(thread coordination, `reduction`, scheduling) are **not exercised** on the default lane. To
cover them, tests with an OpenMP-parallel path carry the **`omp`** marker and are
**additionally** run in the dedicated omp lane above: process-level parallelism **disabled**
so the one test process owns the machine, and `OMP_NUM_THREADS` set to the cores available on
the runner so the parallel path actually uses them. This mirrors how `mpi` gets its own
`mpiexec` lane. `omp` is a selection label (the same tests still run serially on the default
lane); it is not a skip.

---

## 4. CI sharding is derived, not labeled

Do not add suites like `stats-dist`, `fit-esmcmc`, `data-cluster-wl` to balance CI. Wall
time is balanced by `meson test --slice ${i}/N` over a `slice: [1..N]` matrix; tests
auto-distribute as they are added. Keep `priority: 10` on the heaviest executables so they
start first. Note `--slice` balances by test **count**, so genuinely huge single
executables must be **split** into multiple executables first (one cost class each).

---

## 5. Naming

- Files: `test_<source_object_or_feature>.{c,py}`, mirroring the source object
  (`ncm_spline.c` → `test_ncm_spline.c` / `test_ncm_spline.py`).
- Keep basenames unique across the whole Python tree (required by importlib mode).
- Shared fixtures/builders live in `conftest.py` (per-dir conftest for module-local
  fixtures; top-level `tests/python/conftest.py` for cosmo/mset/rng builders). Do not
  copy-paste setup between files.

## 6. Determinism & golden data

- Statistical tests pin an explicit **seeded** RNG (a shared fixture), not the implicit
  global GSL seed. Tolerances are pinned to that seed.
- Golden references use NumCosmo serialization into `data/truth_tables/`, loaded via
  `Ncm.cfg_get_data_filename` and compared with `np.testing.assert_allclose`. No byte-hash
  snapshots, no large hardcoded value arrays.

## 7. Adding a test (checklist)

1. Put it in the directory of the module it covers (create the dir if missing).
2. Default is **unit** — keep it fast and deterministic. If it is random, make it
   `statistical` with a seeded RNG; if it is full integration, make it `acceptance`.
3. If it needs an optional dependency, add the matching capability marker and `--run-*`
   gate; declare any new marker in `pytest.ini`.
4. Do **not** add a CI-shard suite. Reuse fixtures from `conftest.py`.
