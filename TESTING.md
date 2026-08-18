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

### Declare the marker; the directory name is not one

The two halves of an opt-in capability gate read **different things**, and a file can fall
between them:

- `conftest.py` skips on `"<cap>" in item.keywords`. `item.keywords` carries the name of
  every parent node — **including the directory** — so a file under a directory named after
  a capability is treated as gated whether or not it declares anything.
- The capability lane selects with `-m <cap>`, and `-m` evaluates **declared markers only**.
  A directory name is not a marker.

So a file placed in a capability-named directory *without* the marker is skipped in the
default lane (keyword matched, no `--run-*` flag) **and** deselected in the capability lane
(no marker): it runs nowhere. Nothing fails — it reports as an ordinary skip — so this is
silent unless someone looks.

Always declare it at module top, next to the imports:

```python
pytestmark = pytest.mark.powspec
```

To confirm a new file is actually selected by its lane:

```sh
pytest --collect-only -q -m powspec tests/python/<path>   # must list the tests
```

### Parallelism and the `OMP_NUM_THREADS` pin

NumCosmo has three distinct parallel mechanisms, which must not be confused:

- **OpenMP threads** — `#pragma omp parallel [for]`, governed by `OMP_NUM_THREADS`. Used by
  `ncm_fit_esmcmc`, `ncm_fit_mc`, `nc_data_cluster_wl` (when `enable_parallel`),
  `ncm_stats_dist`, `ncm_stats_dist_vkde`. Note: `NcmFitMC`/`NcmFitESMCMC`'s
  `set_use_threads(bool)` only gates whether the OpenMP-parallel code path runs at all; the
  *actual* thread count is always `OMP_NUM_THREADS`. There is **no** separate NumCosmo thread
  pool for these — it is OpenMP. (`NcmFitMCMC` has no threading support at all — it always
  runs single-threaded.)
- **MPI** — separate processes via `mpiexec`; some objects in `perturbations/` branch on both
  OMP and MPI.
- **OpenMP SIMD** — `#pragma omp simd` (e.g. `ncm_sbessel_ode_solver`) is *vectorization*, not
  threading; **unaffected** by `OMP_NUM_THREADS`.

The thread policy is set **per test by meson** (`tests/meson.build`), not via exported env
vars — so a plain `meson test` does the right thing on any machine, and the OpenMP branches
are exercised in the normal run rather than only on a hand-invoked lane:

| Test kind | Scheduling | `OMP_NUM_THREADS` |
|-----------|-----------|-------------------|
| concurrent (default `is_parallel: true`; pytest `-n auto`) | many at once | `1` |
| run-alone OpenMP (`omp` suite ⇒ `is_parallel: false`; non-xdist `py-omp` lane) | one at a time, owns the machine | CPUs available *at test-run time* (all of them) |

Concurrent tests are pinned to a single thread on **every** backend (OMP + OpenBLAS/BLIS/MKL)
to avoid cores² oversubscription and the mixed-OpenBLAS deadlock. Tests with an
OpenMP-parallel path carry the **`omp`** suite tag: meson marks the C ones `is_parallel: false`
and runs them — and the non-xdist `py-omp` pytest lane — through
`tests/scripts/detect_omp_threads.sh`, a wrapper that sets `OMP_NUM_THREADS`/`OMP_THREAD_LIMIT`
to the CPU count available right then (`nproc` on Linux, `sysctl -n hw.ncpu` on macOS, falling
back to `getconf _NPROCESSORS_ONLN`) before `exec`ing the real test/`pytest`, so the parallel
branch (thread coordination, `reduction`, scheduling) actually runs. This is deliberately
evaluated fresh on every `meson test` invocation (not baked in at `meson setup` time) so a
builddir built on one machine can run correctly on a differently-sized one. Because
these tests run alone, this happens inside the ordinary `meson test` invocation — no
`--num-processes=1`, no exported `OMP_NUM_THREADS`, no separate lane required. (`mpi` still
gets its own `mpiexec` lane, pinned to one thread; OpenMP SIMD is unaffected by
`OMP_NUM_THREADS`.)

### Reproducibility and the `-Dflaky_tests` option

The C tests randomize via `g_test_rand*` (glib seeds it per process). Two modes, selected by
the meson option `flaky_tests` (default **false**):

| Mode | `-Dflaky_tests` | g_test seed | Repeats | Run by |
|------|-----------------|-------------|---------|--------|
| **non-flaky** (default) | `false` | fixed `--seed` (passed to every C test) | 1× | every push — bit-reproducible, never flakes |
| **flaky** | `true` | fresh random seed each run | `meson test --repeat=N` | the weekly `weekly_flaky.yml` lane |

Empirically the suite is well-calibrated: every randomized test passes 20/20 at OMP=1, and the
omp tests pass 20–30/20–30 at `OMP_NUM_THREADS=cores`. The fixed seed makes the push lane
deterministic; the weekly randomized + repeated lane is what stresses the seeded/statistical
tolerances across many draws. A fixed seed does **not** determinize an OMP-threaded run (parallel
FP/RNG order varies), so any omp-path statistical assertion must tolerate that spread — e.g.
`test_ncm_fit_esmcmc`'s variance-of-m2lnL check uses reltol 0.6 (not 0.4) for this reason.

---

## 4. CI sharding is derived, not labeled

Do not add suites like `stats-dist`, `fit-esmcmc`, `data-cluster-wl` to balance CI. Wall
time is balanced by `meson test --slice ${i}/N` over a `slice: [1..N]` matrix; tests
auto-distribute as they are added. Keep `priority: 10` on the heaviest executables so they
start first. Note `--slice` balances by test **count**, so genuinely huge single
executables must be **split** into multiple executables first (one cost class each).

---

## 5. Running the suite locally

Plain `meson test -C <builddir>` is the normal case and needs no arguments: meson already
pins the concurrent tests to one thread and runs the `omp` ones alone with the machine
(see §3).

Pick `--num-processes` to match the **physical** core count, not the logical one. The
concurrent tests are each pinned to a single thread, so one process per physical core
saturates the machine; going past that (e.g. using the hyperthread count, which is what
`getconf _NPROCESSORS_ONLN` and `nproc --all` report) oversubscribes the real execution
units and loses throughput on this CPU-bound floating-point work. The same physical-core
count the `omp` wrapper computes is available with:

```sh
awk '/^physical id/ { p = $NF } /^core id/ { print p "," $NF }' /proc/cpuinfo | sort -u | wc -l
```

### Reproducing the CI coverage job

Configure a coverage builddir with the same options CI uses, then capture with `lcov`
driven by `.lcovrc` (which carries the `g_error`/`g_assert*`/`G_OBJECT_WARN` exclusions —
reading raw `gcov` instead ignores them and overstates the gap):

```sh
meson setup Coverage -Dbuildtype=debug -Db_coverage=true -Dnumcosmo_py=true \
      -Dfftw-planner=estimate -Dnumcosmo_debug=disabled
meson compile -C Coverage

lcov --config-file .lcovrc --no-external --capture --initial \
     --directory Coverage --directory numcosmo --directory tests \
     --base-directory "$PWD/Coverage" --output-file cov-base.info

# One lane at a time -- see the warning below.
for s in c c-acceptance c-statistical python py-acceptance py-omp py-powspec py-xcor py-app; do
    meson test -C Coverage --timeout-multiplier 0 --num-processes=<physical-cores> --suite "$s"
done

lcov --config-file .lcovrc --no-external --capture \
     --directory Coverage --directory numcosmo --directory tests \
     --base-directory "$PWD/Coverage" --output-file cov-tests.info
lcov --config-file .lcovrc --add-tracefile cov-base.info --add-tracefile cov-tests.info \
     --output-file cov-full.info
lcov --config-file .lcovrc --remove cov-full.info '*/external/*' '*/tools/*' \
     --output-file numcosmo-coverage.info
```

`.gcda` counters accumulate across runs inside one builddir, so running the lanes in
sequence and capturing once at the end is equivalent to CI's per-job tracefile merge.
Delete stale `.gcda` (`find Coverage -name '*.gcda' -delete`) before a fresh measurement;
editing a source file mid-run regenerates its `.gcno` and makes `lcov` fail that file with
`stamp mismatch with notes file`.

**Run the lanes one at a time.** CI gives each lane its own runner. Every `py-*` lane is a
*single* meson test that internally spawns `pytest -n auto`, so passing several `--suite`
flags to one `meson test` lets meson start multiple lanes concurrently, each with a full set
of xdist workers holding a NumCosmo instance — enough to exhaust memory on a workstation.

### Selecting a subset

Tests are selected by `--suite`, not by test name. Opt-in capability lanes additionally need
their `--run-*` flag when invoked through `pytest` directly:

```sh
meson test -C <builddir> --suite c              # C lane
meson test -C <builddir> --suite python         # default (fast) python lane
pytest --run-powspec -m powspec tests/python    # one capability lane, directly
```

---

## 6. Naming

- Files: `test_<source_object_or_feature>.{c,py}`, mirroring the source object
  (`ncm_spline.c` → `test_ncm_spline.c` / `test_ncm_spline.py`).
- Keep basenames unique across the whole Python tree (required by importlib mode).
- Shared fixtures/builders live in `conftest.py` (per-dir conftest for module-local
  fixtures; top-level `tests/python/conftest.py` for cosmo/mset/rng builders). Do not
  copy-paste setup between files.

## 7. Determinism & golden data

- Statistical tests pin an explicit **seeded** RNG (a shared fixture), not the implicit
  global GSL seed. Tolerances are pinned to that seed.
- Golden references use NumCosmo serialization into `data/truth_tables/`, loaded via
  `Ncm.cfg_get_data_filename` and compared with `np.testing.assert_allclose`. No byte-hash
  snapshots, no large hardcoded value arrays.

## 8. Adding a test (checklist)

1. Put it in the directory of the module it covers (create the dir if missing).
2. Default is **unit** — keep it fast and deterministic. If it is random, make it
   `statistical` with a seeded RNG; if it is full integration, make it `acceptance`.
3. If it needs an optional dependency, add the matching capability marker and `--run-*`
   gate; declare any new marker in `pytest.ini`.
4. If the file lives under a capability-named directory, declare `pytestmark` for that
   capability — the directory name alone gates it without selecting it, and the file then
   runs in no lane at all. Verify with `pytest --collect-only -m <cap>`.
5. Do **not** add a CI-shard suite. Reuse fixtures from `conftest.py`.
