# NumCosmo directory restructure — migration plan

Branch: `directory-restructure`. Target: a major-version boundary (public include
paths change). Tests follow sources as a hard invariant.

## 1. Principles

1. Directory tree mirrors the two GObject namespaces: `ncm/` (foundation) and `nc/`
   (cosmology). **GI namespaces (`NumCosmoMath`, `NumCosmo`) are independent of
   directory layout — the public introspection API does not change.**
2. Each abstract base lives *with* its concrete subclasses (no base-here/impl-there).
3. Vendored third-party code is quarantined under `external/`.
4. Every source move is paired with its test move (`tests/c/...`, `tests/python/...`).
5. Include-path rewrites are scripted, never hand-edited (~1,260 internal +
   ~1,040 public-style references).

## 2. Public include paths — DECIDED: clean break

Moving headers changes their installed path (`install_headers(preserve_path:true)`),
e.g. `<numcosmo/math/ncm_vector.h>` → `<numcosmo/ncm/algebra/ncm_vector.h>`.

**Decision: clean break at the new major version. No compat shims.**

Mitigation that makes this safe: the *documented* entry points are the umbrella
headers `numcosmo.h` and `numcosmo-math.h` — downstream code is instructed to include
only those, never individual subdir headers. So anyone following the documented usage
is unaffected; only code that reached past the umbrellas into private subdir paths
breaks, and that was never supported. Action: keep the two umbrella headers as the
stable public surface, regenerate their `#include` lists for the new paths, and note
the break in ChangeLog + a migration note.

## 3. Target tree

```
numcosmo/
  ncm/            # foundation (was math/)
    core/         # cfg, serialize, util, timer, rng, memory_pool, c (constants),
                  #   func_eval, function_cache, iset, pln1d, dtuple, obj_array,
                  #   cblas, flapack, lapack, gsl_blas_types
    algebra/      # vector, matrix, quaternion, nnls
    spline/       # spline*, spline2d*, ode_spline
    specfunc/     # sf_sbessel, sf_spherical_harmonics, mpsf_*, sbessel_integrator*,
                  #   sbessel_ode_solver, binsplit*
    integration/  # integral1d*, integral_nd, integrate, diff
    stats/        # stats_vec, stats_dist*, bootstrap, function_sample_set
    model/        # model*, mset*, sparam, vparam, reparam*  (the model/MSet system)
    fit/          # fit*, likelihood, lh_ratio*, prior*, mset_trans_kern*, mset_catalog
    data/         # data, dataset, data_gauss*, data_dist*, data_poisson, data_funnel,
                  #   data_rosenbrock, data_gaussmix2d   (generic data framework)
    powspec/      # powspec, powspec_corr3d, powspec_filter, powspec_sphere_proj,
                  #   powspec_spline2d
    fftlog/       # fftlog*
    sphere/       # sphere_map, sphere_nn, spectral
    mpi/          # mpi_job*
    dynamics/     # csq1d  (generic time-dependent oscillator / quantum-mode evolver;
                  #   the namespace-agnostic ENGINE, no cosmology in it)
  nc/
    background/   # hicosmo(+priors), distance, scalefactor,
                  #   + ALL model/nc_hicosmo_* subclasses kept together (one
                  #   polymorphic family): lcdm, de* (xcdm/cpl/jbp/wspline/+reparams),
                  #   gcg, idem2, kinematic q(z) models (qspline/qconst/qlinear/qrbf),
                  #   and quantum/bounce backgrounds (Vexp, qgrw, qgw)
    quantum/      # standalone quantum/early-universe objects that are NOT HICosmo
                  #   subclasses: hiqg_1d (minisuperspace QG), de_cont (contracting-
                  #   universe perturbations). See "tangential physics" note below.
    recomb/       # recomb, recomb_cbe, recomb_seager
    reion/        # hireion, hireion_camb(+reparam_tau)
    primordial/   # hiprim + model/nc_hiprim_* (atan, bpl, expc, power_law, sbpl, two_fluids)
    powspec/      # powspec_ml*, powspec_mnl* + transfer_func* + window* + growth_func
    cmb/          # planck_fi*, cbe(+precision)
    perturbations/# nc_hipert* (unchanged set)
    supernova/    # snia_dist_cov  (or leave at nc/ root — see notes)
    lss/
      halo/       # density_profile*, cm_* (concentration-mass), bias*, multiplicity*,
                  #   mass_function, mass_summary, position, catalog(+generator,
                  #   +member_generator)
      cluster/    # cluster_mass*, cluster_redshift*, cluster_photoz_gauss*,
                  #   cluster_abundance, cluster_pseudo_counts,
                  #   reduced_shear_cluster_mass, cor_cluster_cmb_lens_limber
      galaxy/     # (today's galaxy/ whole) galaxy_sd_*, galaxy_wl_obs, galaxy_hod*,
                  #   galaxy_acf, galaxy_selfunc
      wl/         # wl_surface_mass_density, reduced_shear_calib(+wtg)
    xcor/         # nc_xcor* (unchanged set)
    data/         # Nc likelihood data: data_bao*, data_cmb*, data_cluster*,
                  #   data_snia*, data_hubble*, data_dist_mu, data_planck_lkl, data_xcor
  external/       # vendored: class, plc, levmar, libcuba, lintegrate, sundials,
                  #   toeplitz, numerics(=misc: cubature, Faddeeva, kdtree, libqp,
                  #   LowRankQP, rcm, nnls, cqp/mg_pdip/initial_point)
```

### Tangential / "not exactly cosmology" physics — classification rule
The split is **engine vs application**, not "physics vs cosmology":
- Generic machinery with no cosmology baked in → `ncm/` (foundation). `ncm_csq1d`
  (a time-dependent oscillator evolver) qualifies — its current Ncm membership is
  correct; it gets `ncm/dynamics/`. Future generic mode-function/WKB/quantum-evolution
  solvers join it there (it may start as a near-solo bucket — that's fine, it names
  the category).
- Cosmology-specific model that is a subclass → lives with its base. The quantum/bounce
  backgrounds (Vexp, qgrw, qgw) are `NcHICosmo` subclasses, so they stay in
  `nc/background/` with the rest of the family — do NOT fragment a polymorphic family
  by sub-theme.
- Cosmology-specific but standalone niche object → a named `nc/` bucket: `nc/quantum/`
  for `hiqg_1d` (minisuperspace QG) and `de_cont`.
- CSQ1D-based perturbation modes (`nc_hipert_adiab/gw/em`, `itwo_fluids`) stay in
  `nc/perturbations/` — being CSQ1D systems is an implementation detail; they are
  perturbation theory.

### Other judgement calls (flag for review)
- `ncm/dynamics/` could instead fold into `ncm/integration/` if a one-file subdir feels
  thin. Low stakes.
- `nc/quantum/` naming — alternatives: `nc/early_universe/`, `nc/minisuperspace/`.
- `nc/supernova/` for a single `snia_dist_cov` model may be over-granular; alternative
  is leaving it at `nc/` root next to its data object in `nc/data/`.
- `transfer_func*`, `window*`, `growth_func` moved from `lss/` into `nc/powspec/` as
  linear-theory ingredients (they feed the matter power spectrum, used beyond LSS).
  If you'd rather keep them as the halo-model foundation, they go to `nc/lss/halo/`.
- `mset_catalog` (MCMC chain store) placed in `ncm/fit/`; it is distinct from the new
  `ncm/.../catalog` data primitive (which is `ncm_catalog`, in `ncm/data/` or `core/`).

## 4. Mechanical procedure (per move-group)

For each leaf group (one new directory), one commit:

1. `git mv` each source/header/.c/.h into the new dir.
2. Move paired tests: `tests/c/<old>/test_*.c` and `tests/python/<old>/test_*.py`
   to the mirrored new test dir.
3. Rewrite includes repo-wide with a scripted pass (both forms):
   - `"oldsubdir/foo.h"` → `"newsubdir/foo.h"`
   - `<numcosmo/oldsubdir/foo.h>` → `<numcosmo/newsubdir/foo.h>`
   Use a generated sed script keyed off the actual moved filenames (not blanket
   dir renames — files fan out to several destinations).
4. Update `numcosmo/meson.build`: the `ncm_sources`/`ncm_headers`/`nc_sources`/
   `nc_headers` `files()` lists (subdir-qualified paths), and the umbrella headers
   `numcosmo.h` / `numcosmo-math.h`.
5. Update `tests/c/meson.build` + `tests/meson.build` registrations.
6. Build (`ninja -C Optimized`) + run the affected C and Python suites. Commit.

Keep each commit independently building so the series is bisectable.

Consider splitting `numcosmo/meson.build` into per-subdir `meson.build` via `subdir()`
(as the vendored dirs already do) once files are in place — optional cleanup phase.

## 5. Phasing (low-risk order: leaves first, namespace roots last)

1. **external/** — move vendored dirs (`class, levmar, libcuba, lintegrate, plc,
   sundials, toeplitz, misc→numerics`). Pure relocation; updates the top `subdir()`
   calls + `include_directories` (sundials_inc, libcuba_inc, lintegrate_inc) +
   internal `#include "misc/..."`. No first-party code restructured yet.
2. **nc/ cosmology grouping** — background, recomb, reion, primordial, powspec, cmb,
   supernova; fold `model/` and the loose top-level `nc_*.c` into them. `perturbations/`
   and `xcor/` just gain the `nc/` prefix.
3. **nc/lss/** — split current `lss/` into `halo|cluster|wl`; bring `galaxy/` under
   `lss/` whole; move transfer/window/growth into `nc/powspec/`.
4. **nc/data/** — relocate `data/` under `nc/data/`.
5. **ncm/ subdivision** — the big one (142 files → ~13 buckets, 871 internal refs).
   Do bucket-by-bucket, each its own commit, so failures localize.
6. **Umbrella + cleanup** — finalize `numcosmo.h`/`numcosmo-math.h`, pkg-config,
   `.gir` source paths, docs include paths; optional per-subdir meson.build;
   optional compat-shim phase (section 2).

Each phase: build + full `meson test` green before moving on. Regenerate `.pyi`
stubs only if any public signature changed (layout moves alone don't).

## 6. Things to touch beyond sources (don't forget)
- `numcosmo/meson.build` source/header lists + umbrella headers.
- `tests/c/meson.build`, `tests/meson.build`.
- `docs/` — Doxygen/Quarto input paths if they enumerate subdirs.
- `.gir` generation source lists (driven by the header lists — should follow).
- `numcosmo_export.sh.in`, pkg-config `.pc` (include dir is unchanged: still
  `numcosmo/`; only sub-paths change).
- `numcosmo_py` C-level imports are via GI namespace, unaffected.
```
