# NumCosmo changes to the embedded Planck Likelihood Code (PLC/clik)

This directory is a **vendored, trimmed copy** of the Planck Likelihood Code
(`clik`). It is not a verbatim upstream checkout: NumCosmo carries a small set of
local additions, and selectively back-ports individual upstream fixes as the
native port matures.

This file tracks those divergences so the tree stays auditable.

## Provenance

- **Baseline**: `plc_3.0` — the Planck 2018 public likelihood release. This is
  the version used to reproduce the Planck 2018 published results.
- **Upstream reference for back-ports**: `clik` `16.0b1`
  (`clik_16.0b1-12-gcde855c6debd`, HEAD `9c6a194`, 2025-01-22), kept locally at
  `~/Projects/clik`. Upstream is far ahead of `plc_3.0` (new likelihoods:
  SPT3G, mspec, camspec revisions, …); NumCosmo does **not** sync wholesale.

## Policy

Only the subset of PLC actually exercised by NumCosmo's native port is kept in
sync, and only when a change is a genuine correctness fix. Upstream additions
that are new likelihoods or infrastructure NumCosmo does not use are ignored.
Bugs that affect reproduction of the Planck 2018 published results are kept
(and reproduced on the native side behind an opt-in flag) rather than silently
corrected — see below.

---

## Local NumCosmo additions (not in upstream)

These are NumCosmo-only helpers, present since before this changelog and required
by `nc_data_planck_lkl.c` to validate a likelihood against its stored
`check`/`check_param` anchor at load time:

- `clik_get_check_param()` in `clik.c` (added 2015-10-28).
- `clik_lensing_get_check_param()` in `clik_lensing.c` + declaration in
  `clik_lensing.h`.

Both assemble the fiducial `cl + nuisance` parameter vector and return the
stored check value; they must be preserved across any upstream back-port.

---

## Back-ported upstream fixes

### 2026-07-23 — lensing: read the first-order correction for all files

Back-ported upstream clik commit `44e638b` (K. Benabed, 2021-05-05,
*"Correct a nasty bug in lensing.c"*) into `clik_lensing.c` (`_clik_lensing_init`,
itype=4 loader).

- **Bug**: the loader gated reading of the `cors` first-order (N0)
  renormalization on `hascl != 0`. For the **CMB-marginalized** lensing file
  (`hascl` all zero) this dropped the φφ correction entirely, leaving the
  marginalized likelihood inconsistent with the full one.
- **Fix**: always size `nlt` to include the φφ block and read `cors` whenever
  the file ships a `cors` node (`cldf_haskey`).
- **Effect**: for the fiducial cosmology the CMB-marginalized `m2lnL` changes
  from 15.16 (buggy) to 8.95, now consistent with the full file's 9.44. The full
  file is unaffected.
- The local `clik_lensing_get_check_param()` helper was kept; only the loader
  hunk was taken.

Native counterpart: `NcDataPlanckLensing` / `build_lensing()` read `cors` for
both files and match the fixed clik to machine precision.

---

## Upstream changes deliberately NOT applied

- **commander (`gibbs/comm_gauss_br_mod_v3.f90`)** — the low-ℓ gibbs code
  computes `Dl = Cl·ℓ(ℓ+1)/2/π` with a **single-precision π**, a ~2.8e-8 scale
  on every `Dl` (~2.1e-6 in `m2lnL`). This is *unchanged* upstream (0 diffs in
  `16.0b1`) and is part of the Planck 2018 published result, so it is kept
  as-is. The native `NcDataPlanckCommander` uses correct double π by default and
  reproduces this artifact only under the opt-in `clik-pi-compat` property.

- **`plik/smica.c`** — upstream only swaps an OpenMP `schedule(dynamic)` for
  `schedule(auto)` (gcc9/Intel-2019 compatibility). No effect on results; not
  taken.

- **`lklbs.c`** — upstream changes are `sprintf` idiom tweaks and a new
  `options_table` API (used by SPT3G, which NumCosmo does not build). No change
  to the Cl → likelihood math; not taken.

- **`cmbonly/plik_cmbonly.f90` (plik_lite)** — unchanged upstream (0 diffs).
