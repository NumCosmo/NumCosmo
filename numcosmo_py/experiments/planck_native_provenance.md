# Native Planck 2018 likelihoods — provenance & attribution

NumCosmo provides independent reimplementations of the Planck likelihoods,
including support for alternative resampled representations. The underlying
likelihood methodology and original data products are due to the Planck
Collaboration. Analyses using these likelihoods should cite the corresponding
Planck publications in addition to NumCosmo.

The `numcosmo_py.experiments.planck_{commander,simall,smica,lite,lensing}`
converters and `generate_planck18_native` produce NumCosmo-native Planck 2018
likelihood objects. They are **value-added subsets** derived from the public
Planck 2018 (PR3) `plc_3.0` likelihood release, reprocessed into self-contained
serialized NumCosmo objects (the band-powers, binning, (inverse) covariances,
correction tables and Gaussianization splines are read out of the clik files and
stored directly, with no dependence on the clik/PLC library at use time).

## Source data

- **Origin**: Planck Legacy Archive, Planck 2018 baseline likelihood package
  `plc_3.0` (`baseline/plc_3.0/...`). NumCosmo does not modify the numerical
  content of the released files.
- **Blocks used by `generate_planck18_native`**:
  - low-ℓ EE — `low_l/simall/.../simall_100x143_offlike5_EE_Aplanck_B.clik` (SimAll)
  - low-ℓ TT — `low_l/commander/commander_dx12_v3_2_29.clik` (Commander gibbs_gauss)
  - high-ℓ TT / TTTEEE — `hi_l/plik/plik_rd12_HM_v22{,b}_{TT,TTTEEE}.clik` (SMICA plik)
  - lensing (optional) — `lensing/smicadx12_..._consext8.clik_lensing`

## Reproducibility (data-reduction path)

The converters *are* the reduction script: given a local `plc_3.0` tree, they
regenerate the native objects deterministically. To rebuild from the primary
archive:

1. Obtain `plc_3.0` from the Planck Legacy Archive (the same package clik uses).
   NumCosmo resolves it via `find_baseline_file` under the baseline data dir.
2. Run the converters (`build_commander`, `build_simall`, `build_smica_tt` /
   `build_smica_ttteee`, `build_lensing`, `build_plik_lite`) or
   `generate_planck18_native`.

## Downloadable release (value-added products)

The curated native likelihoods are distributed as self-contained serialized
objects (`planck_native_pr3_<block>.gvar` -- `pr3` is the Planck data release
they are reduced from, the 2018 `plc_3.0`/`R3.00` package, so a set built from
another release sits alongside them), hosted on the shared NumCosmo data-asset
release `datafile-release-v1.0.0` alongside the weak-lensing catalogs, the SNIa
covariances and the Planck clik baseline tarball
(`numcosmo_py.experiments.planck_native_release`):

- **Rebuild** the artifacts from a local `plc_3.0` tree with `build_release()` or
  `numcosmo generate planck18-release <out_dir>` (the data-reduction step).
- **Consume** them with `load_planck_release(id)` (downloads and caches from the
  `RELEASE_TAG` GitHub release), or `generate planck18 --native --from-release`.

Only these lightweight value-added `.gvar` products are hosted — never the raw
maps. Users without the release can always reproduce them from the PLA with the
converters above.

## Two clik bugs reproduced by design

The native code matches clik to machine precision, including two genuine clik
artifacts (documented in `numcosmo/external/plc/CHANGELOG.NumCosmo.md`):

- Commander's single-precision π in the `Dℓ` conversion — reproduced only under
  the opt-in `clik-pi-compat`; the native default uses correct double π.
- The lensing CMB-marginalized loader dropped the φφ N0 correction; NumCosmo
  applies the upstream fix (clik `44e638b`), so the native marginalized lensing
  differs from raw `plc_3.0`.

## Required attribution

The reimplementations here are NumCosmo's own; the likelihood methodology and
the original data products they are built from are the Planck Collaboration's.
Any use or redistribution must therefore cite the Planck Collaboration and the
relevant release papers *in addition to* NumCosmo:

- Planck Collaboration, *Planck 2018 results. V. CMB power spectra and
  likelihoods*, A&A 641, A5 (2020).
- Planck Collaboration, *Planck 2018 results. VIII. Gravitational lensing*,
  A&A 641, A8 (2020) (when the lensing block is used).
- Planck Collaboration, *Planck 2018 results. I. Overview*, A&A 641, A1 (2020).

Planck data products are public domain (ESA open-access policy); attribution and
documented provenance are the conditions for reuse and redistribution.
