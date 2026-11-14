# Continuous integration

This guide is for maintainers who change the CI workflows. It covers the conda
environment the jobs install, the cache of downloaded data files, and the Read the Docs
check on pull requests. A contributor who adds a dependency needs only "Adding a
dependency" in [CONTRIBUTING.md](CONTRIBUTING.md).

## Conda environment

The conda jobs in `.github/workflows/build_check.yml` install a committed lock file
instead of solving `environment.yml`. A lock file is an explicit list of package URLs,
already solved, so no solver runs on the runner:

```
.github/locks/conda-<subdir>-py<version>-<mpi>.lock
```

There is one file per CI target: `linux-64` or `osx-arm64`, a Python version, and
`openmpi` or `mpich`. The target list is `TARGETS` in `.github/scripts/conda_locks.py`
and must match the `build-miniforge` matrices. Each lock is `environment.yml` with Python
pinned, the MPI package swapped, and `libfabric-devel` added.

The `setup-miniforge` action picks the lock that matches the runner and the job's
`python-version` and `mpi` inputs. It uses the lock only if the `environment.yml` sha256
recorded in the lock matches the current file; otherwise it warns and solves. A missing
or stale lock therefore makes the job slower, not failed.

The environment is not cached between runs. Installing a lock downloads and links a
fixed package list, which takes about as long as restoring a multi-gigabyte cache entry,
and leaves the repository's cache quota to the data-file cache.

### Regenerating the locks

Editing `environment.yml` invalidates every lock, and the `lint` job fails until they
are regenerated:

```bash
pip install conda-lock          # or conda install -c conda-forge conda-lock
python .github/scripts/conda_locks.py generate
python .github/scripts/conda_locks.py check
git add environment.yml .github/locks
```

`generate` solves every target on the local machine, macOS included, with one process
per target, since the solver (libsolv) is single-threaded. The three commands take less
than a minute. `--platform`, `--mpi` and `--jobs` restrict or serialize the run. Never
edit a lock file by hand: the sha256 comment on its second line records the
`environment.yml` it came from, and `check` compares it.

A solve that takes minutes instead of seconds is usually a version conflict the solver
is backtracking around. `mamba create --dry-run` on the same specifications prints the
conflict. The `c-compiler` and `fortran-compiler` pins in `environment.yml` prevent one
such conflict: `pyccl` depends on `camb`, which needs the gcc 14 toolchain, and with the
compilers free to move to gcc 15 a solve took hours instead of 20 seconds.

### Weekly lock update

`.github/workflows/conda_lock_update.yml` regenerates all locks weekly, and on
`workflow_dispatch`, installs the Linux ones as a smoke test, and opens a pull request
with the diff. Pull requests opened by `GITHUB_TOKEN` do not start workflows, so close
and reopen that PR to run the build matrix on the new environments before merging it.
The CI environment changes only when such a PR is merged.

Local development environments do not use the locks: `conda env create -f
environment.yml` solves as usual (see [docs/install.qmd](docs/install.qmd)).

## Data-file cache

Several data sets are not in the repository: the Planck baseline (`plc_3.0`), the SNIa
covariance catalogs, the curated weak-lensing catalogs and the native Planck likelihood
objects. NumCosmo downloads each one from a tagged GitHub release the first time it is
requested and keeps it in the NumCosmo user data directory (see `ncm_cfg_init`). The
files total about 470 MB; without a cache every CI job downloads its share on every run.

`.github/actions/data-cache` sets `NUMCOSMO_HOME=$HOME/.local/share/numcosmo` for the
job and restores only the downloaded files into it, never the FFTW wisdom and FFTLog
tables the library writes next to them. The cache key has two parts, and only they make
an entry stale:

- The release tag. The key contains the `datafile-release-vX.Y.Z` tag read from the
  sources that name it, so changing the tag there is enough.
- `DATA_CACHE_VERSION` in `.github/workflows/build_check.yml`. Increase it for any change
  the tag does not capture: a file replaced under an unchanged tag, which the library
  would never download again, and a new asset under an unchanged tag, which every job
  would download on every run, since on a cache hit the writer job neither prefetches nor
  saves. `.github/workflows/weekly_flaky.yml` carries the same value.

There are no `restore-keys`. A partial match would restore a file whose name is
unchanged but whose contents the new release replaced, and the library downloads only
missing files, so the outdated copy would stay.

One job per OS writes the entry (`build-gcc-ubuntu` and `build-gcc-macos`); every other
job restores it. The writer first runs `.github/scripts/prefetch_data.py`, which
downloads every asset, not only those its own tests read, because the rest of the matrix
receives exactly what the writer stored until the key changes. The script reads the
asset list from the enumerations that define the assets, so a catalog added to
`NcDataSNIAId`, `NcGalaxyWLObsCatalogId` or `PlanckReleaseId` needs no change to it.

A cache miss is not a failure: every job downloads what it needs, and the next run finds
the entry.

## Read the Docs check on pull requests

Read the Docs does not build the documentation. The `docs-artifacts` GitHub Actions job
builds the site, uploads it as `numcosmo-site-<SHA>`, and restarts the Read the Docs
build, which downloads that artifact (`docs/download_docs_site.py`).

The check `docs/readthedocs.org:numcosmo` therefore fails on every push to a PR at
first: the Read the Docs webhook fires on the push and finds no artifact yet. It passes
once `docs-artifacts` finishes and restarts it, and Read the Docs then serves a preview
of the PR. Judge the documentation from the `docs-artifacts` job. Do not compare with
merged PRs, whose check passes because their artifact was uploaded earlier.
