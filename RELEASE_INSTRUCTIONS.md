# Releasing NumCosmo

The steps below release version `X.Y.Z`.

## 1. Check master

- CI on master is green and every open pull request is merged or left for the next release.
- Every data file the tests use is on the `datafile-release-v1.0.0` release.
  `python3 .github/scripts/prefetch_data.py` prints their paths; it skips files already in the
  local data directory, so compare the names with
  `gh release view datafile-release-v1.0.0 --json assets --jq '.assets[].name'`.
- `python3 .github/scripts/check_deps.py feedstock` passes. Packages it reports go into the
  conda-forge recipe in step 4.

## 2. Bump the version in a branch

```bash
git checkout master && git pull
git checkout -b release-X.Y.Z
```

Set the version in these three files:

| file | line |
|---|---|
| `meson.build` | `version: 'X.Y.Z',` |
| `pyproject.toml` | `version = 'X.Y.Z'` |
| `tests/python/ncm/core/test_cfg.py` | `VERSION_MAJOR`, `VERSION_MINOR`, `VERSION_MICRO` |

Everything else (the library soname, the C version macros, the documentation) is derived from
`meson.build`.

Regenerate the ChangeLog and add a `[vX.Y.Z]` line below `[Current]`, so the commits since the
last tag fall under the new version:

```bash
bash generate-changelog.sh > ChangeLog.md
```

Commit as `Release version X.Y.Z`, open a pull request and merge it once CI passes.

## 3. Tag and publish

Tag the merged commit on master and push the tag:

```bash
git checkout master && git pull
git tag -a -s vX.Y.Z -m "New version vX.Y.Z."
git push origin vX.Y.Z
```

Build the tarball from a configured build directory and publish the GitHub release with it:

```bash
meson dist -C build --no-tests
gh release create vX.Y.Z build/meson-dist/numcosmo-X.Y.Z.tar.xz --verify-tag --title "NumCosmo vX.Y.Z"
```

## 4. Update conda-forge

After the GitHub release, the conda-forge bot opens a version pull request on
`conda-forge/numcosmo-feedstock`. To get it at once, open an issue on the feedstock titled
`@conda-forge-admin, please update version`. In that pull request:

- the version and the `sha256` match the published tarball, and the build number is 0;
- the `run` requirements include the packages `check_deps.py feedstock` reported in step 1;
- the test import line in `recipe/meta.yaml` lists every `numcosmo_py` subpackage, and the test
  list in `recipe/build.sh` names only tests present in this version.

Merge it once its CI passes.
