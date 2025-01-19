## Creating New Releases

The following steps should be followed to create a new release of NumCosmo:

- Update the version number in `meson.build`, `pyproject.toml`, and `tests/test_py_cfg.py`.
- Update the changelog using the command using the shell, e.g., `bash generate-changelog.sh > ChangeLog.md`.
- Adjust ChangeLog.md so to add \[vX.X.X\] below \[Current\] to reflect the new version.
- Commit changes to `meson.build`, `pyproject.toml`, `tests/test_py_cfg.py`, and `ChangeLog.md`.
- Create a new tag, e.g., `git tag -a -s v0.24.0 -m "New version v0.24.0."`.
- Push changes to GitHub, e.g., `git push --follow-tags origin main`.
- Create a new release on GitHub, e.g., https://github.com/numcosmo/numcosmo/releases/new
  - Create the release package, for example, if the build directory is `build`, use `meson dist -C build --no-tests`
  - Upload the release package to GitHub.
  