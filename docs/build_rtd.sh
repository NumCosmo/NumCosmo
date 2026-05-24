#!/usr/bin/env bash

set -e

echo "Repo=$READTHEDOCS_GITHUB_REPOSITORY"
echo "SHA=$READTHEDOCS_GIT_COMMIT_HASH"

echo "Testing GH artifact API access..."
curl -i \
  https://api.github.com/repos/NumCosmo/NumCosmo/actions/artifacts \
  || true

source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file environment.yml --prune
conda activate numcosmo_developer
meson setup build -Ddocumentation=true -Db_lto=false
python docs/download_quarto_freeze.py || true
meson compile -C build
