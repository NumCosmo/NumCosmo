#!/usr/bin/env bash

set -e

source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file environment.yml --prune
conda activate numcosmo_developer
meson setup build -Ddocumentation=true -Db_lto=false
python docs/download_quarto_freeze.py || true
meson compile -C build
