#!/usr/bin/env bash

set -e

source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file environment.yml --prune
conda activate numcosmo_developer
meson setup build -Ddocumentation=true -Db_lto=false
echo "=== Running freeze downloader ==="
python -u docs/download_quarto_freeze.py || echo "Downloader failed"
echo "=== Freeze downloader finished ==="
meson compile -C build
