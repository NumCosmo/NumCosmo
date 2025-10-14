#!/usr/bin/env bash

set -e
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file devel_environment.yml --prune
conda activate numcosmo_developer
export JUPYTER_RUNTIME_DIR=$(pwd)/.jupyter_runtime
mkdir -p "$JUPYTER_RUNTIME_DIR"
meson setup build -Ddocumentation=true -Db_lto=false
meson compile -C build
