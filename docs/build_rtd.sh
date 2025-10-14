#!/usr/bin/env bash

set -e
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file devel_environment.yml --prune
conda activate numcosmo_developer
meson setup build -Ddocumentation=true -Db_lto=false
ls build/docs
quarto render build/docs --output-dir test-site --execute-debug --log-level debug --no-cache --execute-daemon-restart --debug
echo "Quarto site built"
meson compile -C build
