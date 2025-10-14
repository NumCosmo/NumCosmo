#!/usr/bin/env bash

set -e
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate numcosmo_developer
meson setup build -Ddocumentation=true -Db_lto=false
meson compile -C build
