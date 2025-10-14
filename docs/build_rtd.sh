#!/usr/bin/env bash

set -e
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate base
meson setup build -Ddocumentation=true -Db_lto=false
meson compile -C build
