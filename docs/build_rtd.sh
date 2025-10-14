#!/usr/bin/env bash

set -e
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env update --name numcosmo_developer --file devel_environment.yml --prune
conda activate numcosmo_developer
wget https://github.com/quarto-dev/quarto-cli/releases/download/v1.6.42/quarto-1.6.42-linux-amd64.tar.gz && mkdir ~/opt && tar -C ~/opt -xvzf quarto-1.6.42-linux-amd64.tar.gz
export PATH="$HOME/opt/quarto-1.6.42/bin:$PATH"
meson setup build -Ddocumentation=true -Db_lto=false
meson compile -C build
