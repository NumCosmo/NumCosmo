#!/usr/bin/env bash
set -e

find numcosmo-site \
  -name "*.ipynb" \
  -print0 | \
while IFS= read -r -d '' nb; do
    jupyter nbconvert --to script "$nb"
done
