#!/bin/sh
# Regenerate the introspection stubs. Safe to run from any directory.
set -eu

cd "$(dirname "$0")"

gen() {
    python ./generate_stubs.py "$1" 1.0 > "$2.tmp"
    # A failed generation yields an empty file; never overwrite a good stub with it.
    [ -s "$2.tmp" ] || { echo "generate_stubs.py $1: no output" >&2; rm -f "$2.tmp"; exit 1; }
    mv "$2.tmp" "$2"
}

gen NumCosmoMath ncm.pyi
gen NumCosmo nc.pyi

# Only the generated stubs: "black ." here would reformat whatever tree it is run from.
black --quiet --target-version py314 ncm.pyi nc.pyi
