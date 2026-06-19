#!/bin/sh
#
# Canonical formatter for docs/references.bib.
#
# bibtex-tidy is the authoritative formatter for this file (the bibliography is no
# longer hand- or JabRef-formatted in the repo). This script is the single source
# of truth for the options, shared by developers, the pre-commit hook, and CI
# (.github/workflows/bib_lint.yml).
#
#   docs/tidy_references.sh            tidy the file in place
#   docs/tidy_references.sh --check    fail (non-zero) if the file is not tidy
#
# Requires bibtex-tidy (https://github.com/FlamingTempura/bibtex-tidy):
#   npm install -g bibtex-tidy@1.14.0
# Pin the version: different bibtex-tidy releases can format differently, which
# would defeat the point of a canonical form.

set -eu

DIR=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
BIB="$DIR/references.bib"

MODE=tidy
[ "${1:-}" = "--check" ] && MODE=check

# Keep these options in sync with the version pinned in CI.
set -- \
    --omit=file,owner,timestamp,__markedentry \
    --curly --numeric --align=14 \
    --sort=key \
    --duplicates=key,doi --merge=combine \
    --no-escape --sort-fields --trailing-commas

# bibtex-tidy requires the input file to precede --modify.
if [ "$MODE" = check ]; then
    TMPDIR_=$(mktemp -d)
    trap 'rm -rf "$TMPDIR_"' EXIT
    cp "$BIB" "$TMPDIR_/references.bib"
    bibtex-tidy "$TMPDIR_/references.bib" "$@" --modify --quiet >/dev/null 2>&1
    if diff -u "$BIB" "$TMPDIR_/references.bib"; then
        echo "references.bib is tidy."
    else
        echo "ERROR: references.bib is not tidy. Run: docs/tidy_references.sh" >&2
        exit 1
    fi
else
    bibtex-tidy "$BIB" "$@" --modify
    echo "references.bib tidied."
fi
