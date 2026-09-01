#!/usr/bin/env bash
#
# Documentation style check for comments, docstrings and .qmd pages.
#
# Two checks, described in CONTRIBUTING.md, section "Writing style for comments
# and documentation":
#
#   1. Blocking. Evaluative or promotional wording in lines this branch adds.
#   2. Advisory. Theory-page links (<a href="../../theory/...">, [[numcosmo|Sym]])
#      that this branch removes without adding an equivalent back. Reported, not
#      fatal: a link legitimately disappears when its page or symbol does.
#
# Usage:
#   .github/scripts/check_doc_style.sh            # added lines vs the merge base
#   .github/scripts/check_doc_style.sh --all      # every tracked file
#   .github/scripts/check_doc_style.sh --base REF # diff against REF instead of origin/master

set -euo pipefail

BASE="origin/master"
MODE="diff"

while [ $# -gt 0 ]; do
  case "$1" in
    --all) MODE="all"; shift ;;
    --base) BASE="$2"; shift 2 ;;
    -h|--help) sed -n '2,18p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

cd "$(git rev-parse --show-toplevel)"

# Case-insensitive, word-boundary where it matters. Kept as one alternation so
# the reported line shows which term matched.
WORDS='shipped|blows? up|blowup|comfortably|blindly|genuinely|decisive|safe zone|sweet spot|nets out|load-bearing|at length|worth knowing|needlessly|wins big'

# Legitimate uses of the same words. "the tables shipped with NumCosmo" is
# distribution, not self-praise. A line carrying "doc-style: allow" is skipped
# outright, for the case a term really is the right technical word.
ALLOW='shipped (with|in) |doc-style: allow'

# Files whose prose this applies to. numcosmo/external is upstream code.
paths=(
  ':(glob)numcosmo/**/*.c' ':(glob)numcosmo/**/*.h'
  ':(glob)numcosmo_py/**/*.py' ':(glob)docs/**/*.qmd'
  ':(exclude)numcosmo/external/**'
)

status=0

# ---------------------------------------------------------------- word check
if [ "$MODE" = "all" ]; then
  hits=$(git grep -nEI "$WORDS" -- "${paths[@]}" | grep -Ev "$ALLOW" || true)
else
  merge_base=$(git merge-base HEAD "$BASE" 2>/dev/null || echo "")
  if [ -z "$merge_base" ]; then
    echo "check_doc_style: cannot resolve $BASE, checking the working tree instead" >&2
    merge_base="HEAD"
  fi
  # Added lines only, with the file name carried down from the +++ header.
  hits=$(git diff -U0 "$merge_base" -- "${paths[@]}" | awk -v words="$WORDS" '
    /^\+\+\+ b\// { file = substr($0, 7); next }
    /^@@/ {
      # @@ -a,b +c,d @@ -> next added line is at c
      split($3, p, ",")
      lineno = p[1] + 0
      if (lineno < 0) lineno = -lineno
      next
    }
    /^\+/ && !/^\+\+\+/ {
      body = substr($0, 2)
      if (tolower(body) ~ tolower(words)) printf "%s:%d:%s\n", file, lineno, body
      lineno++
      next
    }
    /^ / { lineno++ }
  ' | grep -Ev "$ALLOW" || true)
fi

if [ -n "$hits" ]; then
  echo "Evaluative wording found (CONTRIBUTING.md, \"Writing style for comments and documentation\"):"
  echo ""
  echo "$hits" | sed 's/^/  /'
  echo ""
  echo "State the fact instead. If a term is genuinely the right technical word here,"
  echo "adjust WORDS in .github/scripts/check_doc_style.sh and say why in the commit message."
  status=1
fi

# ---------------------------------------------------------------- link check
if [ "$MODE" != "all" ]; then
  removed=$(git diff -U0 "$merge_base" -- "${paths[@]}" \
    | grep '^-' | grep -v '^---' \
    | grep -oE 'href="\.\./\.\./theory/[^"]+"|\[\[numcosmo[a-z-]*\|[^]]+\]\]' \
    | sort | uniq -c | sed 's/^ *//' || true)
  added=$(git diff -U0 "$merge_base" -- "${paths[@]}" \
    | grep '^+' | grep -v '^+++' \
    | grep -oE 'href="\.\./\.\./theory/[^"]+"|\[\[numcosmo[a-z-]*\|[^]]+\]\]' \
    | sort | uniq -c | sed 's/^ *//' || true)

  lost=$(comm -23 \
    <(echo "$removed" | sed 's/^[0-9]* //' | sort -u) \
    <(echo "$added"   | sed 's/^[0-9]* //' | sort -u) || true)

  if [ -n "$lost" ]; then
    echo "Note: cross-references this branch removes and does not add back:"
    echo ""
    echo "$lost" | sed 's/^/  /'
    echo ""
    echo "Advisory only. Re-add them unless the target page or symbol is also gone."
  fi
fi

if [ "$status" -eq 0 ]; then
  echo "check_doc_style: clean"
fi
exit "$status"
