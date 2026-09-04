#!/usr/bin/env bash
#
# Documentation style check for comments, docstrings and .qmd pages.
#
# Three checks, described in CONTRIBUTING.md, section "Writing style for
# comments and documentation":
#
#   1. Blocking. Evaluative or promotional wording in lines this branch adds.
#   2. Blocking. Narration of the code's own history in lines this branch adds.
#      A comment states what the code does now; how it came to do that goes in
#      dev-notes/ or an area history document.
#   3. Advisory. Theory-page links (<a href="../../theory/...">, [[numcosmo|Sym]])
#      that this branch removes without adding an equivalent back. The pipe may
#      be escaped as \| , which is how a link is written inside a table. Reported, not
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
    -h|--help) sed -n '2,19p' "$0"; exit 0 ;;
    *) echo "unknown argument: $1" >&2; exit 2 ;;
  esac
done

cd "$(git rev-parse --show-toplevel)"

# Case-insensitive, word-boundary where it matters. Kept as one alternation so
# the reported line shows which term matched.
WORDS='shipped|blows? up|blowup|comfortably|blindly|genuinely|decisive|safe zone|sweet spot|nets out|load-bearing|at length|worth knowing|needlessly|wins big|obvious(ly)?|hot loop|of course|magic|\btraps?\b|\bstory\b|honest(ly)?|clever(ly)?|elegant|it turns out|leverag(e|es|ed|ing)|freshly|simply|carefully|\bfir(e|es|ed|ing)\b'

# Wording that narrates the code's history instead of stating its behaviour.
# Kept separate from WORDS so the message can point at the right fix.
#
# Two exclusions. "used to" and "no longer" on their own are not listed:
# "the buffer used to hold the result" and "freed when no longer needed" are
# the ordinary way to say those things. And a past *interface* is sometimes the
# point -- a migration note saying which property an old serialized file set,
# and what it now maps to, has to state that -- so only past *defects* are
# matched, which no comment needs to recount.
HISTORY='historically|in the past|originally|nowadays|these days|at the time|the old (code|path|version|implementation|way)|used to (abort|crash|hang|fail|refuse|stall|deadlock|leak|be broken|be wrong)'

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

if [ "$MODE" != "all" ]; then
  merge_base=$(git merge-base HEAD "$BASE" 2>/dev/null || echo "")
  if [ -z "$merge_base" ]; then
    echo "check_doc_style: cannot resolve $BASE, checking the working tree instead" >&2
    merge_base="HEAD"
  fi
fi

# Lines matching $1: the whole tree in --all mode, this branch's added lines
# otherwise, with the file name carried down from the +++ header.
scan_lines () {
  if [ "$MODE" = "all" ]; then
    git grep -nEI "$1" -- "${paths[@]}" | grep -Ev "$ALLOW" || true
  else
    git diff -U0 "$merge_base" -- "${paths[@]}" | awk -v words="$1" '
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
    ' | grep -Ev "$ALLOW" || true
  fi
}

# ---------------------------------------------------------------- word check
hits=$(scan_lines "$WORDS")

if [ -n "$hits" ]; then
  echo "Evaluative wording found (CONTRIBUTING.md, \"Writing style for comments and documentation\"):"
  echo ""
  echo "$hits" | sed 's/^/  /'
  echo ""
  echo "State the fact instead. If a term is genuinely the right technical word here,"
  echo "adjust WORDS in .github/scripts/check_doc_style.sh and say why in the commit message."
  status=1
fi

# ------------------------------------------------------------- history check
hits=$(scan_lines "$HISTORY")

if [ -n "$hits" ]; then
  echo "History narration found (CONTRIBUTING.md, \"Writing style for comments and documentation\"):"
  echo ""
  echo "$hits" | sed 's/^/  /'
  echo ""
  echo "State the current behaviour. What the code did before, and why it changed, goes"
  echo "in dev-notes/<topic>.md or an area history document, with a pointer left here."
  status=1
fi

# ---------------------------------------------------------------- link check
if [ "$MODE" != "all" ]; then
  removed=$(git diff -U0 "$merge_base" -- "${paths[@]}" \
    | grep '^-' | grep -v '^---' \
    | grep -oE 'href="\.\./\.\./theory/[^"]+"|\[\[numcosmo[a-z-]*\\?\|[^]]+\]\]' \
    | sed 's/\\|/|/g' \
    | sort | uniq -c | sed 's/^ *//' || true)
  added=$(git diff -U0 "$merge_base" -- "${paths[@]}" \
    | grep '^+' | grep -v '^+++' \
    | grep -oE 'href="\.\./\.\./theory/[^"]+"|\[\[numcosmo[a-z-]*\\?\|[^]]+\]\]' \
    | sed 's/\\|/|/g' \
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
