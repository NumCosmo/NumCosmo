#!/bin/bash
# Author:Andrey Nikishaev, Gunnar Lindholm
#
# Writes the ChangeLog to stdout: the commits between consecutive version tags (v*, newest
# first), then the subjects of the commits up to the oldest tag. With a version argument
# (bash generate-changelog.sh 0.29.0) the commits since the last tag are listed under
# [v0.29.0] below an empty [Current]; without it they are listed under [Current]. Carriage
# returns in old commit messages are removed.
set -eu

NEW_VERSION=${1:-}

echo "CHANGELOG"
echo ----------------------
echo
echo "[Current]"
if [ -n "$NEW_VERSION" ]; then
    echo
    echo "[v$NEW_VERSION]"
fi

NEXT=HEAD
for TAG in $(git tag --list 'v*' --sort=-v:refname); do
    if [ "$NEXT" != HEAD ]; then
        echo
        echo "[$NEXT]"
    fi
    GIT_PAGER=cat git log --no-merges --format=" *%w(80,1,5)%B" "$TAG..$NEXT" | tr -d '\r'
    NEXT=$TAG
done

echo
echo "[$NEXT]"
GIT_PAGER=cat git log --no-merges --format=" *%w(80,1,5)%s" "$NEXT" | tr -d '\r'
