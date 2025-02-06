#!/bin/bash
# Author:Andrey Nikishaev, Gunnar Lindholm
echo "CHANGELOG"
echo ----------------------
git for-each-ref --format='%(taggerdate:raw) %(tag)' refs/tags | sort -n | awk '{print $3}' | tac |grep -v '^$' | while read TAG ; do
     echo
    if [ $NEXT ];then
        echo [$NEXT]
    else
        echo "[Current]"
    fi
    GIT_PAGER=cat git log --no-merges --format=" *%w(80,1,5)%B" $TAG..$NEXT
    NEXT=$TAG
done
FIRST=$(git for-each-ref --format='%(taggerdate:raw) %(tag)' refs/tags | sort -n | awk '{print $3}' | head -1)
echo
echo [$FIRST]
GIT_PAGER=cat git log --no-merges --format=" *%w(80,1,5)%s" $FIRST
