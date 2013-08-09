#!/bin/bash

OUT_LIST=""
(
    for r in 6 11 16 21 ; do
        for t in 1 2 3 ; do
            O="out-$t-$r.txt"
            OUT_LIST="$OUT_LIST $O"
            echo "( time ./be_evolution -resolution $r $t -newton_cd_tol 1e-3 -o out-$t-$r $@ ) >& $O"
        done
    done
) | xargs -n 1 -d '\n' -P 12 bash -c

./make-table.pl $OUT_LIST
