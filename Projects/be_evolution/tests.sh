#!/bin/bash

OUT_LIST=""
rm -f tmp-sim-list.txt
for r in 6 11 16 21 ; do
    for t in 1 2 3 ; do
        OUT_LIST="$OUT_LIST out-$t-$r"
        echo "( time ./be_evolution -resolution $r $t -newton_cd_tol 1e-3 -o out-$t-$r $@ ) >& /dev/null" >> tmp-sim-list.txt
    done
done
xargs -n 1 -d '\n' -P 12 bash -c < tmp-sim-list.txt
rm -f tmp-sim-list.txt

./make-table.pl $OUT_LIST
