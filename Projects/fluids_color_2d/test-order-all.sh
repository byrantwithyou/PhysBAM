#!/bin/bash
rm -rf new_test_order
mkdir new_test_order
for t in 02 03 04 05 06 10 11 16 17 18 19 24 25 28; do
	echo $t
    for b in d n s ; do
        ./test-order.sh new_test_order/conv-$t-$b.png -bc_$b -dt .05 $t &
    done
done
for t in 00 01 07 08 09 12 13 14 20 21 ; do
    ./test-order.sh new_test_order/conv-$t-x.png -dt .05 $t &
done
wait
[ -d ref ] && [ -e ~/bin/src/compare ] && for i in `cd new_test_order ; ls conv*` ; do ~/bin/src/compare new_test_order/$i ref/$i x-$i ; done
times

