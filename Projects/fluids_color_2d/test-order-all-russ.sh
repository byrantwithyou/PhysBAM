#!/bin/bash
rm -rf new_test_order
mkdir new_test_order
maxjobs=9 #Results may vary on your machine.

#set -x verbose #echo on


parallelize () {
    while [ $# -gt 0 ]; do
        count=(`jobs -p`)
        if [ ${#count[@]} -lt $maxjobs ]; then
    		for b in d n s; do
        		./test-order-old.sh new_test_order/conv-$1-$b.png -bc_$b -dt .05 $1 &
        		#shift
    		done
    		shift

			fi
			done
			wait

}


#Numbers of tests with multiple boundary conditions
parallelize 02 03 04 05 06 10 11 16 17 18 19 24 25 28

for t in 00 01 08 09 12 13 14 20 21 ; do
    ./test-order-old.sh new_test_order/conv-$t-x.png -dt .05 $t &
done
wait
[ -d ref ] && [ -e ~/bin/src/compare ] && for i in `cd new_test_order ; ls conv*` ; do ~/bin/src/compare new_test_order/$i ref/$i x-$i ; done
times

