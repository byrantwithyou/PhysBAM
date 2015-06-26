#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

#Doing the visco tests separately so it doesn't take forever

max=9
if [ "x$1" = "x-m" ] ; then
    max=$2
    shift 2
fi

$MASTER &

rm -rf new_test_v
mkdir new_test_v
for t in 252 255 256 257; do
    for b in d n s ; do
        ./test-order-one.sh -b 2 $max new_test_v/conv-$t-$b.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1 -m 1 -kg 1 -bc_$b -dt .05 $t
    done
done
for t in 00 ; do
    ./test-order-one.sh -b 2 $max new_test_v/conv-$t-$b.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1 -m 1 -kg 1 -dt .05 $t
done

$SLAVE -k
wait
