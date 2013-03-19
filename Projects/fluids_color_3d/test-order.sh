#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

max=9
if [ "x$1" = "x-m" ] ; then
    max=$2
    shift 2
fi

$MASTER &

rm -rf new_test_order
mkdir new_test_order
for t in 02 03 04 05 06 10 11 16 17 18 19 24 28 ; do
    for b in d n s ; do
        ./test-order-one.sh -b 2 $max new_test_order/conv-$t-$b.png nice ./fluids_color_2d -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -bc_$b -dt .05 $t
    done
done
for t in 00 01 08 09 12 13 14 20 21 ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-$t-$b.png nice ./fluids_color_2d -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 $t
done

$SLAVE -k
wait
