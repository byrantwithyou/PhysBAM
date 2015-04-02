#!/bin/bash
MASTER="$PHYSBAM/Tools/batch/master"
SLAVE="$PHYSBAM/Tools/batch/slave"

max=20
if [ "x$1" = "x-m" ] ; then
    max=$2
    shift 2
fi

$MASTER &

rm -rf new_test_order
mkdir new_test_order
for t in 00 01 02 03 04 05 ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-$t.png nice ./fluids_stress_2d -dt .05 -resolution 8 -last_frame 1 $t
done
for p in 1 3 7 15 31 63 127 255 19 51 187 ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-6-$p.png nice ./fluids_stress_2d 6 -parameter $p -dt .05 -resolution 8 -last_frame 1 6
done

$SLAVE -k
wait
