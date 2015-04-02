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
for p in 51 55 59 115 179 ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-6-$p.png nice ./fluids_stress_2d 6 -parameter $p -dt .05 -resolution 8 -last_frame 1 6 -no_du -inv_Wi 0
done

$SLAVE -k
wait
