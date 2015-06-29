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
        ./test-order-one.sh -b 2 $max new_test_order/conv-$t-$b.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -bc_$b -dt .05 $t
    done
done
for t in 00 01 08 09 12 13 14 20 21  250 251 252 253 254 255 256 257 258 259  ; do
    ./test-order-one.sh -b 2 $max new_test_order/conv-$t.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 $t
done

U=("u=1;v=1;" "u=t+1;v=t+1;")
P=("p=0;")
S=("s00=1;s10=0;s11=1;" "s00=x*x;s10=x*y;s11=y*y;")
for u in {0..1} ; do
    for p in 0 ; do
        for s in {0..1} ; do
            ./test-order-one.sh -b 2 $max new_test_order/conv-custom-$u-$p-$s.png nice ./fluids_color -resolution 8 -last_frame 1 -s 1.3 -m .8 -kg 1.2 -dt .05 250 -set_u0 "${U[$u]}" -set_S0 "${S[$s]}" -set_p0 "${P[$p]}"
        done
    done
done

$SLAVE -k
wait
