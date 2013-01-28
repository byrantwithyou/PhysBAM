#!/bin/bash
rm -rf new_test_order
mkdir new_test_order
for t in {00..14} {16..21} ; do
    for b in d n s ; do
        (
            T="o-$t-$b"
            L=""
            for r in {1..9} ; do
                (
                    echo -n "$((8*$r)) "
                    nice ./fluids_color_2d -bc_$b -resolution 8 $t -last_frame 1 -dt .05 -steps 1 -mu0 1 -refine $r -o out-$t-$b-$r  | grep error | tail -n 1
                    rm -rf out-$t-$b-$r
                ) > $T-$r &
                L="$L $T-$r"
            done
            wait

            cat $L | sort -n > $T
            rm $L

            gnuplot -p -e "set terminal png ; set output 'new_test_order/conv-$t-$b.png' ; plot '$T' u (log10(\$1)):(log10(\$3)) , -2*x , -2*x-1 , -x-2"
            rm $T
        ) &
    done
done
wait
times

