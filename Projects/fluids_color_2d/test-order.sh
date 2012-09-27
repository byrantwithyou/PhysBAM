#!/bin/bash
for t in {0..12} ; do
    for b in d n s ; do
        (
            T="o-$t-$b"
            L=""
            for r in {1..9} ; do
                (
                    echo -n "$((8*$r)) "
                    nice ./fluids_color_2d -bc_$b -resolution 8 $t -last_frame 1 -dt .05 -steps 1 -mu0 1 -refine $r | grep error | tail -n 1
                ) > $T-$r &
                L="$L $T-$r"
            done
            wait

            cat $L | sort -n > $T
            rm $L

            gnuplot -p -e "set terminal png ; set output 'conv-$t-$b.png' ; plot '$T' u (log10(\$1)):(log10(\$3)) , -2*x , -2*x-1 , -x-2"
        ) &
    done
done
wait
times

