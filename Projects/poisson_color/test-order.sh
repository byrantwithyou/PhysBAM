#!/bin/bash
for t in {1..12}; do
    (
        T="o-$t"
        L=""
        for r in {16..128}; do
            (
                echo -n $r
                nice ./poisson_color -resolution $r $t | grep error | tail -n 1
            ) > $T-$r &
            L="$L $T-$r"
        done
        wait
        
        cat $L | sort -n > $T
        rm $L
        
        gnuplot -p -e "set terminal png ; set output 'conv-$t.png' ; plot '$T' u (log10(\$1)):(log10(\$3)) , -2*x , -2*x-1 , -x-2"
    ) &
done
wait
times

