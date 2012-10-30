#!/bin/bash
for t in {1..11}; do
    (
        T="o-$t"
        L=""
        for r in 16 20 24 28 32 36 40 44 48 52 56 60 64 68 72 76 80 84 88 92 96 100; do
            (
                echo -n $r
                nice ./poisson_color -resolution $r $t -use_preconditioner | grep error | tail -n 1
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

