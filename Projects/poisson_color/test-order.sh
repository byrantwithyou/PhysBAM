#!/bin/bash
for t in {10..10}; do
    (
        T="o-$t"
        L=""
        for r in {40..160}; do
            (
                echo -n $r
                nice ./poisson_color -resolution $r $t -use_preconditioner -double_fine | grep error | tail -n 1
            ) > $T-$r ;
            L="$L $T-$r"
        done
        wait
        
        cat $L | sort -n > $T
        rm $L
        
        gnuplot -p -e "set terminal png ; set output 'conv-$t.png' ; plot '$T' u (log10(\$1)):(log10(\$3)) , -2*x , -2*x-1 , -x-2"
    ) ;
done
wait
times

