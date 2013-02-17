#!/bin/bash

out="$1"
shift

M=`mktemp`

(
    for t in "$@" ; do
        echo -n `head -n 1 $t`
        grep error $t | tail -n 1
    done
) | sort -n > M

gnuplot -p -e "set terminal png ; set output '$out' ; plot '$M' u (log10(\$1)):(log10(\$3)) title 'L-inf error' , -2*x , -2*x-1 , -x-2;"
rm "$@" $M
