#!/bin/bash

# test-order.sh output args
out="$1"
shift
L=""
for r in {1..9} ; do
    T=`mktemp`
    O=`mktemp -d`
    (
        echo -n "$((8*$r)) "
        nice ./fluids_color_2d -resolution 8 -last_frame 1 -refine $r -s 1.3 -m .8 -kg 1.2 -o $O "$@"  | grep error | tail -n 1
        rm -rf $O
    ) > $T &
    L="$L $T"
done
wait

T=`mktemp`
cat $L | sort -n > $T
gnuplot -p -e "set terminal png ; set output '$out' ; plot '$T' u (log10(\$1)):(log10(\$3)) title 'L-inf error' , -2*x , -2*x-1 , -x-2"
rm $T $L
