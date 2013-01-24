#!/bin/bash

F=""
for r in 8 12 16 24 32 48 64 96 128 192 256 384 512 ; do
    T=`tempfile`
    F="$F $T"
    echo $r > $T
    projection -resolution $r "$@" >> $T &
done
wait
T=`tempfile`
cat $F | grep -v 'iterations\|Simulation' | tr '\n' '@' | sed 's/@ *\([up]\)/ \1/g' | tr @ '\n' | sort -n > $T
E1=`tempfile --suffix=.eps`
P1=`tempfile --suffix=.png`
gnuplot <<EOF
set terminal postscript color eps
set output '$E1'
set logscale xy
plot '$T' u 1:3 title 'L-inf p' , '$T' u 1:4 title 'L-1 p' , '$T' u 1:6 title 'L-inf u' , '$T' u 1:7 title 'L-1 u' , x**(-2) title 'order 2' ls 1 , x**(-1) title 'order 1' ls 2
EOF
inkscape -D -e $P1 -h 1200 $E1
convert $P1 -flatten plot-`echo "$@" | sed 's/ /-/g'`.png
rm -f $P1 $E1 $F $T
