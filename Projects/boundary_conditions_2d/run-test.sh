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
cat $T
gnuplot <<EOF
set terminal postscript color eps
set output 'plot.eps'
set size 1,1
plot '$T' u (log10(\$1)):(log10(\$3)) title 'L-inf p' , '$T' u (log10(\$1)):(log10(\$4)) title 'L-1 p' , '$T' u (log10(\$1)):(log10(\$6)) title 'L-inf u' , '$T' u (log10(\$1)):(log10(\$7)) title 'L-1 p' , -2*x title 'order 2' ls 1 , -x title 'order 1' ls 2
EOF
inkscape -D -e plot2.png -h 1200 $PWD/plot.eps
convert plot2.png -flatten plot.png
rm -f plot2.png plot.eps $F $T
