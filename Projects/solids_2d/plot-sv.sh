#!/bin/bash

for f in {1..1000} ; do
    F=`printf "$1" $f`
    [ -e $F ] || break

    gnuplot <<EOF
set terminal png size 640,480
set output "plot-$f.png"
set xrange [$2:$3]
set yrange [$4:$5]
plot "$F"
EOF
done
ffmpeg -i plot-%d.png plot-sv.ogv

