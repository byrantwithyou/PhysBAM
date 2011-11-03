#!/bin/bash

declare -i i=1

dir=${1%/}

if [ -d $dir ]
then
    while [ -f $dir/SV_$i ]
    do
        echo "set size square; set terminal gif nocrop enhanced font verdana 12 size 800,800; set output \""$dir/plot_$i".gif\"; plot[-1:10][-5:5] 0.8/x,0.6/x,0.4/x,0.2/x, \""$dir"/SV_"$i"\"" | gnuplot
        i=i+1
    done
else
    echo \"$1\" does not name a directory
fi