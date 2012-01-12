#!/bin/bash

function work()
{
    d=$1
    shift
    name=`echo "$@" | sed 's/\.txt./-/g' | sed 's/txt$/png/'`
    eps-cap/compos3.pl $d/data "$@" > $d.tex
    latex $d.tex >/dev/null ; dvips $d.dvi ; convert $d.ps $d/frames/$name
#    rm -f $d.tex $d.dvi $d.ps
}

if false ; then
    extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -use_corotated -o fig-corot-with-inv -parameter 2

    for d in fig-corot-with-inv ; do
        mkdir -p $d/frames
        for f in `( cd $d/data ; echo *txt )` ; do
            work $d/data/$f
        done
    done
else
    work fig-corot-with-inv 001.txt 020.txt 100.txt 200.txt
fi
