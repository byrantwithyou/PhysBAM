#!/bin/bash
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -use_corotated -o fig-corot >&/dev/null &
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -o fig-neo >&/dev/null &
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -use_corotated -o fig-corot-inv -parameter 1 >&/dev/null &
extrap_2d -stiffen .1 -dampen 0 31 -image_size 30 -sigma_range 3 -dt 1 -last_frame 250 -poissons_ratio .45 -cutoff .001 -ether_drag 5 -use_corotated -o fig-stretch >&/dev/null &
wait

for d in fig-* ; do
    mkdir -p $d/frames
    for f in `( cd $d/data ; echo *txt )` ; do
        eps-cap/compos.pl < $d/data/$f > $d.tex
        latex $d.tex >/dev/null ; dvips $d.dvi ; convert $d.ps $d/frames/${f/txt/png}
        rm -f $d.tex $d.dvi $d.ps
    done &
done
wait
