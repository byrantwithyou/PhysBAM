#!/bin/bash
#extrap_3d -framerate 9.45 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 960 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null
extrap_3d -framerate 24 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 240 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null
#extrap_3d -framerate 9.47 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 960 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null

for d in fig-tet-corot ; do
    mkdir -p $d/frames
    for f in `( cd $d/data ; echo *txt )` ; do
        compos.pl $f < $d/data/$f > $d.tex
        latex $d.tex
        dvips $d.dvi
        echo ==== $f ====
        convert $d.ps $d/frames/${f/txt/png}
        rm -f $d.tex $d.dvi $d.ps
    done
done

