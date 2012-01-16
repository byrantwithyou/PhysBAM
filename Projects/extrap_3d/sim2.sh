#!/bin/bash
#extrap_3d -framerate 9.45 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 240 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null
#extrap_3d -framerate 25 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 240 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null
#extrap_3d -framerate 9.47 -parameter 10 -stiffen .2 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 0.01 -last_frame 240 -poissons_ratio .45 -cutoff .001 -ether_drag 0 -o fig-tet-corot -use_corotated >&/dev/null

for d in fig-tet-corot ; do
    mkdir -p $d/frames-white
#    for f in `( cd $d/data ; echo *txt )` ; do
    for f in 179.txt ; do
#        compos2.pl $f < $d/data/$f > $d.tex
        compos2.pl -$f < $d/data/$f > $d.tex
        latex $d.tex
        dvips $d.dvi
        echo ==== $f ====
#        convert $d.ps $d/frames-white/${f/txt/png}
        convert $d.ps $d/frames-white/x${f/txt/png}
        rm -f $d.tex $d.dvi $d.ps
    done
done

