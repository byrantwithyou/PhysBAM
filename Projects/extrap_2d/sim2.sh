#!/bin/bash
extrap_2d -framerate 1 -parameter 10 -stiffen 10 -dampen 0 30 -image_size 30 -sigma_range 3 -dt .01 -last_frame 200 -poissons_ratio .45 -cutoff .001 -ether_drag 0.01 -use_corotated -o fig-mattress >&/dev/null

for d in fig-mattress ; do
    mkdir -p $d/frames
    for f in `( cd $d/data ; echo *txt )` ; do
        eps-cap/compos2.pl < $d/data/$f > $d.tex
        latex $d.tex
        dvips $d.dvi
        echo ==== $f ====
        convert $d.ps $d/frames/${f/txt/png}
        rm -f $d.tex $d.dvi $d.ps
    done
done

