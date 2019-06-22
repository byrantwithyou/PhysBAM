#!/bin/bash

NAME=illus-domain

ARGS="../fem -q -illus"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
REFINE=4

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        $ARGS -o $NAME/$c -refine $REFINE ../$c.txt;
        epstopdf $NAME/$c/domain.eps --outfile $NAME/$c-domain.pdf
    done
fi

../fem -q ../simple.txt -illus_meshing -min_corner "0.73 -0.105" -max_corner "0.78 -0.05" -refine 8 -o $NAME/simple-meshing
epstopdf $NAME/simple-meshing/meshing.eps --outfile $NAME/simple-meshing.pdf
epstopdf $NAME/simple-meshing/domain-anno.eps --outfile $NAME/simple-domain-anno.pdf

