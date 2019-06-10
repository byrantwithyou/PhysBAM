#!/bin/bash

NAME=ana-3d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -3d -q -threads 16"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
LO=2
HI=5
ANA="-u 'u=4+2*x*x+x*y-2.2*y+z*x,v=1.2-2.2*x+.7*y-y*y-z*x-3*x*y,w=1+3*x*y+z*y' -p 'p=1+3*x-2*y+z'"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/$c-r$r -refine $r $ANA ../$c.txt;
        done
    done | xargs -P 1 -n 1 -d '\n' bash -c > /dev/null
fi

find $NAME -name log.txt | xargs grep l-inf | sed 's/\(.*\)<print>\([^<]*\).*/\1\2/g'
