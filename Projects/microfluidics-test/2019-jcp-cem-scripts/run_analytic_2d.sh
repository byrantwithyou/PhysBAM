#!/bin/bash

NAME=ana-2d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -q -threads 6"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
LO=2
HI=16
ANA="-u 'u=4+2*x*x+x*y-2.2*y,v=1.2-2.2*x+.7*y-y*y-x*x-3*x*y' -p 'p=2+3*x-2*y'"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/$c-r$r -refine $r $ANA ../$c.txt;
        done
    done | xargs -P 4 -n 1 -d '\n' bash -c > /dev/null
fi

find $NAME -name log.txt | xargs grep l-inf | sed 's/\(.*\)<print>\([^<]*\).*/\1\2/g'
