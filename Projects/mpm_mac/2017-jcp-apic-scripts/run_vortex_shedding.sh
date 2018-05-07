#!/bin/bash

NAME=vortex-shedding

RES=192
ARGS="../mpm_mac 21 -last_frame 3000 -mu 0 -resolution $RES"

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    while true ; do
        #echo $ARGS -no_affine -o $NAME/pic
        echo $ARGS -affine -o $NAME/apic
        echo $ARGS -no_affine -flip 1 -o $NAME/flip
        break
    done | xargs -P 2 -n 1 -d '\n' bash -c
fi

while true ; do
    echo ../../fourier-apic/region-vort $NAME/apic 7 -2 8.8 2 > $NAME/apic.log
    echo ../../fourier-apic/region-vort $NAME/flip 4.6 -2 6.4 2 > $NAME/flip.log
done | xargs -P 2 -n 1 -d '\n' bash -c

for a in apic flip ; do
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+(\(.\+\))/\1/g" > $NAME/$a-vort
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+: \([^ (]\+\).\+/\1/g" > $NAME/$a-time
done

# MATLAB code to extract the frequency.
cp freq.m $NAME
(
    cd $NAME
    matlab -nosplash -nodesktop -r "freq('flip',600);exit"
    matlab -nosplash -nodesktop -r "freq('apic',600);exit"
)

