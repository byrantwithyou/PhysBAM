#!/bin/bash

NAME=vortex-shedding

RES=192
DT=0.01
PPC=4
ARGS="../mpm_mac 21 -last_frame 3000 -mu 0 -resolution $RES -lsproj -no_surface -mls_zero_invalid -reseed -min_dt $DT -max_dt $DT -particles_per_cell $PPC"

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    while true ; do
        echo $ARGS -no_affine -xpic 2 -o $NAME/xpic2
        echo $ARGS -no_affine -xpic 3 -o $NAME/xpic3
        echo $ARGS -no_affine -xpic 5 -o $NAME/xpic5
        echo $ARGS -affine -o $NAME/apic
        break
    done | xargs -P 2 -n 1 -d '\n' bash -c
fi

../../fourier-apic/region-vort $NAME/apic 4.6 -2 6 2 > $NAME/apic.log
../../fourier-apic/region-vort $NAME/xpic2 4.6 -2 6 2 > $NAME/xpic2.log
../../fourier-apic/region-vort $NAME/xpic3 4.6 -2 6 2 > $NAME/xpic3.log
../../fourier-apic/region-vort $NAME/xpic5 4.6 -2 6 2 > $NAME/xpic5.log

for a in apic xpic2 xpic3 xpic5 ; do
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+(\(.\+\))/\1/g" > $NAME/$a-vort
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+: \([^ (]\+\).\+/\1/g" > $NAME/$a-time
done

# Octave code to extract the frequency.
cp freq.m $NAME
(
    cd $NAME
    octave --eval "freq('apic',2000);exit" > apic-fft
    octave --eval "freq('xpic2',2000);exit" > xpic2-fft
    octave --eval "freq('xpic3',2000);exit" > xpic3-fft
    octave --eval "freq('xpic5',2000);exit" > xpic5-fft
)

MAXVORT=8
frames=(1700 1706 1708 1691)
algos=(apic xpic2 xpic3 xpic5)
O=$(readlink -m '../dump_vort')
(
cd $NAME
for a in `seq 0 $((${#algos[@]}-1))` ; do
    x=${algos[$a]}
    f=${frames[$a]}
    $O $NAME/$x -frame $f -vmin 0 -vmax $MAXVORT -o vortex-shedding-$x-$f;
done
)

for a in `seq 0 $((${#algos[@]}-1))` ; do
    x=${algos[$a]}
    f=${frames[$a]}
    sed -e "s/XXXXXX/vortex-shedding-$x-$f.png/g" vortex_shedding_window_vorticity.tex > $NAME/vortex-shedding-$x.tex
done

cat <<EOF > $NAME/SConstruct
import os
import re
env=Environment(ENV = os.environ)
env['PSSUFFIX']=".eps"
r=re.compile(".*\.tex$")
for f in [x for x in os.listdir(".") if r.match(x)]:
    env.PDF(f)
EOF
(
    cd $NAME
    scons -j 8
)
