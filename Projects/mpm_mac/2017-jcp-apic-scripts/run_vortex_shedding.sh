#!/bin/bash

NAME=vortex-shedding

RES=192
FRAME=1600
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

../../fourier-apic/region-vort $NAME/apic 4.6 -2 6 2 > $NAME/apic.log &
../../fourier-apic/region-vort $NAME/flip 4.6 -2 6 2 > $NAME/flip.log &
wait

for a in apic flip ; do
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+(\(.\+\))/\1/g" > $NAME/$a-vort
    grep -o "window vorticity:[^<]\+" $NAME/$a.log | sed "s/.\+: \([^ (]\+\).\+/\1/g" > $NAME/$a-time
done

# MATLAB code to extract the frequency.
cp freq.m $NAME
cp data_zeros.m $NAME
cp gen_period.m $NAME
(
    cd $NAME
    matlab -nosplash -nodesktop -r "freq('flip',1000);exit" > flip-fft
    matlab -nosplash -nojvm -r "gen_period('flip',1000);exit"
    matlab -nosplash -nodesktop -r "freq('apic',1000);exit" > apic-fft
    matlab -nosplash -nojvm -r "gen_period('apic',1000);exit"
)

flip_freq=`grep 'Period' $NAME/flip-fft | sed "s/Period: \(.\+\)./\1/g"`
apic_freq=`grep 'Period' $NAME/apic-fft | sed "s/Period: \(.\+\)./\1/g"`
sed -e "s/AAAAAA/$apic_freq/g; s/FFFFFF/$flip_freq/g" vortex_shedding_periods.tex  > $NAME/vortex-shedding-periods.tex

(
cd $NAME
for a in apic flip ; do
    ../../../opengl_2d/opengl_2d $a -offscreen -start_frame $FRAME -stop_frame $FRAME -so vortex-shedding-$a-$FRAME.png -keys "m6sVvsNs"
    convert -crop 1024x512+0+128 vortex-shedding-$a-$FRAME.png vortex-shedding-$a-$FRAME.png
done
)

sed -e "s/XXXXXX/vortex-shedding-apic-1600.png/g" vortex_shedding_window_vorticity.tex > $NAME/vortex-shedding-apic.tex
sed -e "s/XXXXXX/vortex-shedding-flip-1600.png/g" vortex_shedding_window_vorticity.tex > $NAME/vortex-shedding-flip.tex

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
