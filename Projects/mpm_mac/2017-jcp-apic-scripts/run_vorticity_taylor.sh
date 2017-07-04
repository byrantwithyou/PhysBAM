#!/bin/bash

NAME=vort-taylor

RES=64
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

dt=`perl -e "print (1.0/$RES)"`
DT="-max_dt $dt -min_dt $dt"
ARGS="../mpm_mac 16 -last_frame 10 -frame_dt 1 -mu 0 -scale_mass 3 -resolution $RES -analyze_u_modes -dump_modes_freq $RES $DT -max_ke 15"

opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
opt_name=("default" "regular" "default-cubic" "regular-cubic")

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        echo $ARGS $DT -no_affine ${opt[$o]} -o $NAME/pic-${opt_name[$o]}
        echo $ARGS $DT -affine ${opt[$o]} -o $NAME/apic-${opt_name[$o]}
        echo $ARGS $DT -no_affine ${opt[$o]} -flip 1 -o $NAME/flip-${opt_name[$o]}
    done | xargs -P 8 -n 1 -d '\n' bash -c
fi

for s in pic apic flip ; do
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        ./vort_parse_data.pl < $NAME/$s-${opt_name[$o]}/common/log.txt > $NAME/$s-${opt_name[$o]}-norms.txt
    done
done

for o in ${opt_name[@]} ; do
    sed -e "s/xxx/$o/" vort_taylor_plot.tex  > $NAME/plot-taylor-vort-$o.tex
done

for i in {apic,flip}-{default,regular} ; do
    for j in {0..10} ; do
        convert -scale 512 $NAME/$i/modes-$(($j*$RES))-x.png `printf $NAME/$i-taylor-vort-%02d.png $j`
        convert -scale 512 $NAME/$i/modes-$(($j*$RES))-x.png `printf $NAME/$i-taylor-vort-%02d.eps $j`
    done
done

cat <<EOF > $NAME/SConstruct
import os
import re
env=Environment(ENV = os.environ)
env['PSSUFFIX']=".eps"
r=re.compile(".*\.tex$")
for f in [x for x in os.listdir(".") if r.match(x)]:
    t=env.DVI(f)
    env.PostScript(t)
    env.PDF(t)
EOF
(
    cd $NAME
    scons -j 8
)
