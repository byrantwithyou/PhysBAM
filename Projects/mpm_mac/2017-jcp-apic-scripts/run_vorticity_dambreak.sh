#!/bin/bash

NAME=vort-dambreak

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../mpm_mac 25 -last_frame 10 -frame_dt 1 -clamp -mu 0 -scale_mass 3 -analyze_u_modes -max_ke 15"

opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
opt_name=("default" "regular" "default-cubic" "regular-cubic")

LO=32
HI=256
SKIP=32

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        for r in `seq $HI -$SKIP $LO` ; do
            dt=`perl -e "print (1.0/$r)"`
            DT="-max_dt $dt -min_dt $dt"
            echo $ARGS $DT -no_affine ${opt[$o]} -resolution $r -dump_modes_freq $r -o $NAME/pic-${opt_name[$o]}-$r
            echo $ARGS $DT -affine ${opt[$o]} -resolution $r -dump_modes_freq $r -o $NAME/apic-${opt_name[$o]}-$r
            #echo $ARGS $DT -no_affine ${opt[$o]} -flip 1 -resolution $r -dump_modes_freq $r -o $NAME/flip-${opt_name[$o]}-$r
        done
    done | xargs -P 16 -n 1 -d '\n' bash -c
fi

for s in pic apic ; do
    for r in `seq $HI -$SKIP $LO` ; do
        for o in `seq 0 $((${#opt[@]}-1))` ; do
            ./vort_parse_data.pl < $NAME/$s-${opt_name[$o]}-$r/common/log.txt > $NAME/$s-${opt_name[$o]}-$r/norms.txt
        done
    done
done

RES=64
for o in ${opt_name[@]} ; do
    sed -e "s/xxx/$o/" -e "s/rrr/$RES/" vort_dambreak_plot.tex  > $NAME/plot-dambreak-vort-$o.tex
    sed -e "s/xxx/$o/" -e "s/rrr/$RES/" vel_dambreak_plot.tex  > $NAME/plot-dambreak-vel-$o.tex
done

for i in {pic,apic}-{default,regular} ; do
    for j in {0..10} ; do
        convert -scale 512 $NAME/$i-$RES/modes-$(($j*$RES))-x.png `printf $NAME/$i-dambreak-vort-%02d.png $j`
        convert -scale 512 $NAME/$i-$RES/modes-$(($j*$RES))-x.png `printf $NAME/$i-dambreak-vort-%02d.eps $j`
    done
done

colors=('black' 'red' 'blue' 'green' 'cyan' 'purple' 'yellow' 'orange')
template='\\addplot [mark=none,solid,color=COLOR] table[x=time,y=VVV]{ALGO-RES\/norms.txt};'
for i in {pic,apic}-{default,regular} ; do
    for v in {vel,vort} ; do
        content=''
        legend_content=''
        nline=0
        for r in `seq $HI -$SKIP $LO` ; do
            line=`echo "$template" | sed -e "s/VVV/$v/g; s/ALGO/$i/g; s/RES/$r/g; s/COLOR/${colors[$nline]}/g"`
            content="$content$line\n"
            legend_content="$legend_content$r,"
            nline=$((nline+1))
        done
        sed -e "s/XXXXXX/$content/g; s/LLLLLL/$legend_content/g" ${v}_dambreak_refinement_plot.tex > $NAME/${i}-${v}-refinement-plot.tex
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