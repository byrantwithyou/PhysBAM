#!/bin/bash

NAME=vort-inlet3

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

INLETV=0.2
STOP=80
T=200
ARGS="../mpm_mac -3d 28 -last_frame $T -frame_dt 1 -mu 0 -scale_mass 3 -T $INLETV -T $STOP -reseed -analyze_energy_vort -clamp -only_log -solver_tolerance 1e-2 -warm_start"

#opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
#opt_name=("default" "regular" "default-cubic" "regular-cubic")
opt=("")
opt_name=("default")

LO=64
HI=64
SKIP=32
RES=64

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        for r in `seq $HI -$SKIP $LO` ; do
            dt=`perl -e "print (1.0/$r)"`
            DT="-max_dt $dt -min_dt $dt"
            echo $ARGS $DT -no_affine ${opt[$o]} -resolution $r -o $NAME/pic-${opt_name[$o]}-$r
            echo $ARGS $DT -affine ${opt[$o]} -resolution $r -o $NAME/apic-${opt_name[$o]}-$r
            echo $ARGS $DT -no_affine ${opt[$o]} -flip 1 -resolution $r -o $NAME/flip-${opt_name[$o]}-$r
        done
    done | xargs -P 16 -n 1 -d '\n' bash -c
fi

for s in pic apic flip ; do
    for r in `seq $HI -$SKIP $LO` ; do
        for o in `seq 0 $((${#opt[@]}-1))` ; do
            ./dambreak_parse_data.pl < $NAME/$s-${opt_name[$o]}-$r/common/log.txt > $NAME/$s-${opt_name[$o]}-$r/norms.txt
        done
    done
done

for o in ${opt_name[@]} ; do
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Vorticity/; s/DDD/vort/; s/XMAX/$T/; s/YMAX/9.0/;" vort_dambreak3_plot.tex  > $NAME/plot-inlet-vort-$o.tex
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Kinetic energy/; s/DDD/ke/; s/XMAX/$T/; s/YMAX/0.03/;" dambreak_plot.tex  > $NAME/plot-inlet-ke-$o.tex
    #sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Total energy/; s/DDD/te/; s/XMAX/10/; s/YMAX/0.3/;" te_dambreak_plot.tex  > $NAME/plot-inlet-te-$o.tex
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
