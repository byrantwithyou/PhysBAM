#!/bin/bash

NAME=vort-dambreak

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../mpm_mac 25 -last_frame 10 -frame_dt 1 -clamp -mu 0 -scale_mass 3 -analyze_energy_vort -extrap p"

opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
opt_name=("default" "regular" "default-cubic" "regular-cubic")

LO=32
HI=192
SKIP=32
RES=96

MOVIE_RES=64
FRAME_LO=500
FRAME_HI=600
FRAME_SKIP=25

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        for r in `seq $HI -$SKIP $LO` ; do
            dt=`perl -e "print (1.0/$r)"`
            DT="-max_dt $dt -min_dt $dt"
            PVORT=""
            if [ "X$r" = "X$MOVIE_RES" ]; then
                PVORT="-particle_vort"
            fi
            echo $ARGS $DT $PVORT -no_affine ${opt[$o]} -resolution $r -dump_modes_freq $r -o $NAME/pic-${opt_name[$o]}-$r
            echo $ARGS $DT $PVORT -affine ${opt[$o]} -resolution $r -dump_modes_freq $r -o $NAME/apic-${opt_name[$o]}-$r
            echo $ARGS $DT $PVORT -no_affine ${opt[$o]} -flip 1 -resolution $r -dump_modes_freq $r -o $NAME/flip-${opt_name[$o]}-$r
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

for s in pic apic flip ; do
    for o in ${opt_name[@]} ; do
        folder=$NAME/movie-$s-$o-$MOVIE_RES
        mkdir -p $folder
        for t in `seq $FRAME_LO $FRAME_SKIP $FRAME_HI` ; do
            plot="\\\\dbplot{..\\/$s-$o-$MOVIE_RES\\/pvort-$t}\\n"
            sed -e "s/XXX/$plot/;" dambreak_movie_plot.tex > $folder/frame-$t.tex
        done
        cat <<EOF > $folder/SConstruct
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
            cd $folder
            scons -j 8
        )
    done
done

for o in ${opt_name[@]} ; do
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Vorticity/; s/DDD/vort/" vort_dambreak_plot.tex  > $NAME/plot-dambreak-vort-$o.tex
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Kinetic energy/; s/DDD/ke/" dambreak_plot.tex  > $NAME/plot-dambreak-ke-$o.tex
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Total energy/; s/DDD/te/" dambreak_plot.tex  > $NAME/plot-dambreak-te-$o.tex
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
