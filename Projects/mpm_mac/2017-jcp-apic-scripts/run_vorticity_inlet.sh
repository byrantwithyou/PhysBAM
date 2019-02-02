#!/bin/bash

NAME=vort-inlet

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

INLETV=0.2
STOP=80
ARGS="../mpm_mac 28 -frame_dt 1 -mu 0 -scale_mass 3 -T $INLETV -T $STOP -reseed -analyze_energy_vort -clamp -solver_tolerance 1e-2 -warm_start"

#opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
#opt_name=("default" "regular" "default-cubic" "regular-cubic")
opt=("")
opt_name=("default")

LO=32
HI=64
SKIP=32
RES=64
T=200
MOVIE_T=200

MOVIE_RES=64
FRAME_SKIP=150

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for o in `seq 0 $((${#opt[@]}-1))` ; do
        for r in `seq $HI -$SKIP $LO` ; do
            dt=`perl -e "print (1.0/$r)"`
            DT="-max_dt $dt -min_dt $dt"
            PVORT=""
            last_frame=$T
            if [ "X$r" = "X$MOVIE_RES" ]; then
                last_frame=$MOVIE_T
                PVORT="-particle_vort"
            fi
            LF="-last_frame $last_frame"
            echo $ARGS $DT $LF $PVORT -no_affine ${opt[$o]} -resolution $r -o $NAME/pic-${opt_name[$o]}-$r
            echo $ARGS $DT $LF $PVORT -affine ${opt[$o]} -resolution $r -o $NAME/apic-${opt_name[$o]}-$r
            echo $ARGS $DT $LF $PVORT -no_affine ${opt[$o]} -flip 1 -resolution $r -o $NAME/flip-${opt_name[$o]}-$r
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

DS=../dump_strlns
cp plot.gplt $NAME
cp plot_open.gplt $NAME
for a in pic apic flip ; do
    for f in `seq 80 10 200` ; do
        $DS $NAME/$a-default-$RES/ -frame $f -o $NAME/dump-$a-default-$RES-$f -arrow_sep 18;
        grep COORD $NAME/dump-$a-default-$RES-$f/common/log.txt | sed 's/.*COORD \([^<]*\).*/\1/g' > $NAME/coord-$a-default-$RES-$f;
        grep ARROW $NAME/dump-$a-default-$RES-$f/common/log.txt | sed 's/.*ARROW \([^<]*\).*/\1/g' > $NAME/arrow-$a-default-$RES-$f;
        (
            cd $NAME
            if [ $f = 80 ] ; then
                gnuplot -e "fn='$a-default-$RES-$f'" plot_open.gplt
            else
                gnuplot -e "fn='$a-default-$RES-$f'" plot.gplt
            fi
        )
    done
done

for o in ${opt_name[@]} ; do
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Vorticity/; s/DDD/vort/; s/XMAX/$T/; s/YMAX/2.5/;" vort_dambreak_plot.tex  > $NAME/plot-inlet-vort-$o.tex
    sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Kinetic energy/; s/DDD/ke/; s/XMAX/$T/; s/YMAX/0.025/;" dambreak_plot.tex  > $NAME/plot-inlet-ke-$o.tex
    #sed -e "s/xxx/$o/; s/rrr/$RES/; s/LLL/Total energy/; s/DDD/te/; s/XMAX/$T/; s/YMAX/0.02/;" te_dambreak_plot.tex  > $NAME/plot-inlet-te-$o.tex
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
