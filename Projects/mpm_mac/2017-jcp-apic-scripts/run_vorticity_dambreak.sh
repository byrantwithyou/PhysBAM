#!/bin/bash

NAME=vort-dambreak

FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../mpm_mac 25 -frame_dt 1 -clamp -mu 0 -scale_mass 3 -analyze_energy_vort -extrap p"

opt=("" "-regular_seeding" "-order 3" "-regular_seeding -order 3")
opt_name=("default" "regular" "default-cubic" "regular-cubic")

LO=32
HI=96
SKIP=32
RES=96
T=10
MOVIE_T=15

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

for s in pic apic flip ; do
    for o in ${opt_name[@]} ; do
        folder=$NAME/movie-$s-$o-$MOVIE_RES
        mkdir -p $folder
        last=`perl -e "print (4*$FRAME_SKIP)"`
        for t in `seq 0 $FRAME_SKIP $last` ; do
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

for s in apic pic flip ; do
    o="default"
    folder=$NAME/fullmovie-$s-$o-$MOVIE_RES
    mkdir -p $folder
    last=`perl -e "print ($MOVIE_T*$MOVIE_RES)"`
    for t in `seq 0 $last` ; do
        plot="\\\\dbplot{..\\/$s-$o-$MOVIE_RES\\/pvort-$t}\\n"
        sed -e "s/XXX/$plot/;" dambreak_movie_plot.tex > $folder/frame-$(printf "%03d" $t).tex
    done
    cat <<EOF > $folder/SConstruct
import os
import re
env=Environment(ENV = os.environ)
r=re.compile(".*\.tex$")
for f in [x for x in os.listdir(".") if r.match(x)]:
    t=env.DVI(f)
    env.PDF(t)
EOF
        (
            cd $folder
            scons -j 16
        )

    (
        cd $folder
        for t in `seq 0 $last` ; do
            tname=$(printf "%03d" $t)
            echo convert -density 300 -quality 100 -flatten -sharpen 0x1.0 -scale 1024x768 frame-$tname.pdf frame-$tname.png
        done | xargs -P 16 -n 1 -d '\n' bash -c
        ffmpeg -y -i "frame-%03d.png" -r 30 -c:v libx264 -crf 20 -pix_fmt yuv420p movie.mov
    )
done

cp dambreak_movie_colorbar.tex $NAME/movie-colorbar.tex

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
