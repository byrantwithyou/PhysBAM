#!/bin/bash

NAME=eig-p2p-irregular-seeding

RES=128
LOWRES=32
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

SEED=42

ARGS="../../fourier-apic/fourier_mac_p -pressure -dump_particles -dump_eigenvalues -seed $SEED >&/dev/null"

#np=("-irreg 1" "-irreg 4" "-irreg 9" "-irreg 16" "-irreg 25" "-irreg 64")
#np_name=("ppc1" "ppc4" "ppc9" "ppc16" "ppc25" "ppc64")
np=("-irreg 1" "-irreg 4" "-irreg 9")
np_name=("ppc1" "ppc4" "ppc9")

order=("-order 1" "-order 2" "-order 3")
order_name=("linear" "quadratic" "cubic")

xpic_order=(1 2 3 4 5 10 15 20 25)

px=(0.00001 0.499999 0.25 0.44140625 0.22265625 0.34765625)
py=(0.50390625 0.98828125 0.75390625 0.76171875 0.58203125 0.546875)

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for r in $RES $LOWRES ; do
        res=""
        if [ "$r" -eq "$LOWRES" ] ; then
            res="-$r"
        fi
        for i in `seq 0 $((${#px[@]}-1))` ; do
            x=${px[$i]}
            y=${py[$i]}
            for o in `seq 0 $((${#order[@]}-1))` ; do
                coord=$x-$y
                flip_folder=$NAME/flip-single-$coord-${order_name[$o]}$res
                echo $ARGS -flip 1 -px $x -py $y ${order[$o]} -v $flip_folder -resolution $r -size $r
                for xo in ${xpic_order[@]} ; do
                    xpic_folder=$NAME/xpic$xo-single-$coord-${order_name[$o]}$res
                    echo $ARGS -xpic $xo -px $x -py $y ${order[$o]} -v $xpic_folder -resolution $r -size $r
                done
            done
        done | xargs -P 16 -n 1 -d '\n' bash -c

        for p in `seq 0 $((${#np[@]}-1))` ; do
            for o in `seq 0 $((${#order[@]}-1))` ; do
                flip_folder=$NAME/flip-${np_name[$p]}-${order_name[$o]}$res
                echo $ARGS -flip 1 ${np[$p]} ${order[$o]} -v $flip_folder -resolution $r -size $r
                for xo in ${xpic_order[@]} ; do
                    xpic_folder=$NAME/xpic$xo-${np_name[$p]}-${order_name[$o]}$res
                    echo $ARGS -xpic $xo ${np[$p]} ${order[$o]} -v $xpic_folder -resolution $r -size $r
                done
            done
        done | xargs -P 16 -n 1 -d '\n' bash -c
    done
fi

colors=('black' 'red' 'blue' 'green' 'cyan' 'purple')
template='\\addplot [mark=none,solid,color=COLOR] table[x=x,y=eig]{flip-single-X-Y-xxx\/all-nnn-eig-aaa.txt};'
location_template='\\addplot [mark=*, mark options={COLOR}] coordinates{(X,Y)};'
content=''
location_content=''
nline=0
for i in `seq 0 $((${#px[@]}-1))` ; do
    x=${px[$i]}
    shift_x=`echo $x+0.5 | bc`
    y=${py[$i]}
    line=`echo "$template" | sed -e "s/X/$x/g; s/Y/$y/g; s/COLOR/${colors[$nline]}/g"`
    location_line="`echo $location_template | sed -e "s/X/$shift_x/g; s/Y/$y/g; s/COLOR/${colors[$nline]}/g"`"
    nline=$((nline+1))
    content="$content$line\n"
    location_content="$location_content$location_line\n"
done

sed -e "s/XXXXXX/$location_content/g" eig_irregular_seeding_ppc1_particle_plot.tex > $NAME/eig-irregular-seeding-ppc1-particle.tex
sed -e "s/XXXXXX/$content/g; s/nnn/0/g" eig_irregular_seeding_ppc1_plot.tex > $NAME/eig_irregular_seeding_flip_ppc1.template
for xo in ${xpic_order[@]} ; do
    sed -e "s/XXXXXX/$content/g; s/nnn/0/g; s/flip/xpic$xo/g;" eig_irregular_seeding_ppc1_plot.tex > $NAME/eig_irregular_seeding_xpic${xo}_eig0_ppc1.template
    sed -e "s/XXXXXX/$content/g; s/nnn/1/g; s/flip/xpic$xo/g;" eig_irregular_seeding_ppc1_plot.tex > $NAME/eig_irregular_seeding_xpic${xo}_eig1_ppc1.template
done

Y=(x,0 0,x x,x)
A=(x y xy)

for a in 0 1 2 ; do
    for o in ${order_name[@]} ; do
        sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" $NAME/eig_irregular_seeding_flip_ppc1.template > $NAME/eig-irregular-seeding-flip-ppc1-$o-${A[$a]}.tex
        for xo in ${xpic_order[@]} ; do
            sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" $NAME/eig_irregular_seeding_xpic${xo}_eig0_ppc1.template > $NAME/eig-irregular-seeding-xpic$xo-eig0-ppc1-$o-${A[$a]}.tex
            sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" $NAME/eig_irregular_seeding_xpic${xo}_eig1_ppc1.template > $NAME/eig-irregular-seeding-xpic$xo-eig1-ppc1-$o-${A[$a]}.tex
        done
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
