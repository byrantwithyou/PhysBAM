#!/bin/bash

NAME=eig-irregular-seeding

RES=512
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

SEED=42

ARGS="../../fourier-apic/fourier_mac -resolution 16 -size $RES -dump_particles -dump_eigenvalues -seed $SEED >&/dev/null"


np=("-irreg 1" "-irreg 4" "-irreg 9" "-irreg 16" "-irreg 25" "-irreg 64")
np_name=("ppc1" "ppc4" "ppc9" "ppc16" "ppc25" "ppc64")

order=("-order 1" "-order 2" "-order 3")
order_name=("linear" "quadratic" "cubic")

px=(0.00001 0.499999 0.25 0.44140625 0.22265625 0.34765625)
py=(0.50390625 0.98828125 0.75390625 0.76171875 0.58203125 0.546875)

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for i in `seq 0 $((${#px[@]}-1))` ; do
        x=${px[$i]}
        y=${py[$i]}
        for o in `seq 0 $((${#order[@]}-1))` ; do
            coord=$x-$y
            pic_folder=$NAME/pic-single-$coord-${order_name[$o]}
            apic_folder=$NAME/apic-single-$coord-${order_name[$o]}
            echo $ARGS -px $x -py $y ${order[$o]} -v $pic_folder -o $pic_folder/eigen.png
            echo $ARGS -affine -px $x -py $y ${order[$o]} -v $apic_folder -o $apic_folder/eigen.png
        done
    done | xargs -P 8 -n 1 -d '\n' bash -c

    for p in `seq 0 $((${#np[@]}-1))` ; do
        for o in `seq 0 $((${#order[@]}-1))` ; do
            pic_folder=$NAME/pic-${np_name[$p]}-${order_name[$o]}
            apic_folder=$NAME/apic-${np_name[$p]}-${order_name[$o]}
            echo $ARGS ${np[$p]} ${order[$o]} -v $pic_folder -o $pic_folder/eigen.png
            echo $ARGS -affine ${np[$p]} ${order[$o]} -v $apic_folder -o $apic_folder/eigen.png
        done
    done | xargs -P 8 -n 1 -d '\n' bash -c
fi

colors=('black' 'red' 'blue' 'green' 'cyan' 'purple')
template='\\addplot [mark=none,solid,color=COLOR] table[x=x,y=eig]{apic-single-X-Y-xxx\/eig-aaa.txt};'
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
sed -e "s/XXXXXX/$content/g" eig_irregular_seeding_ppc1_plot.tex > $NAME/eig_irregular_seeding_apic_ppc1.template
sed -e "s/apic/pic/g" $NAME/eig_irregular_seeding_apic_ppc1.template > $NAME/eig_irregular_seeding_pic_ppc1.template

Y=(x,0 0,x x,x)
A=(x y xy)

for a in 0 1 2 ; do
    for o in ${order_name[@]} ; do
        sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" $NAME/eig_irregular_seeding_apic_ppc1.template > $NAME/eig-irregular-seeding-apic-ppc1-$o-${A[$a]}.tex
        sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" $NAME/eig_irregular_seeding_pic_ppc1.template > $NAME/eig-irregular-seeding-pic-ppc1-$o-${A[$a]}.tex
        sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/; s/xxx/$o/g" eig_irregular_seeding_plot.tex  > $NAME/eig-irregular-seeding-$o-${A[$a]}.tex
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
