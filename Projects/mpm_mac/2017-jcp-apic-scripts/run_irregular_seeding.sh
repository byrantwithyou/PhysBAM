#!/bin/bash

NAME=eig-irregular-seeding

RES=64
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

SEED=42

ARGS="../../fourier-apic/fourier_mac -resolution $RES -size $RES -dump_particles -dump_eigenvalues -seed $SEED"


np=("-irreg 1" "-irreg 4" "-irreg 9" "-irreg 16" "-irreg 25" "-irreg 64")
np_name=("ppc1" "ppc4" "ppc9" "ppc16" "ppc25" "ppc64")

order=("-order 2" "-order 3")
order_name=("quadratic" "cubic")

px=($(seq 0 0.1 0.5))
#py=($(seq 0 0.1 0.5))
py=($(seq 0.1 0.1 0.1))

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for x in ${px[@]}; do
        for y in ${py[@]}; do
            for o in `seq 0 $((${#order[@]}-1))` ; do
                coord=$x-$y
                pic_folder=$NAME/pic-single-$coord-${order_name[$o]}
                apic_folder=$NAME/apic-single-$coord-${order_name[$o]}
                echo $ARGS -px $x -py $y ${order[$o]} -v $pic_folder -o $pic_folder/eigen.png
                echo $ARGS -affine -px $x -py $y ${order[$o]} -v $apic_folder -o $apic_folder/eigen.png
            done
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

colors=('black' 'red' 'blue' 'green' 'purple' 'cyan')
template='\\addplot [mark=none,solid,color=COLOR] table[x=waven,y=eig]{apic-single-X-Y-xxx\/eig-aaa.txt};'
content=''
legend=''
nline=0
for x in ${px[@]}; do
    for y in ${py[@]}; do
        line=`echo $template | sed -e "s/X/$x/g; s/Y/$y/g; s/COLOR/${colors[$nline]}/g"`
        nline=$((nline+1))
        legend="$legend($x,$y)\\\\\\\\"
        content=$content$line'\n'
    done
done
legend_entries="legend entries={$legend}"

sed -e "s/XXXXXX/$content/g; s/LEGEND/$legend_entries/g" eig_irregular_seeding_ppc1_plot.tex > $NAME/eig_irregular_seeding_apic_ppc1.template
sed -e "s/apic/pic/g" $NAME/eig_irregular_seeding_apic_ppc1.template > $NAME/eig_irregular_seeding_pic_ppc1.template

for a in x y xy ; do
    sed -e "s/aaa/$a/g; s/xxx/quadratic/g" $NAME/eig_irregular_seeding_apic_ppc1.template > $NAME/eig-irregular-seeding-apic-ppc1-quadratic-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/cubic/g" $NAME/eig_irregular_seeding_apic_ppc1.template > $NAME/eig-irregular-seeding-apic-ppc1-cubic-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/quadratic/g" $NAME/eig_irregular_seeding_pic_ppc1.template > $NAME/eig-irregular-seeding-pic-ppc1-quadratic-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/cubic/g" $NAME/eig_irregular_seeding_pic_ppc1.template > $NAME/eig-irregular-seeding-pic-ppc1-cubic-$a.tex

    sed -e "s/aaa/$a/g; s/xxx/quadratic/g" eig_irregular_seeding_plot.tex  > $NAME/eig-irregular-seeding-quadratic-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/cubic/g" eig_irregular_seeding_plot.tex  > $NAME/eig-irregular-seeding-cubic-$a.tex
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
