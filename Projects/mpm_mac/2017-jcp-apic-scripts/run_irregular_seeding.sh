#!/bin/bash

NAME=eig-irregular-seeding

RES=64
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../../fourier-apic/fourier_mac -resolution $RES -size $RES -dump_particles -dump_eigenvalues"

default_seed=0

np=("-irreg 1" "-irreg 4" "-irreg 9" "-irreg 25" "-irreg 64")
np_name=("ppc1" "ppc4" "ppc9" "ppc25" "ppc64")

order=("-order 2" "-order 3")
order_name=("quadratic" "cubic")

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for p in `seq 0 $((${#np[@]}-1))` ; do
        for o in `seq 0 $((${#order[@]}-1))` ; do
            pic_folder=$NAME/pic-${np_name[$p]}-${order_name[$o]}
            apic_folder=$NAME/apic-${np_name[$p]}-${order_name[$o]}
            echo $ARGS ${np[$p]} ${order[$o]} -v $pic_folder -o $pic_folder/eigen.png
            echo $ARGS -affine ${np[$p]} ${order[$o]} -v $apic_folder -o $apic_folder/eigen.png
        done
    done | xargs -P 8 -n 1 -d '\n' bash -c
fi

for a in x y xy ; do
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
