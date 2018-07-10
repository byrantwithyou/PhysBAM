#!/bin/bash

NAME=eig-regular-seeding

RES=512
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../../fourier-apic/fourier_mac -resolution 16 -size $RES -dump_particles -dump_eigenvalues >&/dev/null"

np=("-ppd 1" "-ppd 2" "-ppd 3" "-ppd 4" "-ppd 5" "-ppd 8")
np_name=("ppc1" "ppc4" "ppc9" "ppc16" "ppc25" "ppc64")

order=("-order 1" "-order 2" "-order 3")
order_name=("linear" "quadratic" "cubic")

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
    sed -e "s/aaa/$a/g; s/xxx/linear/g" eig_regular_seeding_plot.tex  > $NAME/eig-regular-seeding-linear-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/quadratic/g" eig_regular_seeding_plot.tex  > $NAME/eig-regular-seeding-quadratic-$a.tex
    sed -e "s/aaa/$a/g; s/xxx/cubic/g" eig_regular_seeding_plot.tex  > $NAME/eig-regular-seeding-cubic-$a.tex
done

Y=(x,0 0,x x,x)
A=(x y xy)
XM=(.5 .125)
YN=(-.04 .8)
FORGET=('%' '%' '')
for a in 0 1 2 ; do
    sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/;s/XM/${XM[0]}/;s/YN/${YN[0]}/" eig_regular_seeding_plot_combined.tex  > $NAME/eig-regular-seeding-combined-${A[$a]}-0.tex
    sed -e "s/aaa/${A[$a]}/g;s/bbb/${Y[$a]}/;s/XM/${XM[1]}/;s/YN/${YN[1]}/;s/FORGET/${FORGET[$a]}/" eig_regular_seeding_plot_combined_loglog.tex  > $NAME/eig-regular-seeding-combined-${A[$a]}-1.tex
done

cp eig_regular_seeding_plot_legends.tex $NAME/eig-regular-seeding-plot-legends.tex

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
