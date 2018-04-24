#!/bin/bash

NAME=eig-p2p-regular-seeding

RES=128
LOWRES=32
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

ARGS="../../fourier-apic/fourier_mac_p -pressure -dump_particles >&/dev/null"

np=("-ppd 4")
np_name=("ppc16")

order=("-order 1" "-order 2" "-order 3")
order_name=("linear" "quadratic" "cubic")

xpic_order=(1 2 5 10 15 20 25)

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for r in $RES $LOWRES ; do
        for p in `seq 0 $((${#np[@]}-1))` ; do
            for o in `seq 0 $((${#order[@]}-1))` ; do
                flip_folder=$NAME/flip-${np_name[$p]}-${order_name[$o]}-$r
                echo $ARGS -flip 1 ${np[$p]} ${order[$o]} -v $flip_folder -resolution $r -size $r
                for x in ${xpic_order[@]} ; do
                    xpic_folder=$NAME/xpic$x-${np_name[$p]}-${order_name[$o]}-$r
                    echo $ARGS -xpic $x ${np[$p]} ${order[$o]} -v $xpic_folder -resolution $r -size $r
                done
            done
        done
    done | xargs -P 16 -n 1 -d '\n' bash -c
fi

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
