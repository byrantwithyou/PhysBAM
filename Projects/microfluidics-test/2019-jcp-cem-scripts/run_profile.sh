#!/bin/bash

NAME=profile

export OPENBLAS_NUM_THREADS=1
PWD=$(readlink -f .)
CMD=$PWD/../"fem_profile -q"
BIN=$PWD/../../../build/nocona/profile/Projects/microfluidics-test/fem
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(grid10 grid10 fab-0 fab-0)
config=("-refine 6" "-3d" "-refine 6" "-3d -refine 2")

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for i in `seq 0 $((${#tests[@]}-1))` ; do
        mkdir $NAME/${tests[$i]}-r$i
        echo "(cd $NAME/${tests[$i]}-r$i; \
            $CMD -o . ${config[$i]} $PWD/../${tests[$i]}.txt; \
            gprof $BIN > report.txt;)"
    done | xargs -P 8 -n 1 -d '\n' bash -c > /dev/null
fi

for f in `find $NAME -name report.txt`;
do
    sed -i "s/PhysBAM:://g; \
        s/VECTOR<double, 2>/TV2/g; \
        s/VECTOR<double, 3>/TV3/g" $f;
done
