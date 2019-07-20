#!/bin/bash

NAME=scale-cores

export OPENBLAS_NUM_THREADS=1
ARGS="../fem_profile -q -timing -mu 8.9e-4"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(grid20 rgrid0 voronoi-s4)
R2=8;
R3=3;

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for i in `seq 0 4` ; do
            th=`perl -e "{print 2**$i}"`;
            $ARGS -refine $R2 -threads $th ../$c.txt -o $NAME/$c-t$th > /dev/null;
            mv timing-* $NAME/$c-t$th;
            $ARGS -refine $R3 -threads $th ../$c.txt -3d -o $NAME/$c-3d-t$th > /dev/null;
            mv timing-* $NAME/$c-3d-t$th;
        done
    done
fi

cat <<EOF > $NAME/stats.csv
Name,Total dofs,Blocks,Tasks,Avg dofs/block
EOF
for c in ${tests[@]} ; do
    tot=`grep "dofs 2d" $NAME/$c-t1/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
    blks=`grep ">blocks:" $NAME/$c-t1/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
    jbs=`grep "jobs:" $NAME/$c-t1/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
    dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
    rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
    rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
    echo "$c,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
done
for c in ${tests[@]} ; do
    tot=`grep "dofs 3d" $NAME/$c-3d-t1/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
    blks=`grep ">blocks:" $NAME/$c-3d-t1/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
    jbs=`grep "jobs:" $NAME/$c-3d-t1/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
    dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
    rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
    rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
    echo "$c-3d,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
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
    scons
)
