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

for c in ${tests[@]} ; do
    cat <<EOF > $NAME/$c.txt
threads prep solving
EOF
    for i in `seq 0 4` ; do
        th=`perl -e "{print 2**$i}"`;
        echo "$th `./timing_parse.pl < $NAME/$c-t$th/common/log.txt`" >> $NAME/$c.txt
    done

    cat <<EOF > $NAME/$c-3d.txt
threads prep solving
EOF
    for i in `seq 0 4` ; do
        th=`perl -e "{print 2**$i}"`;
        echo "$th `./timing_parse.pl < $NAME/$c-3d-t$th/common/log.txt`" >> $NAME/$c-3d.txt
    done
done

sed -e 's/XXXX/grid20/g' -e 's/YYYY/rgrid0/g' -e 's/ZZZZ/voronoi-s4/g' \
    -e 's/EEEE/-1/g' -e 's/CCCC/10^4/g' timing_plot.tex  > $NAME/plot.tex 
sed -e 's/XXXX/grid20-3d/g' -e 's/YYYY/rgrid0-3d/g' -e 's/ZZZZ/voronoi-s4-3d/g' \
    -e 's/EEEE/-1/g' -e 's/CCCC/10^5/g' timing_plot.tex  > $NAME/plot-3d.tex 

for c in ${tests[@]} ; do
    cat <<EOF > $NAME/waiting-$c.txt
threads waiting solving percentage
EOF
    for i in `seq 0 4` ; do
        th=`perl -e "{print 2**$i}"`;
        echo "$th `cat $NAME/$c-t$th/timing-*.txt | ./waiting_parse.pl`" >> $NAME/waiting-$c.txt
    done

    cat <<EOF > $NAME/waiting-$c-3d.txt
threads waiting solving percentage
EOF
    for i in `seq 0 4` ; do
        th=`perl -e "{print 2**$i}"`;
        echo "$th `cat $NAME/$c-3d-t$th/timing-*.txt | ./waiting_parse.pl`" >> $NAME/waiting-$c-3d.txt
    done
done

sed -e 's/XXXX/waiting-grid20/g' -e 's/YYYY/waiting-rgrid0/g' -e 's/ZZZZ/waiting-voronoi-s4/g' \
    waiting_percentage_plot.tex  > $NAME/plot-waiting.tex 
sed -e 's/XXXX/waiting-grid20-3d/g' -e 's/YYYY/waiting-rgrid0-3d/g' -e 's/ZZZZ/waiting-voronoi-s4-3d/g' \
    waiting_percentage_plot.tex  > $NAME/plot-waiting-3d.tex 

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
