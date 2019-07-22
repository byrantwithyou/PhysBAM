#!/bin/bash

NAME=scale-refine

export OPENBLAS_NUM_THREADS=1
ARGS="../fem_profile -q -timing -mu 8.9e-4 -threads 16"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(grid20 rgrid0 voronoi-s4)
R2_a=9
R2_b=18
R3_a=2
R3_b=5

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $R2_b -1 $R2_a`; do
            $ARGS -refine $r ../$c.txt -o $NAME/$c-r$r > /dev/null;
            mv timing-* $NAME/$c-r$r;
        done

        for r in `seq $R3_b -1 $R3_a` ; do
            $ARGS -refine $r ../$c.txt -3d -o $NAME/$c-3d-r$r > /dev/null;
            mv timing-* $NAME/$c-3d-r$r;
        done
    done
fi

cat <<EOF > $NAME/stats.csv
Name,Resolution,Total dofs,Blocks,Tasks,Avg dofs/block
EOF
for c in ${tests[@]} ; do
    for r in $R2_a $R2_b ; do
        tot=`grep "dofs 2d" $NAME/$c-r$r/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
        blks=`grep ">blocks:" $NAME/$c-r$r/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
        jbs=`grep "jobs:" $NAME/$c-r$r/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
        dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
        rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
        rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
        res=`perl -e "{print 2*$r;}"`
        echo "$c,$res,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
    done
done
for c in ${tests[@]} ; do
    for r in $R3_a $R3_b ; do
        tot=`grep "dofs 3d" $NAME/$c-3d-r$r/common/log.txt | sed 's/.*total: \([^<]*\).*/\1/g'`
        blks=`grep ">blocks:" $NAME/$c-3d-r$r/common/log.txt | sed 's/.*blocks: \([^<]*\).*/\1/g'`
        jbs=`grep "jobs:" $NAME/$c-3d-r$r/common/log.txt | sed 's/.*jobs: \([^<]*\).*/\1/g'`
        dofspb=`perl -e "{printf('%.1f',$tot/$blks);}"`
        rndtot=`perl -e "{printf('%.1f M',$tot/1e6);}"`
        rndjbs=`perl -e "{printf('%.1f K',$jbs/1e3);}"`
        res=`perl -e "{print 2*$r;}"`
        echo "$c-3d,$res,$rndtot,$blks,$rndjbs,$dofspb" >> $NAME/stats.csv
    done
done

for c in ${tests[@]} ; do
    cat <<EOF > $NAME/$c.txt
res prep solving
EOF
    for r in `seq $R2_a $R2_b`; do
        res=`perl -e "{print 2*$r;}"`
        echo "$res `./timing_parse.pl < $NAME/$c-r$r/common/log.txt`" >> $NAME/$c.txt
    done

    cat <<EOF > $NAME/$c-3d.txt
res prep solving
EOF
    for r in `seq $R3_a $R3_b`; do
        res=`perl -e "{print 2*$r;}"`
        echo "$res `./timing_parse.pl < $NAME/$c-3d-r$r/common/log.txt`" >> $NAME/$c-3d.txt
    done
done

sed -e "s/DDDD/18:36/g" timing_vs_res_plot.tex > $NAME/plot.tex
sed -e "s/DDDD/4:10/g" timing_vs_res_plot.tex > $NAME/plot-3d.tex
for i in `seq 0 $((${#tests[@]}-1))` ; do
    for dim in "" "-3d" ; do
        c=${tests[$i]}$dim
        rm fit.log
        gnuplot -e "fit a*x+b \"$NAME/$c.txt\" u (log10(\$1)):(log10(\$3)) via a,b" 2>/dev/null
        a=`grep Final -A 4 fit.log| grep "^a" | sed 's/a *= \([^ ]*\).*/\1/g'`
        b=`grep Final -A 4 fit.log| grep "^b" | sed 's/b *= \([^ ]*\).*/\1/g'`
        order=`perl -e "{printf('%.2f',$a);}"`
        b10=`perl -e "{print 10**$b}"`
        sed -i -e "s/XXXX$i/$c/g" -e "s/EEEE$i/$a/g" -e "s/CCCC$i/$b10/g" -e "s/OOOO$i/$order/g" $NAME/plot$dim.tex
    done
done

for c in ${tests[@]} ; do
    cat <<EOF > $NAME/waiting-$c.txt
res waiting solving percentage
EOF
    for r in `seq $R2_a $R2_b`; do
        res=`perl -e "{print 2*$r;}"`
        echo "$res `cat $NAME/$c-r$r/timing-*.txt | ./waiting_parse.pl`" >> $NAME/waiting-$c.txt
    done

    cat <<EOF > $NAME/waiting-$c-3d.txt
res waiting solving percentage
EOF
    for r in `seq $R3_a $R3_b`; do
        res=`perl -e "{print 2*$r;}"`
        echo "$res `cat $NAME/$c-3d-r$r/timing-*.txt | ./waiting_parse.pl`" >> $NAME/waiting-$c-3d.txt
    done
done

sed -e 's/XXXX/waiting-grid20/g' -e 's/YYYY/waiting-rgrid0/g' -e 's/ZZZZ/waiting-voronoi-s4/g' \
    waiting_percentage_vs_res_plot.tex  > $NAME/plot-waiting.tex 
sed -e 's/XXXX/waiting-grid20-3d/g' -e 's/YYYY/waiting-rgrid0-3d/g' -e 's/ZZZZ/waiting-voronoi-s4-3d/g' \
    waiting_percentage_vs_res_plot.tex  > $NAME/plot-waiting-3d.tex 

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
