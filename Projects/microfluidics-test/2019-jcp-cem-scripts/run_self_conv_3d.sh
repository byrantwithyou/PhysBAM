#!/bin/bash

NAME=self-conv-3d
SOL_FOLDER=sol-3d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -3d -q -threads 16 -mu 8.9e-4"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
LO=2
HI=5
SOL=8

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            $ARGS -o $NAME/$c-r$r -refine $r -check_sol $SOL_FOLDER/$c-r$SOL/sol.gz ../$c.txt;
        done
    done
fi

for c in ${tests[@]} ; do
    cat <<EOF > $NAME/$c-p.txt
res linf l2
EOF
    cat <<EOF > $NAME/$c-v.txt
res linf l2
EOF
    for r in `seq $LO 1 $HI` ; do
        res=`bc <<< "2*$r"`
        maxsol=`grep 'sol max' $SOL_FOLDER/$c-r$SOL/common/log.txt | sed "s/.*sol maxv \([^ ]*\) maxp \([^<]*\).*/\1 \2/g"`
        read maxv maxp <<< "$maxsol"
        err=`grep 'l-2' $NAME/$c-r$r/common/log.txt | sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/\1 \2 \3 \4/g"`;
        read vlinf vl2 plinf pl2 <<< "$err"
        echo "$res `perl -e \"{print $plinf/$maxp}\"` `perl -e \"{print $pl2/$maxp}\"`" >> $NAME/$c-p.txt
        echo "$res `perl -e \"{print $vlinf/$maxv}\"` `perl -e \"{print $vl2/$maxv}\"`" >> $NAME/$c-v.txt
    done

    sed -e 's/LLLL/p/g' -e "s/XXXX/$c-p/g" -e "s/TTTT/$c pressure/g" conv_plot.tex  > $NAME/plot-$c-p.tex
    sed -e 's/LLLL/\\mathbf{v}/g' -e "s/XXXX/$c-v/g" -e "s/TTTT/$c velocity/g" conv_plot.tex  > $NAME/plot-$c-v.tex

    V=`./conv_regression.pl $NAME/$c-v.txt $NAME/$c-v-regres.txt`
    read ova ovb <<< "$V"
    sed -i -e "s/IIII/$ova/g" -e "s/EEEE/$ovb/g" $NAME/plot-$c-v.tex

    P=`./conv_regression.pl $NAME/$c-p.txt $NAME/$c-p-regres.txt`
    read opa opb <<< "$P"
    sed -i -e "s/IIII/$opa/g" -e "s/EEEE/$opb/g" $NAME/plot-$c-p.tex
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
