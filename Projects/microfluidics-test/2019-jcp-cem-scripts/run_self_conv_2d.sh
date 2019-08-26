#!/bin/bash

NAME=self-conv-2d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -q -threads 16 -mu 8.9e-4"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
names=(wide grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
LO=2
HI=16
SOL=32

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        $ARGS -o $NAME/$c-r$SOL -refine $SOL -dump_sol ../$c.txt > /dev/null
    done
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/$c-r$r -refine $r -check_sol $NAME/$c-r$SOL/sol.gz ../$c.txt;
        done
    done | xargs -P 4 -n 1 -d '\n' bash -c > /dev/null

    $ARGS -o $NAME/simple-error -dump_sol -refine $HI ../simple.txt;
    ../sol_viewer $NAME/simple-error/sol.gz -res 1024 -ref $NAME/simple-r$SOL/sol.gz -o $NAME/simple-plot-error
fi

for i in `seq 0 $((${#tests[@]}-1))` ; do
    c=${tests[$i]}
    name=${names[$i]}
    cat <<EOF > $NAME/$c-p.txt
res linf l2
EOF
    cat <<EOF > $NAME/$c-v.txt
res linf l2
EOF
    for r in `seq $LO 1 $HI` ; do
        res=`bc <<< "2*$r"`
        maxsol=`grep 'sol max' $NAME/$c-r$SOL/common/log.txt | sed "s/.*sol maxv \([^ ]*\) maxp \([^<]*\).*/\1 \2/g"`
        read maxv maxp <<< "$maxsol"
        err=`grep 'l-2' $NAME/$c-r$r/common/log.txt | sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/\1 \2 \3 \4/g"`;
        read vlinf vl2 plinf pl2 <<< "$err"
        echo "$res `perl -e \"{print $plinf/$maxp}\"` `perl -e \"{print $pl2/$maxp}\"`" >> $NAME/$c-p.txt
        echo "$res `perl -e \"{print $vlinf/$maxv}\"` `perl -e \"{print $vl2/$maxv}\"`" >> $NAME/$c-v.txt
    done

    sed -e 's/LLLL/p/g' -e "s/XXXX/$c-p/g" -e "s/TTTT/$name pressure/g" conv_plot.tex  > $NAME/plot-$c-p.tex
    sed -e 's/LLLL/\\mathbf{v}/g' -e "s/XXXX/$c-v/g" -e "s/TTTT/$name velocity/g" conv_plot.tex  > $NAME/plot-$c-v.tex

    V=`./conv_regression.pl $NAME/$c-v.txt $NAME/$c-v-regres.txt`
    read ova ovb <<< "$V"
    sed -i -e "s/IIII/$ova/g" -e "s/EEEE/$ovb/g" $NAME/plot-$c-v.tex

    P=`./conv_regression.pl $NAME/$c-p.txt $NAME/$c-p-regres.txt`
    read opa opb <<< "$P"
    sed -i -e "s/IIII/$opa/g" -e "s/EEEE/$opb/g" $NAME/plot-$c-p.tex
done

(
    cd $NAME/simple-plot-error;
    X0=750;
    Y0=80;
    X1=790;
    Y1=120;
    S0=40;
    S1=40;

    convert -crop ${S0}x${S1}+$X0+$Y0 dp.png dp-large.png
    convert -crop ${S0}x${S1}+$X0+$Y0 dv.png dv-large.png
    convert -crop ${S0}x${S1}+$X0+$Y0 err-p.png err-p-large.png
    convert -crop ${S0}x${S1}+$X0+$Y0 err-v.png err-v-large.png

    convert -fill none -stroke red -strokewidth 2 -draw "rectangle $X0,$Y0,$X1,$Y1" dp.png dp-anno.png
    convert -fill none -stroke red -strokewidth 2 -draw "rectangle $X0,$Y0,$X1,$Y1" dv.png dv-anno.png
    convert -fill none -stroke red -strokewidth 2 -draw "rectangle $X0,$Y0,$X1,$Y1" err-p.png err-p-anno.png
    convert -fill none -stroke red -strokewidth 2 -draw "rectangle $X0,$Y0,$X1,$Y1" err-v.png err-v-anno.png
)

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
