#!/bin/bash

NAME=conv-2d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -q -threads 6"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
LO=2
HI=16
ANA="-u 'u=sin(x)*y+cos(y)+x*y,v=cos(x)*cos(y)+sin(y)*x+x*x-1' -p 'p=sin(x+y+1)'"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/$c-r$r -refine $r $ANA ../$c.txt;
        done
    done | xargs -P 4 -n 1 -d '\n' bash -c > /dev/null
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
        grep 'l-2' $NAME/$c-r$r/common/log.txt |\
            sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/$res \3 \4/g" >> $NAME/$c-p.txt
        grep 'l-2' $NAME/$c-r$r/common/log.txt |\
            sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/$res \1 \2/g" >> $NAME/$c-v.txt
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
