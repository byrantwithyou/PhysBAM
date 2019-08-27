#!/bin/bash

NAME=conv-3d

export OPENBLAS_NUM_THREADS=1
ARGS="../fem -3d -q -threads 12"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
names=(wide grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
pinv=(0 0 0 0 1 0)
LO=2
HI=5
ANA="-u 'u=sin(14*x)*y+cos(15*y)*z+x*y,v=cos(14*x)*cos(16*y)+sin(15*y)*x+x*x+y*z-1,w=sin(17*z)*y+cos(15*x)*z' -p 'p=sin(15*x+14*y+1)+cos(16*z)'"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for i in `seq 0 $((${#tests[@]}-1))` ; do
        c=${tests[$i]}
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/$c-r$r -refine $r $ANA ../$c.txt
        done
        p=${pinv[$i]}
        if [ "X$p" = "X1" ] ; then
            for r in `seq $HI -1 $LO` ; do
                echo $ARGS -pinv -o $NAME/$c-pinv-r$r -refine $r $ANA ../$c-pinv.txt;
            done
        fi
    done | xargs -P 1 -n 1 -d '\n' bash -c > /dev/null
fi


for i in `seq 0 $((${#tests[@]}-1))` ; do
    c=${tests[$i]}
    name=${names[$i]}
    q=${pinv[$i]}

    S=()
    if [ "X$q" = "X1" ] ; then
        S=("-" "-pinv-")
    else
        S=("-")
    fi

    for k in ${S[@]} ; do
        cat <<EOF > $NAME/$c${k}p.txt
res linf l2
EOF
        cat <<EOF > $NAME/$c${k}v.txt
res linf l2
EOF
        for r in `seq $LO 1 $HI` ; do
            res=`bc <<< "2*$r"`
            grep 'l-2' $NAME/$c${k}r$r/common/log.txt |\
                sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/$res \3 \4/g" >> $NAME/$c${k}p.txt
            grep 'l-2' $NAME/$c${k}r$r/common/log.txt |\
                sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/$res \1 \2/g" >> $NAME/$c${k}v.txt
        done
    done

    if [ "X$q" = "X1" ] ; then
        cp conv_plot_pinv.tex $NAME/plot-$c-p.tex
        cp conv_plot_pinv.tex $NAME/plot-$c-v.tex
    else
        cp conv_plot.tex $NAME/plot-$c-p.tex
        cp conv_plot.tex $NAME/plot-$c-v.tex
    fi

    sed -i -e 's/LLLL/p/g' -e "s/XXXX/$c-p/g" -e "s/TTTT/$name pressure/g" $NAME/plot-$c-p.tex
    sed -i -e 's/LLLL/\\mathbf{v}/g' -e "s/XXXX/$c-v/g" -e "s/TTTT/$name velocity/g" $NAME/plot-$c-v.tex

    V=`./conv_regression.pl $NAME/$c-v.txt $NAME/$c-v-regres.txt`
    read ova ovb <<< "$V"
    sed -i -e "s/IIII/$ova/g" -e "s/EEEE/$ovb/g" $NAME/plot-$c-v.tex

    P=`./conv_regression.pl $NAME/$c-p.txt $NAME/$c-p-regres.txt`
    read opa opb <<< "$P"
    sed -i -e "s/IIII/$opa/g" -e "s/EEEE/$opb/g" $NAME/plot-$c-p.tex

    if [ "X$q" = "X1" ] ; then
        sed -i -e "s/YYYY/$c-pinv-p/g" $NAME/plot-$c-p.tex
        sed -i -e "s/YYYY/$c-pinv-v/g" $NAME/plot-$c-v.tex

        V=`./conv_regression.pl $NAME/$c-pinv-v.txt $NAME/$c-pinv-v-regres.txt`
        read ova ovb <<< "$V"
        sed -i -e "s/MMMM/$ova/g" -e "s/FFFF/$ovb/g" $NAME/plot-$c-v.tex

        P=`./conv_regression.pl $NAME/$c-pinv-p.txt $NAME/$c-pinv-p-regres.txt`
        read opa opb <<< "$P"
        sed -i -e "s/MMMM/$opa/g" -e "s/FFFF/$opb/g" $NAME/plot-$c-p.tex
    fi
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
