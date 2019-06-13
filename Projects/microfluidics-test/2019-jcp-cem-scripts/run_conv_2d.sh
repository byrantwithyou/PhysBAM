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

names=("" "" "ul-inf" "ul-2" "pl-inf" "pl-2")
for c in ${tests[@]} ; do
    cat <<EOF > $NAME/$c.txt
EOF
    for r in `seq $LO 1 $HI` ; do
        grep 'l-2' $NAME/$c-r$r/common/log.txt |\
            sed "s/.*l-inf \([^ ]*\).*l-2 \([^ ]*\).*l-inf \([^ ]*\).*l-2 \([^<]*\).*/$r \1 \2 \3 \4/g" >> $NAME/$c.txt
    done

    for i in 2 3 4 5 ; do
        gnuplot -e 'set terminal pdf' \
            -e "set output \"$NAME/$c-${names[$i]}.pdf\"" \
            -e "fit a*x+b \"$NAME/$c.txt\" u (log10(\$1)):(log10(\$$i)) via a,b" \
            -e "plot \"$NAME/$c.txt\" u (log10(\$1)):(log10(\$$i)),a*x+b title sprintf(\"%.2f*x+%.2f\",a,b)" \
            2>/dev/null
    done
done
