#!/bin/bash

NAME=conv-2d

export OPENBLAS_NUM_THREADS=6
ARGS="../fem -q"
FULL=1 # Set to 1 for a full rebuild; 0 to skip rerunning the simulations

tests=(0 1 2 3 13 14)
LO=2
HI=24
ANA="-u 'u=sin(x)*y+cos(y)+x*y,v=cos(x*y)+sin(y)*x+x*x-1' -p 'p=sin(x+y+1)'"

if [ "X$FULL" = "X1" ] ; then
    rm -rf $NAME
    mkdir -p $NAME
    for c in ${tests[@]} ; do
        for r in `seq $HI -1 $LO` ; do
            echo $ARGS -o $NAME/comp-$c-r$r -refine $r $ANA ../comp-$c.txt
        done
    done | xargs -P 4 -n 1 -d '\n' bash -c > /dev/null
fi

names=("" "" "ul-inf" "ul-2" "pl-inf" "pl-2")
for c in ${tests[@]} ; do
    cat <<EOF > $NAME/comp-$c.txt
EOF
    for r in `seq $LO 1 $HI` ; do
        grep 'l-2' $NAME/comp-$c-r$r/common/log.txt |\
            sed "s/.*f \([^ ]*\).*2 \([^ ]*\).*f \([^ ]*\).*2 \([^<]*\).*/$r \1 \2 \3 \4/g" >> $NAME/comp-$c.txt
    done

    for i in 2 3 4 5 ; do
        gnuplot -e 'set terminal pdf' \
            -e "set output \"$NAME/comp-$c-${names[$i]}.pdf\"" \
            -e "fit a*x+b \"$NAME/comp-$c.txt\" u (log10(\$1)):(log10(\$$i)) via a,b" \
            -e "plot \"$NAME/comp-$c.txt\" u (log10(\$1)):(log10(\$$i)),a*x+b title sprintf(\"%.2f*x+%.2f\",a,b)" \
            2>/dev/null
    done
done
