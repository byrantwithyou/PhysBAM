#!/bin/bash

NAME=conv-taylor

ARGS="../mpm_mac 16 -last_frame 1 -frame_dt 1 -mu 2 -scale_mass 3"

LO=32
HI=256
SKIP=16
rm -rf $NAME
mkdir -p $NAME

for r in `seq $HI -$SKIP $LO` ; do
    dt=`perl -e "print (1.0/$r)"`
    DT="-max_dt $dt -min_dt $dt"
    echo $ARGS $DT -no_affine -o $NAME/pic-$r -resolution $r
    echo $ARGS $DT -affine -o $NAME/apic-$r -resolution $r
    echo $ARGS $DT -no_affine -flip 1 -o $NAME/flip-$r -resolution $r
done | xargs -P 16 -n 1 -d '\n' bash -c

for t in grid particle ; do
    for s in pic apic flip ; do
        rm -f $NAME/$t-$s.txt
        echo "res linf l2" >> $NAME/$t-$s.txt
        for r in `seq $LO $SKIP $HI` ; do
            echo -n "$r " >> $NAME/$t-$s.txt
            grep "$t velocity error" $NAME/$s-$r/common/log.txt >> $NAME/$t-$s.txt
        done
        sed -i -e 's/ .*L-inf: / /' -e 's/ L-2://' -e 's/ N: .*//' $NAME/$t-$s.txt
    done
done
sed -e 's/xxx/pic/' -e 's/XXX/PIC/' conv_taylor_plot.tex  > $NAME/plot-pic.tex 
sed -e 's/xxx/apic/' -e 's/XXX/APIC/' conv_taylor_plot.tex  > $NAME/plot-apic.tex 
sed -e 's/xxx/flip/' -e 's/XXX/FLIP/' conv_taylor_plot.tex  > $NAME/plot-flip.tex 

for t in grid particle ; do
    for s in pic apic flip ; do
        A=`./conv_regression.pl conv-taylor/$t-$s.txt conv-taylor/$t-$s-regres.txt`
        read oa ob <<< "$A"
        sed -i -e "s/$t-order-inf/$oa/" -e "s/$t-order-2/$ob/" $NAME/plot-$s.tex 
    done
done

for s in pic apic flip ; do
    pdflatex $NAME/plot-$s.tex 
done
