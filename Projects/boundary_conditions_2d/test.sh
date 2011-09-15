#!/bin/bash

R="4 12 20 28 36 44 52 60 68 76 84 92 100 116 132 148 164"

#EXTRA_ARGS="-viscosity 10 -boundary_offset .0 -dt 1 -framerate 200 -last_frame 1 -levelset_viscosity 4"
EXTRA_ARGS="$@"

for r in $R ; do
    rm -rf o-$r
#    nice boundary_conditions_2d -resolution $r -o o-$r $EXTRA_ARGS 2>/dev/null | grep "W@Z-" | sed 's/W@Z-//' | sed 's/[()]//g' > nt-nop-$r &
    nice test_bed -o o-$r $EXTRA_ARGS -resolution $r 2>/dev/null | grep "W@Z-" | sed 's/W@Z-//' | sed 's/[()]//g' > nt-nop-$r &
done
wait
#pattern=`printf "%s\\|" "$@"`"WILDCARD"
firstres=`echo $R | awk '{print $1}'`
for l in `seq 1 $(wc -l < nt-nop-$firstres)` ; do
    name=`head -n $l < nt-nop-$firstres | tail -n 1 | awk '{print $1}'`
#    [ "$pattern" != "WILDCARD" ] && ! echo "$name" |grep "$pattern"  && continue
    echo
    echo ===== "$name" =====
    echo
    for r in $R ; do
        head -n $l < nt-nop-$r | tail -n 1
    done | tee nt-in-$l | order-of-accuracy 2 3
done
#cat nt-nop-{?,??,???} | 
