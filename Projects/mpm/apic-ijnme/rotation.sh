#!/bin/bash

ARGS="1 -max_dt 1e-2 -newton_tolerance 1e-4 -regular_seeding -print_stats -last_frame 120"
if ! : ; then
../mpm -affine -midpoint -o rotation-apic $ARGS >/dev/null &
../mpm -midpoint -o rotation-pic $ARGS >/dev/null &
../mpm -flip 1 -midpoint -o rotation-flip $ARGS >/dev/null &

../mpm -affine -o rotation-apic-be $ARGS >/dev/null &
../mpm -o rotation-pic-be $ARGS >/dev/null &
../mpm -flip 1 -o rotation-flip-be $ARGS >/dev/null &

../mpm -affine -symplectic_euler -o rotation-apic-symp $ARGS -max_dt 1e-3 >/dev/null &
../mpm -symplectic_euler -o rotation-pic-symp $ARGS -max_dt 1e-3 >/dev/null &
../mpm -flip 1 -symplectic_euler -o rotation-flip-symp $ARGS -max_dt 1e-3 >/dev/null &
wait
fi

for i in {apic,pic,flip}{,-be,-symp} ; do
    grep angular rotation-$i/common/log.txt | grep 'particle state angular' | awk '{print $5 " " $7}' | sed 's/[()]//g' > am-$i.txt
done

./plot-am.m \
apic "APIC Midpoint Rule"  r- \
apic-be "APIC Backward Euler"  r-- \
apic-symp "APIC Symplectic Euler"  r-. \
pic "PIC Midpoint Rule"  g- \
pic-be "PIC Backward Euler"  g-- \
pic-symp "PIC Symplectic Euler"  g-. \
flip "FLIP Midpoint Rule"  b- \
flip-be "FLIP Backward Euler"  b-- \
flip-symp "FLIP Symplectic Euler"  b-.
