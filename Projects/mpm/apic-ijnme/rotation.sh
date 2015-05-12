#!/bin/bash

ARGS="1 -max_dt 1e-2 -newton_tolerance 1e-4 -regular_seeding -print_stats -last_frame 120"
if ! : ; then
../mpm -affine -midpoint -o rotation-apic $ARGS >/dev/null &
../mpm -affine -o rotation-apic-be $ARGS >/dev/null &
../mpm -midpoint -o rotation-pic $ARGS >/dev/null &
../mpm -o rotation-pic-be $ARGS >/dev/null &
../mpm -flip 1 -midpoint -o rotation-flip $ARGS >/dev/null &
../mpm -flip 1 -o rotation-flip-be $ARGS >/dev/null &
../mpm -flip .99 -midpoint -o rotation-flip-99 $ARGS >/dev/null &
../mpm -flip .99 -o rotation-flip-99-be $ARGS >/dev/null &
../mpm -flip .95 -midpoint -o rotation-flip-95 $ARGS >/dev/null &
../mpm -flip .95 -o rotation-flip-95-be $ARGS >/dev/null &
wait
fi

for i in {apic,pic,flip,flip-99,flip-95}{,-be} ; do
    grep angular rotation-$i/common/log.txt | grep 'particle state angular' | awk '{print $5 " " $7}' | sed 's/[()]//g' > am-$i.txt
done

./plot-am.m \
apic "APIC Midpoint Rule"  r- \
apic-be "APIC Backward Euler"  r-- \
pic "PIC Midpoint Rule"  k- \
pic-be "PIC Backward Euler"  k-- \
flip "FLIP Midpoint Rule"  c- \
flip-be "FLIP Backward Euler"  c-- \
flip-99 "FLIP-99 Midpoint Rule"  g- \
flip-99-be "FLIP-99 Backward Euler" g-- \
flip-95 "FLIP-95 Midpoint Rule"  b- \
flip-95-be "FLIP-95 Backward Euler" b--
