#!/bin/bash
# ./render-frame.sh <pov-file.pov> <frame> <output-dir> <H> <W>

F=`printf "%05d" $2`
SCRIPT=$1
B=${SCRIPT%.pov}
FRAME=$2
OUTPUT=$3
HEIGHT=$4
WIDTH=$5

nice ../../pov_render/pov_render $SCRIPT $OUTPUT/$B-$F.pov $FRAME
nice povray +Q8 +I$OUTPUT/$B-$F.pov -O$OUTPUT/$B-$F.png Antialias=On Display=Off +H$HEIGHT +W$WIDTH

