#!/bin/bash
# ./render.sh <pov-file.pov> <frame_first> <frame_last> <output-dir> <H> <W>

SCRIPT=$1
FIRST=$2
LAST=$3
OUTPUT=$4
HEIGHT=$5
WIDTH=$6

if [ ! -d $OUTPUT ]; then
    mkdir $OUTPUT
fi

for f in `seq $FIRST $LAST` ; do
    echo ./render-frame.sh $SCRIPT $f $OUTPUT $HEIGHT $WIDTH
done | xargs -P 16 -n 1 -d '\n' bash -c

