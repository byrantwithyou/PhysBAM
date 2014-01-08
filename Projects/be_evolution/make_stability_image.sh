#!/bin/bash

SIZE=480

rm -rf image-ours image-vanilla
mkdir image-ours image-vanilla

DEF_OPTS="-2d 102 -newton_it 500 -stiffen 10 -o $PWD/ignore-out -nolog -v 0"
VAN_OPTS="-kry_tol 1e-5 -use_vanilla_newton -mr"
OUR_OPTS=
SLAVE="/home/craig/PhysBAM/Tools/batch/slave"
EXEC="$PWD/be_evolution"

for i in {000..120} ; do
    dt=`perl -e '$a="'$i'";$b=2**(-(120-$a)/10);print "$b\n";'`
    $SLAVE -- nice $EXEC -image_size "$SIZE $SIZE" $DEF_OPTS $OUR_OPTS -image_file $PWD/image-ours/image_$i.png -raw_image_file $PWD/image-ours/raw_image_$i.txt -dt $dt
    $SLAVE -- nice $EXEC -image_size "$SIZE $SIZE" $DEF_OPTS $VAN_OPTS -image_file $PWD/image-vanilla/image_$i.png -raw_image_file $PWD/image-vanilla/raw_image_$i.txt -dt $dt
done


