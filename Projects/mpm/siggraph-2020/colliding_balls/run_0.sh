#!/bin/bash

LF=200

(

CMD="../../mpm -3d 11 -max_dt 1 -symplectic_euler -last_frame $LF -no_implicit_plasticity -threads 1 -scale_E 100 -float -resolution 64 -use_reflect -use_reflect_friction -friction 0.1 -frame_dt .002 -verbose_cfl"

function F()
{
    O="out_$1"
    shift
    rm -rf $O
    $CMD -o $O $@ | grep frame
    mkdir $O/geo
    ../../../../Tools/mpm2partio/mpm2partio -3d $O $O/geo/geo.%04d.bgeo
}

F no_b -no_affine_cfl &
F pic -no_affine &
F no_f_c &
F no_c -strong_cfl -cfl_F 0.2 &
F no_f -sound_cfl -cfl_c .9 &
F no_f_13 -sound_cfl -cfl_c 1.3 &
F no_f_15 -sound_cfl -cfl_c 1.5 &
F no_f_20 -sound_cfl -cfl_c 2.0 &
F ok -sound_cfl -cfl_c .9 -strong_cfl -cfl_F 0.2 &

wait

)
