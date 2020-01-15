#!/bin/bash

DIR="out_sphere_bc"

RES=32

CMD="../mpm -3d 76 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 200 -cfl_F 0.2 -fooT1 5 -scale_E 10000 -resolution $RES -use_reflect_friction -float"

$CMD -o $DIR"_no_mu" -separate -friction 0
$CMD -o $DIR"_mu" -separate -friction .3 
$CMD -o $DIR"_slip" -slip 
$CMD -o $DIR"_stick" -stick 
