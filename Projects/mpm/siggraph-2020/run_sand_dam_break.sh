#!/bin/bash

DIR="out_sand_dam_break"

RES=128

CMD="../mpm -3d 34 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 30 -cfl_F 0.2 -no_implicit_plasticity -fooT1 10 -fooT2 1000 -fooT3 3 -scale_E .1 -resolution $RES -use_reflect_friction -sand_color -float"

$CMD -o "sand_dam_break/"$DIR"_hi_mu" -friction 1 
$CMD -o "sand_dam_break/"$DIR"_lo_mu" -friction .45 
