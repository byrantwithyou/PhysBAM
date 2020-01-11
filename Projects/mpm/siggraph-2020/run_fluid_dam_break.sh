#!/bin/bash

DIR="out_fluid_dam_break"

CMD="../mpm 22 -strong_cfl -cfl_F 0.2 -sound_cfl -symplectic_euler -use_reflect -last_frame 400 -dilation_only -friction 0 -separate -use_reflect_friction -scale_E 10 -single_particle_cfl -cfl_c .9 -float"

$CMD -cfl_p 0.9 -o $DIR"_ours"
$CMD -cfl_p 1.3 -o $DIR"_max_stable"
$CMD -cfl_p 1.4 -o $DIR"_min_unstable"
$CMD -cfl_p 2.0 -o $DIR"_unstable"

