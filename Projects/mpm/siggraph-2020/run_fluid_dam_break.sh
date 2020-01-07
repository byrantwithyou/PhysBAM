#!/bin/bash

DIR="out_fluid_dam_break"

../mpm 22 -strong_cfl -cfl_F 0.2 -sound_cfl -symplectic_euler -use_reflect -last_frame 300 -dilation_only -friction 0 -separate -use_reflect_friction -o $DIR"_not_stable" -max_dt 1