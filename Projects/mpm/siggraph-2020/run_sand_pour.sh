#!/bin/bash

DIR="out_sand_pour"

../mpm -3d 44 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 192 -o $DIR -cfl_F 0.2 -no_implicit_plasticity -friction .5  -fooT3 4 -resolution 64 -fooT2 35  -use_reflect_friction -sand_color -scale_E .01 -float

