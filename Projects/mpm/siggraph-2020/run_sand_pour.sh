#!/bin/bash

DIR="out_sand_pour"

../mpm -3d 44 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 144 -o $DIR -cfl_F 0.2 -no_implicit_plasticity -friction .1  -fooT3 1 -resolution 144 -fooT2 35  -use_reflect_friction -sand_color -scale_E .1

