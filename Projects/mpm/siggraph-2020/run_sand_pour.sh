#!/bin/bash

DIR="out_sand_pour"

../mpm -3d 44 -strong_cfl -sound_cfl -symplectic_euler -use_exp_F -max_dt 1 -use_reflect -threads 8 -last_frame 20 -o $DIR -cfl_F 0.2 -no_implicit_plasticity -friction .1  -fooT3 1 -resolution 10 -fooT2 35  -use_reflect_friction

