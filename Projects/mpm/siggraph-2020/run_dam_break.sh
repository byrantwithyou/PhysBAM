#!/bin/bash

DIR="out_dam_break"

./mpm -3d 34 -strong_cfl -sound_cfl -symplectic_euler -use_exp_F -max_dt 1 -use_reflect -threads 8 -last_frame 30 -o $DIR -cfl_F 0.2 -no_implicit_plasticity -friction .1 -fooT1 10 -fooT2 1000 -fooT3 3 -scale_E 10
