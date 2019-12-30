#!/bin/bash

DIR="out_fluid_pour"

./mpm -3d 49 -strong_cfl -sound_cfl -symplectic_euler -use_exp_F -max_dt 1 -use_reflect -threads 8 -last_frame 100 -o $DIR -cfl_F 0.2 -dilation_only -separate
