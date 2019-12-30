#!/bin/bash

DIR="out_cube_pour"

./mpm -3d 73 -strong_cfl -sound_cfl -symplectic_euler -use_exp_F -max_dt 1 -use_reflect -threads 8 -last_frame 200 -o $DIR -cfl_F 0.2 -friction .1 -dump_collisions
