#!/bin/bash

DIR="out_fluid_pour"

../mpm -3d 49 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 280 -o $DIR -cfl_F 0.2 -dilation_only -separate
