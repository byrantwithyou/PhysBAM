#!/bin/bash

DIR="out_collision_stability"

../mpm 76 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 200 -o $DIR -cfl_F 0.2 -friction .1 -use_reflect_friction -float -scale_E 10 -resolution 128 -fooT1 0.33