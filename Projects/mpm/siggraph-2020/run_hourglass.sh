#!/bin/bash

DIR="out_hourglass"

../mpm -3d 36 -max_dt 1 -sound_cfl -symplectic_euler -strong_cfl -cfl_F 0.2 -last_frame 240 -no_implicit_plasticity -threads 8 -dump_collisions -scale_E .1 -float -o $DIR -sand_color -resolution 128
