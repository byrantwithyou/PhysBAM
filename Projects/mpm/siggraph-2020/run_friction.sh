#!/bin/bash

DIR="out_friction"

../mpm -3d 74 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -o "friction/"$DIR -last_frame 340 -scale_E .5 -resolution 64 -friction .1 -threads 8 -use_reflect_friction -cfl_F 0.2 -frame_dt 0.01

