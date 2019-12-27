#!/bin/bash

DIR="out_friction"
EXP="-strong_cfl -sound_cfl -symplectic_euler"

../mpm -3d 74 $EXP -o $DIR -last_frame 100 -use_reflect -scale_E .1 -resolution 32 -friction .1 -threads 8

