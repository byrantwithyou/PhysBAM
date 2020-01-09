#!/bin/bash

DIR="out_single_particle_stability_test"

../mpm 75 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -last_frame 300 -o $DIR -cfl_F 0.2 -no_affine 
