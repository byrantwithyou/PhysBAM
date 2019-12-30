#!/bin/bash

DIR="out_friction"
EXP="-strong_cfl -sound_cfl -symplectic_euler -use_exp_F -max_dt 1 -use_reflect"

../mpm -3d 74 $EXP -o $DIR -last_frame 100 -scale_E .1 -resolution 32 -friction .1 -threads 8

## Fluid Source zero friction, three bc types 

./mpm -3d 49 $EXP -threads 8 -last_frame 100 -o $DIR"Test_3d_49""seperate" -cfl_F 0.2 -dilation_only -separate

./mpm -3d 49 $EXP -threads 8 -last_frame 100 -o $DIR"Test_3d_49""slip" -cfl_F 0.2 -dilation_only -slip

./mpm -3d 49 $EXP -threads 8 -last_frame 100 -o $DIR"Test_3d_49""stick" -cfl_F 0.2 -dilation_only -stick

## Cube Pouring down

./mpm -3d 73 $EXP -threads 8 -last_frame 200 -o $DIR"Test_3d_73" -cfl_F 0.2 -friction .1 -dump_collisions

## Dam Break

./mpm -3d 34 $EXP -threads 8 -last_frame 30 -o $DIR"Test_3d_34" -cfl_F 0.2 -no_implicit_plasticity -friction .1 -fooT1 10 -fooT2 1000 -fooT3 3 -scale_E 10


## Hourglass
./mpm -3d 36 $EXP -threads 8 -last_frame 30 -o $DIR"Test_3d_36" -cfl_F 0.2 -no_implicit_plasticity -friction .1 -dump_collisions 

## Sand Pour

./mpm -3d 44 $EXP -threads 8 -last_frame 20 -o $DIR"Test_3d_44" -cfl_F 0.2 -no_implicit_plasticity -friction .1  -fooT3 1 -resolution 10 -fooT2 35   



