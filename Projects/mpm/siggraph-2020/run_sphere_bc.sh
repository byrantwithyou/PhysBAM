#!/bin/bash

DIR="out_sphere_bc"

DATA="../../../Public_Data"

RES=32
sph=9k


CMD="../mpm -sph $sph -d $DATA -3d 76 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 200 -cfl_F 0.2 -fooT1 5 -scale_E 2500 -resolution $RES -use_reflect_friction"


$CMD -o "sphere_bc/"$DIR"_""no_mu" -separate -friction 0 
$CMD -o "sphere_bc/"$DIR"_""mu" -separate -friction .3 
$CMD -o "sphere_bc/"$DIR"_""slip" -slip  
$CMD -o "sphere_bc/"$DIR"_""stick" -stick  
