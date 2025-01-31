#!/bin/bash

DIR="out_fluid_pour"

CMD="../mpm -3d 49 -fooT1 0.1 -strong_cfl -sound_cfl -symplectic_euler -max_dt 1 -use_reflect -threads 8 -last_frame 200 -cfl_F 0.2 -dilation_only -separate -resolution 64 -single_particle_cfl -cfl_p .99 -float"

 $CMD -o "fluid_pour/"$DIR
 
