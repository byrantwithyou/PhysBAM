#!/bin/bash

CMD="../mpm -symplectic_euler -last_frame 100 75 -use_reflect -scale_E 1 -resolution 32 -friction .1 -3d -threads 8 -float"
$CMD -o vibration/out_vibration -cfl_c .99 -strong_cfl -sound_cfl
$CMD -o vibration/out_vibration_max_stable -max_dt 0.0026 -min_dt .0026
$CMD -o vibration/out_vibration_min_unstable -max_dt 0.0027 -min_dt .0027
$CMD -o vibration/out_vibration_explode -max_dt 0.0029 -min_dt .0029
