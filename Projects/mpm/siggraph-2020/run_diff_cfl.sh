### Results not satisfying

### Both CFL

# ../mpm 11 -max_dt 1 -sound_cfl -symplectic_euler -use_exp_F -strong_cfl -cfl_F 0.2 -last_frame 400 -3d -resolution 64 -frame_dt 0.23 -threads 8 -use_reflect -use_reflect_friction -friction .3 -separate -float

title="Velocity when turning on both CFL conditions"
python3 plot.py Test_3d_11 $title &


### Only Sound CFL

# DIR='Test_3d_11_sndcfl'

title="Velocity when turning on only sound speed CFL"
python3 plot.py Test_3d_11_sndcfl $title &



# ../mpm 11 -o $DIR -max_dt 1 -sound_cfl -symplectic_euler -use_exp_F -last_frame 400 -3d -resolution 64 -frame_dt 0.23 -threads 8 -separate -use_reflect -use_reflect_friction -friction .3 -float

### Only Strong CFL

# DIR='Test_3d_11_strcfl'

# ../mpm 11 -max_dt 1 -symplectic_euler -use_exp_F -strong_cfl -cfl_F 0.2 -last_frame 400 -3d -resolution 64 -frame_dt 0.23 -o $DIR -threads 8 -separate -use_reflect -use_reflect_friction -friction .3 -float

title="Velocity when turning on only strong CFL"
python3 plot.py Test_3d_11_strcfl $title &
