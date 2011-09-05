#spd
./euler_2d_nocona -last_frame 50 -circle -timesplit -resolution 50 -slip 9
./euler_2d_nocona -last_frame 50 -circle -timesplit -resolution 50 -slip 6
./euler_2d_nocona -last_frame 200 -circle -timesplit -resolution 50 -slip -transition_to_incompressible -time_start_transition 1 -time_end_transition 3 6
./euler_2d_nocona -last_frame 50 -circle -resolution 50 -all_verbose -slip -timesplit -last_frame 100 5

#non spd
./euler_2d_nocona -last_frame 50 -circle -timesplit -resolution 50 6

