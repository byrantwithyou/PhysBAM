#fluids only
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  1
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit 1
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 1
#solid-fluids non-spd
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit 2
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 -no_walls 2
#spd
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit -slip 2
./euler_fsi_1d_nocona -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 -slip -no_walls 2

for folder in {Test_1__Resolution_200_explicit,Test_1__Resolution_200_semiimplicit,Test_1__Resolution_200_semiimplicit_transition_incompressible,Test_2__Resolution_200_semiimplicit_mass_1.000000,Test_2__Resolution_200_semiimplicit_transition_incompressible_mass_1.000000,Test_2__Resolution_200_semiimplicit_slip_mass_1.000000,Test_2__Resolution_200_semiimplicit_slip_transition_incompressible_mass_1.000000}; do
    opengl_1d_nocona -double -keys az\!\&V Standard_Tests/$folder;
done

# MPI TESTING

#fluids only
mpirun -np 3 xterm -e ./euler_fsi_1d_nocona -no_preconditioner -all_verbose -xprocs 2 -resolution 200 -last_frame 100  1
mpirun -np 3 xterm -e ./euler_fsi_1d_nocona -no_preconditioner -all_verbose -xprocs 2 -resolution 200 -last_frame 100  -timesplit 1
mpirun -np 3 xterm -e ./euler_fsi_1d_nocona -no_preconditioner -all_verbose -xprocs 2 -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 1
#spd
mpirun -np 3 xterm -e ./euler_fsi_1d_nocona -no_preconditioner -all_verbose -xprocs 2 -resolution 200 -last_frame 100  -timesplit -slip 2
mpirun -np 3 xterm -e ./euler_fsi_1d_nocona -no_preconditioner -all_verbose -xprocs 2 -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 -slip -no_walls 2

for folder in {Test_1__Resolution_200_explicit,Test_1__Resolution_200_semiimplicit,Test_1__Resolution_200_semiimplicit_transition_incompressible,Test_2__Resolution_200_semiimplicit_slip_mass_1.000000,Test_2__Resolution_200_semiimplicit_slip_transition_incompressible_mass_1.000000}; do
    fluids_mpi_merge_nocona -double -compressible -solid_fluid -np 3 -1d Standard_Tests/$folder;
    opengl_1d_nocona -double -keys az\!\&V Standard_Tests/$folder;
done
