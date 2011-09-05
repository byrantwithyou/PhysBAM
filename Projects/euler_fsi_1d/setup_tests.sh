#fluids only
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  1
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit 1
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 1
#solid-fluids non-spd
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit 2
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 -no_walls 2
#spd
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit -slip 2
./euler_fsi_1d_nocona_debug -resolution 200 -last_frame 100  -timesplit -transition_to_incompressible -time_start_transition .6 -time_end_transition .9 -slip -no_walls 2

if [ -d Standard_Tests/datum ]; then
    echo "Standard_Tests/datum directory already exists"
else
    echo "creating Standard_Tests/datum directory";
    mkdir Standard_Tests/datum
fi

for folder in {Test_1__Resolution_200_explicit,Test_1__Resolution_200_semiimplicit,Test_1__Resolution_200_semiimplicit_transition_incompressible,Test_2__Resolution_200_semiimplicit_mass_1.000000,Test_2__Resolution_200_semiimplicit_transition_incompressible_mass_1.000000,Test_2__Resolution_200_semiimplicit_slip_mass_1.000000,Test_2__Resolution_200_semiimplicit_slip_transition_incompressible_mass_1.000000}; do
    mv Standard_Tests/$folder Standard_Tests/datum/$folder;
    ln -s datum/$folder Standard_Tests/${folder}_exact;
done
