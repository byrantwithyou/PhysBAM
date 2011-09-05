export outdir=/n/inner/disk2/daily_backups/`date +%Y-%m-%d_%H-%M`
echo $outdir
mkdir -p $outdir
./solids_3d_nocona -fully_implicit -dt 1 -o $outdir/Test_19_param0 -parameter 0 -use_axial -pull_in -combined_collisions -combined_repulsions -no_friction 19

./solids_3d_nocona -fully_implicit -dt 1 -o $outdir/Test_19_param1 -parameter 1 -use_axial -pull_in -combined_collisions -combined_repulsions -no_friction 19

./solids_3d_nocona -fully_implicit -dt 1 -o $outdir/Test_19_param2 -parameter 2 -use_axial -pull_in -combined_collisions -combined_repulsions -no_friction 19

./solids_3d_nocona -fully_implicit -dt 1 -o $outdir/Test_19_param3 -parameter 3 -use_axial -pull_in -combined_collisions -combined_repulsions -no_friction 19

./solids_3d_nocona -fully_implicit -dt 1 -o $outdir/Test_5 -use_axial -pull_in -combined_collisions -combined_repulsions -no_friction 5