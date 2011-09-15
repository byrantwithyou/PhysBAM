./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -o 2d_Tests/Test_1 1 > output.txt
./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -setv -o 2d_Tests/Test_1_setv 1 > output.txt

./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -o 2d_Tests/Test_2 2 > output.txt
./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -setv -o 2d_Tests/Test_2_setv 2 > output.txt

./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -o 2d_Tests/Test_3 3 > output.txt
./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -setv -o 2d_Tests/Test_3_setv 3 > output.txt

./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -o 2d_Tests/Test_4 4 > output.txt
./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -setv -o 2d_Tests/Test_4_setv 4 > output.txt

./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -o 2d_Tests/Test_7 7 > output.txt
./framerate_rigid_deformable_nocona -stiffen 10 -dampen 0 -2d -setv -o 2d_Tests/Test_7_setv 7 > output.txt

./framerate_rigid_deformable_nocona -dampen 0 -2d -side_panels 5 -o 2d_Tests/Test_9 9 > output.txt
./framerate_rigid_deformable_nocona -dampen 0 -setv -2d -side_panels 5 -o 2d_Tests/Test_9_setv 9 > output.txt
