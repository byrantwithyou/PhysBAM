#!/bin/bash
#scons -Q --random -u -j 64 CXX="/opt/icecream/bin/g++"
#rm -rf *png

for res in 25 50 100 200 400; do
    for scheme in analytic extrapolation gfm new_gfm; do
        ./compressible_thin_shells_1d_nocona -strawman -resolution ${res} -${order} -v 0 &
    done;
    wait;
done

##for res in 25 50 100 200 400 800 1600 3200 6400 12800; do
#res=25
#for vel in 1.000000 0.000000 -1.000000; do
#    for order in 1 2 3; do
#        pushd Strawman_Example/Solution_${vel}_Resolution_${res}_Order_${order}
#            for frame in `ls | grep -v common | grep -v png`; do sh ../../generate_1d_slice.sh . ${frame} 10;done
#        popd
#    done;
#    wait;
#done
