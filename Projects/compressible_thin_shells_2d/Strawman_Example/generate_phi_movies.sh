#!/bin/bash
#scons -Q --random -u -j 64 CXX="/opt/icecream/bin/g++"
#rm -rf *png
for phi in extrapolated fixed ghost; do 
    for dir in `ls -d */ | grep Solut`; do
        pushd $dir;
        for frame in `ls -d */ | grep -v common`; do
            cp ${frame}/phi_${phi}.gz ${frame}/density.gz;
        done
        opengl_2d_nocona -keys "<F8><F8>6" -offscreen -w 640 -h 480 -so capture_${phi}.mov . ;
        popd;
    done;
    wait;
done
