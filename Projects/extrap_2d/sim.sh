#!/bin/bash
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -use_corotated -o fig-corot >&/dev/null &
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .4 -cutoff .001 -ether_drag 20 -framerate 48 -o fig-neo >&/dev/null &
extrap_2d -stiffen .1 -dampen 0 100 -image_size 30 -sigma_range 3 -dt 1 -last_frame 500 -poissons_ratio .45 -cutoff .001 -ether_drag 20 -framerate 48 -use_corotated -o fig-corot-inv -parameter 1 >&/dev/null &
wait
cp eps-cap/camera_script fig-corot
cp eps-cap/camera_script fig-neo
cp eps-cap/camera_script fig-corot-inv
$PHYSBAM/Projects/opengl_2d/opengl_2d$SUFFIX fig-corot -offscreen -w 640 -h 480 -so corot-%03d.eps -start_frame 0 -stop_frame 500 -keys 'wW=================' >&/dev/null &
$PHYSBAM/Projects/opengl_2d/opengl_2d$SUFFIX fig-neo -offscreen -w 640 -h 480 -so neo-%03d.eps -start_frame 0 -stop_frame 500 -keys 'wW=================' >&/dev/null &
$PHYSBAM/Projects/opengl_2d/opengl_2d$SUFFIX fig-corot-inv -offscreen -w 640 -h 480 -so corot-inv-%03d.eps -start_frame 0 -stop_frame 500 -keys 'wW=================' >&/dev/null &
wait
rm */*/*eps
mv corot-inv-*eps eps-cap/cap-corot-inv/
mv corot-*eps eps-cap/cap-corot/
mv neo-*eps eps-cap/cap-neo/

(
    cd eps-cap
    for k in neo corot corot-inv ; do for i in {001..500} ; do compos.pl < cap-$k/$k-$i.eps > $k-$i.eps ; convert -quality 100 $k-$i.eps $k-$i.jpg ; echo $i ; done & done
    wait
)
