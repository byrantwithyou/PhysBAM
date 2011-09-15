#!/usr/bin/bash

OUTPUT_DIRECTORY=`mktemp -d`
INPUT=$1
FRAME=$2
OUTPUT_FILE=$3

physbam2gnuplot_nocona -dimension 2 -pressure -start_frame $FRAME -last_frame $FRAME $INPUT -o $OUTPUT_DIRECTORY
X=`particle_dump_nocona $INPUT/$FRAME/rigid_geometry_particles.gz | head -n1 | sed s,"X = \[",,g | sed s,\ .*$,,g`
Y=`particle_dump_nocona $INPUT/$FRAME/rigid_geometry_particles.gz | head -n1 | sed s,"X = \[",,g | sed s,^[^\ ]*\ ,,g | sed s,\ .*$,,g`

pushd $OUTPUT_DIRECTORY
~/bin/gnuplot <<EOF
set contour base;set cntrparam levels incremental 2,(26./53.),132;unset surface;set view map
set table 'contour.dat'
splot "pressure.$FRAME" using 1:2:3 with lines
unset table

set terminal png enhanced size 4000,800 font "Vera,24"
set output "$OUTPUT_FILE"

unset key
set xrange [0:1]
set yrange [0:.2]
set size ratio .2
set xtics 0,.25,1
set ytics 0,.1,.2
set mxtics .05
set mytics .02
set parametric
set border -1 lw 3
plot [0:2*pi] 'contour.dat' with lines lw 3, $X+.0505*sin(t),$Y+.0505*cos(t) with filledcurves closed lt 3
EOF

popd
cp $OUTPUT_DIRECTORY/$OUTPUT_FILE .
rm -rf $OUTPUT_DIRECTORY
