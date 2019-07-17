#!/bin/bash

NAME=illus-domain

ARGS="../fem -q"

REFINE=4
# tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)

rm -rf $NAME
mkdir -p $NAME

FILLED="$ARGS -illus -illus_fill -illus_zoom"
$FILLED -illus_size "512 512" -o $NAME/grid20 -min_corner "-0.02 0.18" -max_corner "0.12 0.32" ../grid20.txt;
$FILLED -illus_size "512 512" -o $NAME/rgrid0 -min_corner "1 0.35" -max_corner "1.65 1" ../rgrid0.txt;
$FILLED -illus_size "512 512" -o $NAME/rgrid1 -min_corner "1 0.35" -max_corner "1.65 1" ../rgrid1.txt;
$FILLED -illus_size "512 512" -o $NAME/voronoi-s4 -min_corner "0 0" -max_corner "0.15 0.15" ../voronoi-s4.txt;
$FILLED -illus_size "512 512" -o $NAME/voronoi-s15 -min_corner "-0.1 -0.04" -max_corner "0.05 0.11" ../voronoi-s15.txt;

$ARGS -illus_size "512 150" -illus -illus_zoom -min_corner "0.73 -0.105" -max_corner "0.78 -0.05" -refine 4 -o $NAME/simple ../simple.txt
$ARGS -illus_size "512 512" -illus_zoom -min_corner "0.73 -0.105" -max_corner "0.78 -0.05" -refine 4 -o $NAME/simple-zoom ../simple.txt

for f in $NAME/*/*.eps ;
do
    echo epstopdf $f --outfile "${f%%.*}.pdf"
done | xargs -P 16 -n 1 -d '\n' bash -c > /dev/null
