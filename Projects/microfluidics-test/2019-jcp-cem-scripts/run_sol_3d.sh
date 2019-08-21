#!/bin/bash

NAME=sol-3d

export OPENBLAS_NUM_THREADS=1

# tests=(simple grid20 rgrid0 rgrid1 voronoi-s4 voronoi-s15)
RES=8
PATTERN="/disk2/cache/cache%d"
ARGS="../fem -3d -q -dump_sol -refine $RES -cache $PATTERN -mu 8.9e-4"

rm -rf $NAME
mkdir -p $NAME

c="voronoi-s4";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 16 -cache $PATTERN -cache_size 300 ../$c.txt > ../out.txt

c="voronoi-s15";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 10 -cache $PATTERN -cache_size 400 ../$c.txt > ../out.txt

c="rgrid0";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 10 -cache $PATTERN -cache_size 500 ../$c.txt > ../out.txt

c="rgrid1";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 10 -cache $PATTERN -cache_size 500 ../$c.txt > ../out.txt

c="grid20";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 16 -cache $PATTERN -cache_size 400 ../$c.txt > ../out.txt

c="simple";
rm /disk2/cache/*
$ARGS -o $NAME/$c-r$RES -threads 10 -cache $PATTERN -cache_size 150 ../$c.txt > ../out.txt
