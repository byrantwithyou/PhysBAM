#!/bin/sh
DATA_DIR=A
SIM=B
MOTION=C

python grow.py -s $SIM -l 0.30 -n 50 -e 50 -i $DATA_DIR
python graph.py -s $SIM -d $DATA_DIR
python generate_fixed_points.py -d -s $SIM -k $MOTION $DATA_DIR

#python grow.py -s shaking_straight_long_50_50_smaller_perturb -l 0.30 -n 50 -e 50 -i /solver/vol3/hair1/data/body2
#python graph.py -s shaking_straight_long_50_50_smaller_perturb   -d /solver/vol3/hair1/data/body2
#python generate_fixed_points.py
#usage: generate_fixed_points.py [options] <data directory>
#
#options:
#  -h, --help            show this help message and exit
#  -m MODEL_NAME, --model_name=MODEL_NAME
#  -s SIM_NAME, --sim_name=SIM_NAME
#  -d, --deformable      
#  -k KEYFRAME_MOTION, --keyframe_motion=KEYFRAME_MOTION

