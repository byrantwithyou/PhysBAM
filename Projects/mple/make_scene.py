#!/usr/bin/python

import sys

if len(sys.argv)!=5:
    print "Usage: make.scene.py <input> <output> <surface> <particles>"
    exit(-1)

input_file_name=sys.argv[1]
output_file_name=sys.argv[2]
surface_file_name=sys.argv[3]
particles_file_name=sys.argv[4]

input_file=open(input_file_name,"r")
input_lines=input_file.readlines()

surface_file=open(surface_file_name,"r")
surface_lines=surface_file.readlines()

particles_file=open(particles_file_name,"r")
particles_lines=particles_file.readlines()

output_file=open(output_file_name,"w")

for input_line in input_lines:
    if "#emit" in input_line:
        if "sim_surface" in input_line:
            for surface_line in surface_lines:
                output_file.write(surface_line)
        if "sim_particles" in input_line:
            for particles_line in particles_lines:
                output_file.write(particles_line)
    else:
        output_file.write(input_line)

input_file.close()
surface_file.close()
particles_file.close()
output_file.close()








