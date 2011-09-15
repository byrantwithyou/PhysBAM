#!/usr/bin/python

from subprocess import *
import os
import sys

test_numbers=[2,3,5]
resolutions_test2=[1600]
resolutions_test3=[1600]
resolutions_test5=[1600]
resolution_dictionary=dict()
resolution_dictionary[2]=resolutions_test2
resolution_dictionary[3]=resolutions_test3
resolution_dictionary[5]=resolutions_test5
maxxrange_dictionary=dict()
maxxrange_dictionary[2]=2
maxxrange_dictionary[3]=20
maxxrange_dictionary[5]=3
variables=["pressure"]
masses_dictionary=dict()
masses_dictionary[2]=[1]
masses_dictionary[3]=[6]
masses_dictionary[5]=[1]
frame_numbers_dictionary=dict()
frame_numbers_dictionary[2]=[0,25,50,75,100]
frame_numbers_dictionary[3]=[0,300,600,900,1500,2000]
frame_numbers_dictionary[5]=[0,150,287,320,400]
rigid_dictionary=dict()
rigid_dictionary[2]=True
rigid_dictionary[5]=True
deformable_dictionary=dict()
deformable_dictionary[3]=True
plot_width=6
font_size=20
verbosity=0
xylabel=""
for test_number in test_numbers:
    resolutions=resolution_dictionary[test_number]
    for resolution in resolutions:
        masses=masses_dictionary[test_number]
        for mass in masses:
            for sim_type in ["explicit","semiimplicit"]:
                directory_name="Test_%d__Resolution_%d_%s_mass_%f"%(test_number,resolution,sim_type,mass)
                for variable in variables:
                    frame_numbers=frame_numbers_dictionary[test_number]
                    for frame_number in frame_numbers:
                        plot_file_name="plot_%d_resolution_%d_%s_mass_%s_variable_%s_frame_%d"%(test_number,resolution,sim_type,mass,variable,frame_number)
                        maxxrange=maxxrange_dictionary[test_number]
                        command_string="""plot.py -v %s -f %d --pw %d --maxxrange %d --fontsize %d --verbosity %d %s --double -n %d -c "wf %s q" %s"""%(variable,frame_number,plot_width,maxxrange,font_size,verbosity,xylabel,test_number,plot_file_name,directory_name)
                        if(rigid_dictionary.__contains__(test_number)):
                            command_string+=" --rbody"
                        if(deformable_dictionary.__contains__(test_number)):
                            command_string+=" --deformable"

                        print "running: ",command_string
                        os.system(command_string)
                        #print "directory_name: ",directory_name
                        #print "plot_file_name: ",plot_file_name
