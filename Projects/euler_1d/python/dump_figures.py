#!/usr/bin/python

from subprocess import *
import os
import sys

data_path="/disk2/compressible/sims/smooth_flow_with_cfl_multiplication_factor_unit_domain" #Final results for smoothflow
#data_path="/n/lie/disk2/Compressible/output_python_with_enthalpy_machnumber/" #Final results for sod
#data_path="/n/viscosity/data/kwatra/PhysBAM/Projects/euler_1d/output_python_smoothflow_with_cfl_multiplication_factor"
#data_path="/n/viscosity/data/kwatra/PhysBAM/Projects/euler_1d/output_python_2/"
examples_sod={
        1:"Sod shock tube",
        4:"Lax's shock tube",
        5:"Strong shock tube",
        6:"Two symmetric rarefaction waves",
        7:"Mach 3 shock test",
        8:"High mach flow test",
        9:"Two shocks",
        10:"Interaction of blast waves (Bang Bang)"}

frame_number_values={
        1:15,
        4:12,
        5:5,
        6:15,
        7:9,
        8:10,
        9:15,
        10:10}

#example_kind_values=["sod","smoothflow"]
example_kind_values=["smoothflow"]

test_number_values={
        "sod":examples_sod,
        "smoothflow":{1:""}}

variables=["density","centered_velocities","pressure","machnumber","enthalpy","internal_energy","entropy"]

#directories_format_default="Test_%d__Resolution_400_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit Test_%d__Resolution_1600_eno_order-2_rk_order-3_explicit"
directories_format_default="Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_implicit_rk_density_weighted_ENO_semiimplicit"
plot_file_prefix_default="plot_1601_explicit_vs_401_density_weighted_implicit_rk_semiimplicit_"
font_size_default=40

# Fill Sod data
directories_format_dictionary_sod=dict()
plot_file_prefix_dictionalry_sod=dict()
font_size_dictionary_sod=dict()
for test_number in examples_sod.keys():
    directories_format_dictionary_sod[test_number]=[directories_format_default]
    plot_file_prefix_dictionalry_sod[test_number]=[plot_file_prefix_default]
    font_size_dictionary_sod[test_number]=[font_size_default]

    # Add comparison of the various ENO schemes
    if(test_number==1):
        # Old ENO
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_explicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_1601_explicit_vs_401_old_ENO_")
        font_size_dictionary_sod[test_number].append(font_size_default)
        # New ENO (density weighted velocity average)
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_density_weighted_ENO_explicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_1601_explicit_vs_401_new_ENO_density_weighted_")
        font_size_dictionary_sod[test_number].append(font_size_default)
        # New ENO (usual average)
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_velocity_weighted_ENO_explicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_1601_explicit_vs_401_new_ENO_velocity_weighted_")
        font_size_dictionary_sod[test_number].append(font_size_default)

    # Add comparison of implict inside vs outside rungekutta
    if((test_number==1) | (test_number==5)):
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_401_eno_order-2_rk_order-3_semiimplicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_401_semiimplicit_")
        font_size_dictionary_sod[test_number].append(font_size_default)
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_401_eno_order-2_rk_order-3_implicit_rk_semiimplicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_401_implicit_rk_semiimplicit_")
        font_size_dictionary_sod[test_number].append(font_size_default)

        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_semiimplicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_1601_explicit_vs_401_semiimplicit_")
        font_size_dictionary_sod[test_number].append(font_size_default)
        directories_format_dictionary_sod[test_number].append("Test_%d__Resolution_1601_eno_order-2_rk_order-3_explicit Test_%d__Resolution_401_eno_order-2_rk_order-3_implicit_rk_semiimplicit")
        plot_file_prefix_dictionalry_sod[test_number].append("plot_1601_explicit_vs_401_implicit_rk_semiimplicit_")
        font_size_dictionary_sod[test_number].append(font_size_default)


# Fill Smoothflow data
font_size_smoothflow=20
directories_format_dictionary_smoothflow=dict()
plot_file_prefix_dictionalry_smoothflow=dict()
font_size_dictionary_smoothflow=dict()
directories_format_dictionary_smoothflow[1]=list()
plot_file_prefix_dictionalry_smoothflow[1]=list()
font_size_dictionary_smoothflow[1]=list()
#directories_format_dictionary_smoothflow[1].append(directories_format_default)
#plot_file_prefix_dictionalry_smoothflow[1].append(plot_file_prefix_default)

# Comparing the 3 resolutions of semi-implicit at cfl 3. Baseline explicit at 3200
directories_format_dictionary_smoothflow[1].append("Test_%d__Resolution_3200_eno_order-2_rk_order-3_explicit Test_%d__Resolution_200_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_400_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_800_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit")
plot_file_prefix_dictionalry_smoothflow[1].append("plot_3200_explicit_200_400_800_effective_sound_speed_cfl_3_density_weighted_ENO_semiimplicit_")
font_size_dictionary_smoothflow[1].append(font_size_smoothflow)

# Comparing the 3 resolutions of semi-implicit at cfl 3. Baseline explicit at 3200
directories_format_dictionary_smoothflow[1].append("Test_%d__Resolution_3200_eno_order-2_rk_order-3_explicit Test_%d__Resolution_800_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_1600_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_3200_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit") 
plot_file_prefix_dictionalry_smoothflow[1].append("plot_3200_explicit_800_1600_3200_effective_sound_speed_cfl_3_density_weighted_ENO_semiimplicit_")
font_size_dictionary_smoothflow[1].append(font_size_smoothflow)

# Comparing the 5 resolutions of semi-implicit at cfl 3. Baseline explicit at 3200
directories_format_dictionary_smoothflow[1].append("Test_%d__Resolution_3200_eno_order-2_rk_order-3_explicit Test_%d__Resolution_200_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_400_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_800_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_1600_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_3200_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit")
plot_file_prefix_dictionalry_smoothflow[1].append("plot_3200_explicit_200_400_800_1600_3200_effective_sound_speed_cfl_3_density_weighted_ENO_semiimplicit_")
font_size_dictionary_smoothflow[1].append(font_size_smoothflow)

# Comparing results from multplying resolution and effective sound speed cfl by the same amount
directories_format_dictionary_smoothflow[1].append("Test_%d__Resolution_3200_eno_order-2_rk_order-3_cfl_sound_speed_multiple-6.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_32000_eno_order-2_rk_order-3_cfl_sound_speed_multiple-60.000000_density_weighted_ENO_semiimplicit Test_%d__Resolution_320000_eno_order-2_rk_order-3_cfl_sound_speed_multiple-600.000000_density_weighted_ENO_semiimplicit")
plot_file_prefix_dictionalry_smoothflow[1].append("plot_3200_32000_320000_varying_effective_sound_speed_cfl_density_weighted_ENO_semiimplicit_")
font_size_dictionary_smoothflow[1].append(font_size_smoothflow)


# main dumping

directories_format_dictionary={
        "sod":directories_format_dictionary_sod,
        "smoothflow":directories_format_dictionary_smoothflow}
plot_file_prefix_dictionalry={
        "sod":plot_file_prefix_dictionalry_sod,
        "smoothflow":plot_file_prefix_dictionalry_smoothflow}
font_size_dictionary={
        "sod":font_size_dictionary_sod,
        "smoothflow":font_size_dictionary_smoothflow}

for example_kind in example_kind_values:
    directory_path="%s/%s"%(data_path,example_kind)
    if(example_kind=="sod"):
        plot_width=6
        verbosity=0
        xylabel=""
        variables=["density","centered_velocities","pressure","internal_energy"]
    elif(example_kind=="smoothflow"):
        plot_width=3
        verbosity=3
        xylabel="--xlabel --ylabel"
        variables=["pressure"]
    test_number_dictionary=test_number_values[example_kind]
    test_numbers=test_number_dictionary.keys()
    test_numbers.sort()
    for test_number in test_numbers:
        directories_format_list=directories_format_dictionary[example_kind][test_number]
        plot_file_prefix_list=plot_file_prefix_dictionalry[example_kind][test_number]
        font_size_list=font_size_dictionary[example_kind][test_number]
        list_index=0
        for directories_format in directories_format_list:
            plot_file_prefix=plot_file_prefix_list[list_index]
            font_size=font_size_list[list_index]
            list_index+=1
            directories=directories_format.replace("Test_%d","%s/Test_%d"%(directory_path,test_number))
            print "\n\n"
            print "______________________"
            print "directories=",directories
            print "______________________"
            print "\n\n"
            for variable in variables:
                frame_number=frame_number_values[test_number]
                plot_file_name=plot_file_prefix+"%s_%05d"%(variable,frame_number)
                command_string="""scripts/plot.py -v %s -f %d --pw %d --fontsize %d --minxrange 0 --maxxrange 1 --verbosity %d %s --double --path %s -n %d -c "wf %s q" %s"""%(variable,frame_number,plot_width,font_size,verbosity,xylabel,directory_path,test_number,plot_file_name,directories)
                print "\n\n\n"
                print "running: ",command_string
                os.system(command_string)
                print "\n\n\n"
