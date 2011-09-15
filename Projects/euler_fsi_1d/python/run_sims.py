#!/usr/bin/python

def Run_Euler_1d(binary,output_directory,resolution,eno_scheme,eno_order,rk_order,cfl,cfl_sound_speed_multiple,implicit_rk,timesplit,example_kind,test_number):
    """Runs the euler_1d test with the given parameters"""

    simulation_directory="%s/%s/Test_%d__Resolution_%d"%(output_directory,example_kind,test_number,resolution)
    simulation_directory+="_eno_order-%d"%eno_order
    simulation_directory+="_rk_order-%d"%rk_order
    if(cfl_sound_speed_multiple>0): simulation_directory+="_cfl_sound_speed_multiple-%f"%cfl_sound_speed_multiple
    if(implicit_rk): simulation_directory+="_implicit_rk"
    if(eno_scheme==2): simulation_directory+="_density_weighted_ENO"
    elif(eno_scheme==3): simulation_directory+="_velocity_weighted_ENO"
    if(timesplit): simulation_directory+="_semiimplicit"
    else: simulation_directory+="_explicit"

    command_string=binary
    command_string+=(" -resolution %d"%resolution)
    command_string+=(" -eno_scheme %d"%eno_scheme)
    command_string+=(" -eno_order %d"%eno_order)
    command_string+=(" -rk_order %d"%rk_order)
    command_string+=(" -cfl %f"%cfl)
    if(cfl_sound_speed_multiple>0): command_string+=" -cfl_sound_speed_multiple %f"%cfl_sound_speed_multiple

    if(implicit_rk): command_string+=(" -implicit_rk")
    if(timesplit): command_string+=(" -timesplit")
    command_string+=(" -o %s"%simulation_directory)
    command_string+=(" -%s"%example_kind)
    command_string+=(" %d"%test_number)
    return (command_string,simulation_directory)

def Run_Simulation(command_string,simulation_directory,run_erring_sims,log_file):
    if(os.path.exists("%s/done"%simulation_directory)):
        log_message="Skipping. Simulation in directory %s already complete\n\n"%simulation_directory
        Print_And_Log(log_file,log_message)
        return
    elif(os.path.exists("%s/error"%simulation_directory)):
        log_message="Simulation in directory %s did not complete succesfully in last run\n"%simulation_directory
        Print_And_Log(log_file,log_message)
        if(run_erring_sims):
            log_message="Re-running:\n"
            Print_And_Log(log_file,log_message)
        else:
            log_message="Skipping\n\n"
            Print_And_Log(log_file,log_message)
            return

    command_string_split=command_string.split(" ");
    log_message="running command %s\n"%command_string
    Print_And_Log(log_file,log_message)

    # Run the process
    process=Popen(command_string_split);
    process.wait()
    if(process.returncode==None):
        log_message="Error: process hasn't terminated\n\n"
    elif(process.returncode==0):
        log_message="Process completed succesfully\n\n"
        done_file=open("%s/done"%simulation_directory,"w")
        done_file.close()
    else:
        log_message="Process did not completed succesfully\n\n"
        error_file=open("%s/error"%simulation_directory,"w")
        error_file.close()
    Print_And_Log(log_file,log_message)


def Print_And_Log(log_file,log_message):
    """Prints the log_message and also writes into log_file"""

    print log_message
    log_file.write(log_message)
    log_file.flush()

def Find_Last_Index_Of_File_Sequence(file_format,index_range):
    last_found_index=-1
    for index in index_range:
        file_name=file_format%index
        if(os.path.exists(file_name)):
            last_found_index=index
        else:
            if(last_found_index==-1):
                print("No file in range found")
            return last_found_index
    else:
        print("All files in range exist")
        return index_range[1]

from subprocess import *
import sys
import os

from optparse import OptionParser

usage="usage: %prog [options] <binary> <output_directory>"
parser=OptionParser(usage)
parser.add_option("--run_erring_sims",action="store_true",dest="run_erring_sims",default=False,
        help="Rerun the simulations which did not complete succesfully in the previous run. Default=%default")
parser.add_option("-l","--logfile",action="store",type="string",dest="logfile_name",default="logfile.txt",help="Name of Log file. Default=%default")
(options,args)=parser.parse_args()
if(len(args)!=2):
    parser.error("Should specify 2 arguments -- binary file and output directory")

run_erring_sims=options.run_erring_sims
logfile_name=options.logfile_name
(binary,output_directory)=args

print "run_erring_sims =",run_erring_sims
print "binary=",binary
print "output_directory =",output_directory
print "logfile_name =",logfile_name

examples_sod={
        1:"Sod shock tube",
        4:"Lax's shock tube",
        5:"Strong shock tube",
        6:"Two symmetric rarefaction waves",
        7:"Mach 3 shock test",
        8:"High mach flow test",
        9:"Two shocks",
        10:"Interaction of blast waves (Bang Bang)"}

resolution_values=[401,801,1601]
#resolution_values=[401,801,1601,3201,6401,12801]
#resolution_values=[10,20,40]

eno_scheme_values={
        1:"",
        2:"density_weighted",
        3:"velocity_weighted"}

eno_order_values=[2]

rk_order_values=[3]

cfl_values=[.5]

timesplit_values=[True,False]

example_kind_values=["sod","smoothflow"]

test_number_values={
        "sod":examples_sod,
        "smoothflow":{1:""}}

if(not os.path.exists(output_directory)):
    try:
        os.makedirs(output_directory,0755)
        print "Directory %s created succesfully!"%output_directory
    except:
        print "Can't create directory %s"%output_directory
        sys.exit(0)

if(os.path.exists("%s/%s"%(output_directory,logfile_name))):
    old_log_file_format=output_directory+"/logfile_old_%02d.txt"
    last_index=Find_Last_Index_Of_File_Sequence(old_log_file_format,range(0,100))
    os.system("mv %s/%s %s/logfile_old_%02d.txt"%(output_directory,logfile_name,output_directory,last_index+1))
log_file=open("%s/%s"%(output_directory,logfile_name),"w")

Print_And_Log(log_file,"run_erring_sims = %s\n"%run_erring_sims)
Print_And_Log(log_file,"output_directory = %s\n"%output_directory)
Print_And_Log(log_file,"logfile_name = %s\n"%logfile_name)

for example_kind in example_kind_values:
    test_number_dictionary=test_number_values[example_kind]
    test_numbers=test_number_dictionary.keys()
    test_numbers.sort()
    for test_number in test_numbers:
        test_name=test_number_dictionary[test_number]
        log_message="\n\nRunning %s tests for test number %d (%s):\n\n"%(example_kind,test_number,test_name)
        Print_And_Log(log_file,log_message)
        for timesplit in timesplit_values:
            if(timesplit): implicit_rk_values=[True,False]
            else: implicit_rk_values=[False]
            for implicit_rk in implicit_rk_values:
                for cfl in cfl_values:
                    for rk_order in rk_order_values:
                        for eno_order in eno_order_values:
                            eno_schemes=eno_scheme_values.keys()
                            eno_schemes.sort()
                            for eno_scheme in eno_schemes:
                                eno_scheme_name=eno_scheme_values[eno_scheme]
                                for resolution in resolution_values:
                                    cfl_sound_speed_multiple=0
                                    (command_string,simulation_directory)=Run_Euler_1d(binary,output_directory,resolution,eno_scheme,eno_order,rk_order,cfl,cfl_sound_speed_multiple,implicit_rk,timesplit,example_kind,test_number)
                                    Run_Simulation(command_string,simulation_directory,run_erring_sims,log_file)


# Extra runs for smoothflow
example_kind_values=["smoothflow"]
for example_kind in example_kind_values:
    test_number_dictionary=test_number_values[example_kind]
    test_numbers=test_number_dictionary.keys()
    test_numbers.sort()
    for test_number in test_numbers:
        test_name=test_number_dictionary[test_number]
        log_message="\n\nRunning %s tests for test number %d (%s):\n\n"%(example_kind,test_number,test_name)
        Print_And_Log(log_file,log_message)
        for timesplit in timesplit_values:
            if(timesplit): implicit_rk_values=[True,False]
            else: implicit_rk_values=[False]
            for implicit_rk in implicit_rk_values:
                for cfl in cfl_values:
                    for rk_order in rk_order_values:
                        for eno_order in eno_order_values:
                            eno_schemes=eno_scheme_values.keys()
                            eno_schemes.sort()
                            for eno_scheme in eno_schemes:
                                eno_scheme_name=eno_scheme_values[eno_scheme]

                                # Sims at different resolution with cfl_sound_speed_multiple=3/cfl (i.e. effective sound speed cfl of 3)
                                resolution_values=[200,400,800,1600,3200]
                                if(timesplit): cfl_sound_speed_multiple=3./cfl
                                else: cfl_sound_speed_multiple=0
                                for resolution in resolution_values:
                                    (command_string,simulation_directory)=Run_Euler_1d(binary,output_directory,resolution,eno_scheme,eno_order,rk_order,cfl,cfl_sound_speed_multiple,implicit_rk,timesplit,example_kind,test_number)
                                    Run_Simulation(command_string,simulation_directory,run_erring_sims,log_file)

                                # Sims with scaling resolution and cfl_sound_speed_multiple by the same amount
                                resolution_base=3200
                                if(timesplit): cfl_sound_speed_multiple_base=3./cfl
                                else: cfl_sound_speed_multiple_base=0
                                scaling_values=[2,10,100]
                                for scale in scaling_values:
                                    resolution=resolution_base*scale
                                    cfl_sound_speed_multiple=cfl_sound_speed_multiple_base*scale
                                    (command_string,simulation_directory)=Run_Euler_1d(binary,output_directory,resolution,eno_scheme,eno_order,rk_order,cfl,cfl_sound_speed_multiple,implicit_rk,timesplit,example_kind,test_number)
                                    Run_Simulation(command_string,simulation_directory,run_erring_sims,log_file)

# Close the log file
log_file.close()
