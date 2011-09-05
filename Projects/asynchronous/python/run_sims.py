#!/usr/bin/python

from subprocess import *
import sys
import os
from optparse import OptionParser

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

def Run_Process(command_string,log_file):
    """Runs command_string and returns on the completion of the process"""

    command_string_split=command_string.split(" ");
    log_message="running command %s\n"%command_string
    Print_And_Log(log_file,log_message)

    process=Popen(command_string_split)
    process.wait()

    return process.returncode

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

    returncode=Run_Process(command_string,log_file)

    if(returncode==None):
        log_message="Error: process hasn't terminated\n\n"
    elif(returncode==0):
        log_message="Process completed succesfully\n\n"
        done_file=open("%s/done"%simulation_directory,"w")
        done_file.close()
    else:
        log_message="Process did not completed succesfully\n\n"
        error_file=open("%s/error"%simulation_directory,"w")
        error_file.close()
    Print_And_Log(log_file,log_message)


def Get_Asynchronous_Command(binary,output_directory,test_number,model_number,number_of_spheres,asynchronous_mode,overdamping_fraction,projection_rigidity,make_it_bad,last_frame):
    simulation_directory="%s/Test_%d_model_%d_%s_overdamping_%f"%(output_directory,test_number,model_number,asynchronous_mode,overdamping_fraction)
    if(model_number==1): simulation_directory+="_spheres_%d"%number_of_spheres
    if(asynchronous_mode=="async"): simulation_directory+="_prigid_%f"%projection_rigidity
    if(make_it_bad): simulation_directory+="_makeitbad_%f"%make_it_bad

    command_string=binary
    command_string+=" -model_num %d"%model_number
    if(asynchronous_mode=="async"):
        command_string+=" -async"
    elif(asynchronous_mode=="explicit"): command_string+=" -dt 0"
    elif(asynchronous_mode=="implicit"): command_string+=" -fully_implicit"
    if(model_number==1): command_string+=" -number_of_spheres %d"%number_of_spheres
    command_string+=" -overdamping_fraction %f"%overdamping_fraction
    if(asynchronous_mode=="async"): command_string+=" -projection_rigidity %f"%projection_rigidity
    if(make_it_bad): command_string+=" -make_it_bad %f"%make_it_bad
    command_string+=(" -last_frame %d"%last_frame)
    command_string+=(" -o %s"%simulation_directory)
    command_string+=(" %d"%test_number)

    return (command_string,simulation_directory)

def main():
    usage="usage: %prog [options] <binary> <output_directory>"
    parser=OptionParser(usage)
    parser.add_option("--run_erring_sims",action="store_true",dest="run_erring_sims",default=False,
            help="Rerun the simulations which did not complete succesfully in the previous run. Default=%default")
    parser.add_option("-l","--logfile",action="store",type="string",dest="logfile_name",default="logfile.txt",help="Name of Log file. Default=%default")
    (options,args)=parser.parse_args()
    if(len(args)!=2):parser.error("Should specify 2 arguments -- binary file and output directory")

    run_erring_sims=options.run_erring_sims
    logfile_name=options.logfile_name
    (binary,output_directory)=args

    print "run_erring_sims =",run_erring_sims
    print "binary=",binary
    print "output_directory =",output_directory
    print "logfile_name =",logfile_name

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


    test_number_values=[24]
    model_number_values=[1,2]
    asynchronous_mode_values=["async","explicit","implicit"]
    overdamping_fraction_values=[1]
    projection_rigidity_values=[0,.25,.5,.75,1]
    make_it_bad_values=[0]

    for test_number in test_number_values:
        for model_number in model_number_values:
            if(model_number==1):
                number_of_spheres_values=[1,2]
                last_frame=300
            else:
                number_of_spheres_values=[0]
                last_frame=400
            for number_of_spheres in number_of_spheres_values:
                if(number_of_spheres==2): last_frame=200
                if(model_number==1 & number_of_spheres==1): make_it_bad_values=[0,.01,.001]
                else: make_it_bad_values=[0]
                for make_it_bad in make_it_bad_values:
                    for asynchronous_mode in asynchronous_mode_values:
                        if(asynchronous_mode=="async"): projection_rigidity_values_temp=projection_rigidity_values
                        else: projection_rigidity_values_temp=[0]
                        for overdamping_fraction in overdamping_fraction_values:
                            for projection_rigidity in projection_rigidity_values_temp:
                                (command_string,simulation_directory)=Get_Asynchronous_Command(binary,output_directory,test_number,model_number,number_of_spheres,asynchronous_mode,overdamping_fraction,projection_rigidity,make_it_bad,last_frame)
                                print "%s"%command_string
                                #Run_Simulation(command_string,simulation_directory,run_erring_sims,log_file)

if __name__=="__main__":
    main()
