#!/usr/bin/python
from subprocess import *

def Run_Simulation(command_string,log_fd):
    command_string_split=command_string.split(" ");
    #print "running command %s\n"%command_string

    # Run the process
    process=Popen(command_string_split,stdout=log_fd,stderr=log_fd);
    process.wait()
    if(process.returncode==None):
        #print "Error: process hasn't terminated\n\n"
        return -1
    elif(process.returncode==0):
        #print "Process completed succesfully\n\n"
        return 0
    else:
        #print "Process did not completed succesfully\n\n"
        return 1

def Create_Command_String(binary,test_number,resolution,eno_order,rk_order,cfl):
    command_string=binary
    command_string+=(" -resolution %d"%resolution)
    command_string+=(" -eno_order %d"%eno_order)
    command_string+=(" -rk_order %d"%rk_order)
    command_string+=(" -cfl %f"%cfl)
    command_string+=(" %d"%test_number)

    return command_string
    
def Binary_Search(binary,test_number,resolution,eno_order,rk_order):
    log_fd=open("binary_search_log_test_%d_resolution_%d_eno_order_%d_rk_order_%d_.txt"%(test_number,resolution,eno_order,rk_order),"w");
    lower=.5;
    upper=1;
    epsilon=.01

    last_failed=1
    while((upper-lower)>epsilon):
        middle=(lower+upper)*.5;

        command_string=Create_Command_String(binary,test_number,resolution,eno_order,rk_order,middle)
        if(Run_Simulation(command_string,log_fd)==0):
            lower=middle
        else:
            upper=middle
            last_failed=middle

    return last_failed

def main():
    binary="./euler_fsi_1d_nocona_debug -sod -last_frame 10"
    resolution=200
    test_number=1

    eno_order_values=[1,2,3]
    rk_order_values=[1,2,3]
    resolution_values=[200,400]

    for eno_order in eno_order_values:
        for rk_order in rk_order_values:
            print
            for resolution in resolution_values:
                last_failed=Binary_Search(binary,test_number,resolution,eno_order,rk_order);
                print "test_number=%d, resolution=%d, eno_order=%d, rk_order=%d, last failed cfl=%f"%(test_number,resolution,eno_order,rk_order,last_failed)

if __name__=="__main__":
    main()
