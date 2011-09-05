#!/usr/bin/python
import sys
import os
import time

def Print_And_Log(log_file,log_message):
    """Prints the log_message and also writes into log_file"""

    print log_message
    log_file.write(log_message)
    log_file.flush()

def Loop_Read_Boolean_File_And_Run_Command(command_file_name,run_file_name,log_file_name,sleep_time):
    while(True):
        boolean_filehandle=open(run_file_name,'r')
        line=boolean_filehandle.readline().strip()
        print "line=",line
        boolean_filehandle.close()
        if(line=="run"):
            log_file=open(log_file_name,'a')
            command_filehandle=open(command_file_name,'r')
            command=command_filehandle.readline().strip()
            command_filehandle.close()
            Print_And_Log(log_file,"running command %s\n"%command)
            log_file.close()
            os.system(command)
        time.sleep(sleep_time)

def main():
    if(len(sys.argv)<5):
        print "Usage: %s <command_file_name> <run_file_name> <log_file_name> <sleep_time>"%sys.argv[0]
        sys.exit(1)

    command_file_name=sys.argv[1]
    run_file_name=sys.argv[2]
    log_file_name=sys.argv[3]
    sleep_time=float(sys.argv[4])

    print "******************************************"
    print "command_file_name=",command_file_name
    print "run_file_name=",run_file_name
    print "log_file_name=",log_file_name
    print "sleep_time=",sleep_time
    print "******************************************"

    Loop_Read_Boolean_File_And_Run_Command(command_file_name,run_file_name,log_file_name,sleep_time)


if __name__=="__main__":
    main()


#./convergence.py .03 tmp.dat tmp.png ../Standard_Tests/Test_1__Resolution_%d_%d_semiimplicit/ 50 10 1 2 3
