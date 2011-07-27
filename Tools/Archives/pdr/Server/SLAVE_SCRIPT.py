# file_server_input = "pressure.stanford.edu:/var/data/armadillo/input"
# file_server_common = "pressure.stanford.edu:/var/data/armadillo/common"
# file_server_output = "pressure.stanford.edu:/var/data/armadillo/output"
# server_host = "pressure.stanford.edu"
# server_directory = "/usr/local/pdr"
# command_line_arguments = "-foo -bar"
import os,sys
import stat
import signal
import traceback
import shutil
import time
from random import randint

# unpack server supplied data
binary,platform,jobid,frame=sys.argv
render_pid=-1
killed=0

scp_command="scp -o StrictHostKeyChecking=no -qB"
ssh_command="ssh -o StrictHostKeyChecking=no -xnq"

def Handle_Kill_Signal(signum,frame):
    global render_pid
    global killed
    try:
        if render_pid>0:
            os.kill(render_pid,signal.SIGKILL)
    except:
        pass
    Log("Killed")
    Emit_Fail_Code()
    Cleanup()
    killed=1

def SSH_Command_With_Retry(command):
    Run_Command_With_Retry("%s %s \"%s\""%(ssh_command,server_host,command))

def SSH_Command(command):
    os.system("%s %s \"%s\""%(ssh_command,server_host,command))

def Run_Command_With_Retry(command):
    for i in range(100):
        try:
            status=os.system(command)
            if status==0:
                return
            print "Retrying command"
        except:
            print "Command Failed, Retrying"
        time.sleep(randint(3,10))
    raise NameError,"Failed after 100 retrys"

def Emit_Fail_Code():
    try:
        Report_Status_To_Server(1)
        copy_command=("%s Output/stdout.%s %s/"%(scp_command,frame,file_server_output))
        os.system(copy_command)
    except:
        pass # still want to send essential command
    command=("touch %s/Commands/Failed/%s-%s"%(server_directory,jobid,frame))
    SSH_Command_With_Retry(command)
    
def Emit_Success_Command():
    Report_Status_To_Server()
    command=("touch %s/Commands/Succeed/%s-%s"%(server_directory,jobid,frame))
    SSH_Command_With_Retry(command)

def Cleanup():
    shutil.rmtree("Input")
    shutil.rmtree("Common")
#    shutil.rmtree("Output")

def Make_Directories():
    map(os.mkdir,["Input","Output"])
    os.symlink("Input","Common") # so things that are groups of files some frame dependent and non-frame dependent can work with with the same prefix

def Setup_IO_Redirect():
    # Redirect output
    err_log=out_log=file(("Output/stdout.%s"%frame),'a+')
    dev_null=file('/dev/null','r')
    os.dup2(out_log.fileno(),sys.stdout.fileno())
    os.dup2(err_log.fileno(),sys.stderr.fileno())
    os.dup2(dev_null.fileno(),sys.stdin.fileno())

def Copy_Input_From_Server():
    file_server_input_path=("\"%s/*.%s\""%(file_server_input,frame))
    file_server_input_path_gz=("\"%s/*.%s.gz\""%(file_server_input,frame))
    file_server_common_path=("\"%s/*\""%(file_server_common))
    os.system("%s %s Input/"%(scp_command,file_server_input_path))
    os.system("%s %s Input/"%(scp_command,file_server_input_path_gz))
    os.system("%s %s Common/"%(scp_command,file_server_common_path))

def Copy_Output_To_Server():
    copy_command=("%s Output/* %s/"%(scp_command,file_server_output))
    #print copy_command
    if os.system(copy_command)!=0:
        time.sleep(randint(3,10))
        if os.system(copy_command)!=0:
            time.sleep(randint(3,10))
            if os.system(copy_command)!=0:
                raise NameError,"Failed to copy output to server server"

def Start_Child_Process_And_Wait():
    global render_pid
    global killed
    output_cycle_count=0
    render_pid=os.fork()
    if render_pid<=0: # Child Process
        os.nice(19) # reduce priority to minimum
        Run_Executable()
    else: # Server Process
        Report_Pid_To_Server()
        # Wait for child to die
        Log("Rendering (pid %d)"%render_pid)
        while 1:
            finished_pid,exitcode=os.waitpid(render_pid,os.WNOHANG)
            if finished_pid:
                if exitcode != 0: raise NameError,"Got non-zero error code"
                render_pid=-1
                break
            time.sleep(1)
            output_cycle_count+=1
            if killed==1:
                sys.exit(1)
            if output_cycle_count>30:
                output_cycle_count=0
                Report_Status_To_Server()

def Run_Executable(): # note this is done in child process
    if platform=="Linux":
        executable="./Linux_Executable"
    elif platform=="Windows":
        executable="./Windows_Executable.exe"
    else:
        raise NameError,"Invalid platform"
    # build command line
    args=[executable]
    command_line_arguments_substituted=command_line_arguments.replace("<frame>",frame) # substitute frame in
    args.extend(command_line_arguments_substituted.split(" "))
    #print args
    # Try to spawn the render task
    code=0
    try:
        #print "I am chmoding"
        os.chmod(executable,stat.S_IRUSR|stat.S_IWUSR|stat.S_IXUSR)
        #print "I am execing..."
        os.execvp(executable,args)
    except:
        # failed to execute
        code=1
    #print "I am exiting with code %d"%code
    Log("Executable exited with error code %d"%code)
    sys.exit(code)

def Report_Pid_To_Server():
    global render_pid
    my_pid=os.getpid()
    try:
        command=("echo %s %s > %s/Commands/Slave_Report_Data/%s-%s.pid"%(my_pid,render_pid,server_directory,jobid,frame))
        SSH_Command_With_Retry(command)        
        command=("touch %s/Commands/Slave_Report/%s-%s"%(server_directory,jobid,frame))
        SSH_Command_With_Retry(command)
    except:
        raise NameError,"Failed to give pid to server"

def Report_Status_To_Server(essential=0):
    try:
        sys.stdout.flush()
        copy_command=("%s Output/stdout.%s %s:%s/Commands/Slave_Report_Data/%s-%s.out"%(scp_command,frame,server_host,server_directory,jobid,frame))
        if essential:
            Run_Command_With_Retry(copy_command)
        else:
            os.system(copy_command)
        command=("touch %s/Commands/Slave_Report/%s-%s"%(server_directory,jobid,frame))
        SSH_Command(command)
    except:
        pass

def Log(text):
    print "\n<SLAVE>%s</SLAVE>\n"%text


if __name__=="__main__":
    signal.signal(signal.SIGTERM,Handle_Kill_Signal)
    try:
        Report_Pid_To_Server()
        Make_Directories()
        Setup_IO_Redirect()
        Log("Copying Input");Report_Status_To_Server()
        Copy_Input_From_Server()
        Log("Spawning Process")
        Start_Child_Process_And_Wait()
        print Log("Copying Output");Report_Status_To_Server()
        Copy_Output_To_Server()
        print Log("Success")
        Emit_Success_Command()
    except:
        if killed==0:
            Log("Failure render_pid=%d"%render_pid)
            if render_pid!=-1:
                Log("Killing %d"%render_pid)
                try:
                    os.kill(render_pid,signal.SIGKILL)
                except:pass
            traceback.print_exc(file=sys.stdout)
            Emit_Fail_Code()
    # whatever happens cleanup
    Cleanup()
