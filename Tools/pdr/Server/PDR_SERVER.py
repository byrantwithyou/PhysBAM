#!/usr/bin/python
import time
import sys
import os
import traceback
import PDR_UTILITY

from CREATE_COMMAND import CREATE_COMMAND 
from ADD_SLAVE_COMMAND import ADD_SLAVE_COMMAND
from REMOVE_SLAVE_COMMAND import REMOVE_SLAVE_COMMAND
from CANCEL_COMMAND import CANCEL_COMMAND
from FAILED_COMMAND import FAILED_COMMAND
from SUCCEEDED_COMMAND import SUCCEEDED_COMMAND
from RETRY_COMMAND import RETRY_COMMAND
from SLAVE_REPORT import SLAVE_REPORT
from KILL_COMMAND import KILL_COMMAND
from PDR_STARTER import PDR_STARTER
from PDR_FULL_SCHEDULER import PDR_FULL_SCHEDULER
from PDR_FCFS_SCHEDULER import PDR_FCFS_SCHEDULER

class PDR_SERVER:
    def __init__(self,server_root):
        # get my hostname
        self.my_hostname=os.popen("hostname").read().replace("\n","")
        # define directories
        self.server_root=server_root
        self.state_directory=os.path.join(self.server_root,"State")
        self.jobs_directory=os.path.join(self.server_root,"State","Jobs")
        self.slaves_directory=os.path.join(self.server_root,"State","Slaves")
        # load commands
        self.Instantiate_Commands()
        # make sure our directory structure is here
        self.Initialize_Directory_Structure()
        # find where to start job ids (compute max of remaining jobs and add 1)
        try:
            self.next_jobid=1+reduce(max,map(int,os.listdir(self.jobs_directory)))
        except:
            self.next_jobid=70
        # Open logs
        self.command_log=open("command_log","a+")
        # Don't wait too long
        self.wait_between_cycles=20
        self.job_start_timeout=60*5
        self.jobs_start_throttle=2
        # Load scheduler
        self.scheduler=PDR_FULL_SCHEDULER()

    def Log(self,log,text):
        log.write("%s %s\n"%(time.ctime(),text))
        log.flush()

    def Instantiate_Commands(self):
        self.commands=[]
        self.commands.append(CREATE_COMMAND(self))
        self.commands.append(RETRY_COMMAND(self))
        self.commands.append(CANCEL_COMMAND(self))
        self.commands.append(ADD_SLAVE_COMMAND(self))
        self.commands.append(REMOVE_SLAVE_COMMAND(self))
        self.commands.append(SUCCEEDED_COMMAND(self))
        self.commands.append(FAILED_COMMAND(self))
        self.commands.append(SLAVE_REPORT(self))
        self.commands.append(KILL_COMMAND(self))

    def Process(self):
        while 1:
            print "Running..."
            for i in self.commands:
                i.Handle_Pending_Commands()
            self.Manage_Execution()
            time.sleep(self.wait_between_cycles)

    def Initialize_Directory_Structure(self):
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.server_root,"State"))
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.server_root,"State","Jobs"))
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.server_root,"State","Slaves"))
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.server_root,"Commands"))
        # setup each command's directories
        for i in self.commands:
            i.Initialize_Directory_Structure(self)
                    
    def New_Job_ID(self):
        jobid=self.next_jobid
        self.next_jobid+=1
        return jobid

    def Get_Slave_And_Job_Information(self):
        pending_jobs=[]
        used_slaves={}   # the slave and what job it is on
        unknown_frames=[]
        for job in os.listdir(self.jobs_directory):
            pending_directory=os.path.join(self.jobs_directory,job,"Pending_Frames")
            # Figure out used_slaves and unknown_frames by going through in_progress_frames
            in_progress_directory=os.path.join(self.jobs_directory,job,"In_Progress_Frames")
            try:
                files=os.listdir(os.path.join(self.jobs_directory,job))
                for frame in os.listdir(in_progress_directory):
                    try:
                        used_slaves[open(os.path.join(in_progress_directory,frame)).readline().replace("\n","")]=(job,frame)
                        if not os.path.exists(os.path.join(self.jobs_directory,job,"Status","%s.out"%frame)):
                            unknown_frames.append( (job,frame) )
                    except:
                        print "Failed to get slave on job %s for in progress frame  %s" % (job,frame)
            except:
                print "Failed to get in progress frames for job " + job
            # Get list of unprocessed frames
            try:
                fp=open(os.path.join(self.jobs_directory,job,"Job_Info"))
                fp.readline();fp.readline();priority=int(fp.readline().replace("\n",""));memory_requirement=int(fp.readline().replace("\n",""))
                pending_jobs.append( (int(job),int(priority),int(memory_requirement),map(int,os.listdir(pending_directory))) )
            except:
                print "Failed to get pending frames for job " + job
                traceback.print_exc(file=sys.stdout)
                
                
        return pending_jobs,used_slaves,unknown_frames

    def Manage_Execution(self):
        jobs_started=0
        # get list of frames and list of used slaves
        pending_jobs,used_slaves,unknown_frames=self.Get_Slave_And_Job_Information()
        # timeout
        self.Timeout(unknown_frames)
        # construct and write free slave list by subtracting used slaves
        potential_slave_list=os.listdir(self.slaves_directory)
        free_slaves_strings=[]
        for i in potential_slave_list:
            platform,priority,memory,hostname,num=i.split("__")
            if not used_slaves.has_key("__".join([platform,hostname,num])):
                free_slaves_strings.append(i)
        free_slaves=[]
        for slave in free_slaves_strings:
            platform,priority,memory,host,instance=slave.split("__")
            free_slaves.append((platform+"__"+host+"__"+instance,int(priority),int(memory)))
        open(os.path.join(self.state_directory,"Free_Slaves"),"w").write("\n".join(map(repr,free_slaves))+"\n")
        print "Got %d pending jobs and %d unused slaves"%(len(pending_jobs),len(free_slaves))
        # Spawn the external scheduler
        matches=self.scheduler.Schedule(pending_jobs,free_slaves,self.jobs_start_throttle)
        # Run pairs
        for job_int,frame_int,slave in matches:
            frame=str(frame_int)
            job=str(job_int)
            os.unlink(os.path.join(self.jobs_directory,job,"Pending_Frames",frame))
            open(os.path.join(self.jobs_directory,job,"In_Progress_Frames",frame),"w").write(slave + "\n")
            try:
                os.unlink(os.path.join(self.jobs_directory,job,"Status",frame+".pid"))
                os.unlink(os.path.join(self.jobs_directory,job,"Status",frame+".out"))
            except:
                pass
            starter_thread=PDR_STARTER(self,job,frame,slave)
            starter_thread.start()

    def Timeout(self,unknown_frames):
        for job,frame in unknown_frames:
            job_directory=os.path.join(self.jobs_directory,job)
            try:
                start_time=open(os.path.join(self.jobs_directory,job,"Status",frame+".start")).read()
                if int(time.time())-int(start_time)>self.job_start_timeout:
                    os.system("touch %s"%os.path.join(self.server_root,"Commands","Failed","%s-%s"%(job,frame)))
            except:
                print "Failed to get start time"


if __name__=="__main__":
    server=PDR_SERVER(os.path.abspath(os.getcwd()))   
    server.Process()
