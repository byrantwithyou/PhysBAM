import os
import shutil
import sys
import traceback
import PDR_UTILITY
from RENDER_COMMAND import RENDER_COMMAND

class CREATE_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Create",server)
        self.data_path=os.path.join(self.server.server_root,"Commands","Create_Data")

    def Process(self,submit_id):
        filename=os.path.basename(submit_id)
        job_submit_directory=os.path.join(self.data_path,filename)
        infofilename=os.path.join(job_submit_directory,"Job_Info")
        # try to read information
        try:
            job_directory="-1"
            infolines=[s.replace("\n","") for s in open(infofilename).readlines()]
            print infolines
            jobname,username,priority,memory_requirement,frames=infolines[:5]
            ## Construct the job
            # Get the job id and make the directory
            jobid=self.server.New_Job_ID()
            job_directory=os.path.join(self.server.jobs_directory,str(jobid))
            os.mkdir(job_directory)
            [os.mkdir(os.path.join(job_directory,dirname)) for dirname in ["Pending_Frames","In_Progress_Frames","Done_Frames","Failed_Frames","Status"]]
            # copy executables and script
            copied={}
            for i in ["Windows_Executable.exe","Linux_Executable","Job_Info"]:
                try:
                    shutil.copy(os.path.join(job_submit_directory,i),os.path.join(job_directory,i));copied[i]=1
                except IOError:pass
            fp=open(os.path.join(job_directory,"Slave_Script"),"w")
            fp.write("#!/usr/bin/python\n")
            for line in open(os.path.join(job_submit_directory,"Script_Header"),"r").readlines():
                key,val=line.replace("\n","").split("=")
                fp.write("%s=\"%s\"\n"%(key,val))
            fp.write(open("bin/SLAVE_SCRIPT.py","r").read())
            fp.close()
            # Check to make sure everything we need is there
            if not copied.has_key("Windows_Executable.exe") and not copied.has_key("Linux_Executable"):
                raise NameError, "Didn't get binary"
            # add all frames to pending directory
            pending_frame_dir=os.path.join(job_directory,"Pending_Frames")
            frame_list=filter(lambda x:x!="",frames.split(" "))
            for i in frame_list: open(os.path.join(pending_frame_dir,i),"w")
        except:
            print "Failed to construct job '%s':"%submit_id
            traceback.print_exc(file=sys.stdout)
            if job_directory!="-1":shutil.rmtree(job_directory)
        # get rid of submit directory
        try:
            shutil.rmtree(job_submit_directory)
        except:
            print "Failed to delete submit directory '%s'"%submit_id
        

    def Initialize_Directory_Structure(self,server):
        RENDER_COMMAND.Initialize_Directory_Structure(self,server)
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.data_path))
        
