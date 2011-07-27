import os
from RENDER_COMMAND import RENDER_COMMAND
import PDR_UTILITY

class KILL_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Kill",server)

    def Process(self,id):
        jobid,frame=id.split("-")
        status_file=os.path.join(self.server.jobs_directory,jobid,"Status",frame+".pid")
        in_progress_file=in_progress_directory=os.path.join(self.server.jobs_directory,jobid,"In_Progress_Frames",frame)
        slave_pid,render_pid=open(status_file).read().split(" ")
        slave=open(in_progress_file).read().replace("\n","")
        platform,host,instance=slave.split("__")
        server=os.path.join(self.server.jobs_directory,jobid,"Status",frame+".pid")
        PDR_UTILITY.Remote_Command(host,"kill %s"%slave_pid)

