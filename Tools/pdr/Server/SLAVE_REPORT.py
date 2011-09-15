import os
import PDR_UTILITY
from RENDER_COMMAND import RENDER_COMMAND

class SLAVE_REPORT(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Slave_Report",server)
        self.data_path=os.path.join(self.server.server_root,"Commands","Slave_Report_Data")

    def Process(self,id):
        job,frame=id.split("-")
        job_directory=os.path.join(self.server.jobs_directory,job)
        status_directory=os.path.join(job_directory,"Status")
        pid_file_in=os.path.join(self.data_path,("%s-%s.pid"%(job,frame)))
        out_file_in=os.path.join(self.data_path,("%s-%s.out"%(job,frame)))
        pid_file_out=os.path.join(status_directory,("%s.pid"%frame))
        out_file_out=os.path.join(status_directory,("%s.out"%frame))
        try: os.rename(pid_file_in,pid_file_out)
        except: pass
        try: os.rename(out_file_in,out_file_out)
        except: pass

    def Initialize_Directory_Structure(self,server):
        RENDER_COMMAND.Initialize_Directory_Structure(self,server)
        PDR_UTILITY.Make_Directory_If_Does_Not_Exist(os.path.join(self.data_path))
