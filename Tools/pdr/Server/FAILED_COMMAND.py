import os
from RENDER_COMMAND import RENDER_COMMAND

class FAILED_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Failed",server)

    def Process(self,id):
        jobid,frame=id.split("-")
        pending_file=os.path.join(self.server.jobs_directory,str(jobid),"In_Progress_Frames",str(frame))
        failed_file=os.path.join(self.server.jobs_directory,str(jobid),"Failed_Frames",str(frame))
        os.rename(pending_file,failed_file)

        
