import os
from RENDER_COMMAND import RENDER_COMMAND

class SUCCEEDED_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Succeed",server)

    def Process(self,id):
        jobid,frame=id.split("-")
        pending_file=os.path.join(self.server.jobs_directory,str(jobid),"In_Progress_Frames",str(frame))
        succeeded_file=os.path.join(self.server.jobs_directory,str(jobid),"Done_Frames",str(frame))
        print "moving from '%s' to '%s'" % (pending_file,succeeded_file)
        os.rename(pending_file,succeeded_file)

        
