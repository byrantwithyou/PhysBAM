import os
from RENDER_COMMAND import RENDER_COMMAND

class RETRY_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Retry",server)

    def Process(self,id):
        jobid,frame=id.split("-")
        failed_file=os.path.join(self.server.jobs_directory,str(jobid),"Failed_Frames",str(frame))
        pending_file=os.path.join(self.server.jobs_directory,str(jobid),"Pending_Frames",str(frame))
        print "moving from '%s' to '%s'" % (failed_file,pending_file)
        os.rename(failed_file,pending_file)


