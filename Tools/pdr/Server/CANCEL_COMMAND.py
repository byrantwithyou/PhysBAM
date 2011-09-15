import os
from RENDER_COMMAND import RENDER_COMMAND

class CANCEL_COMMAND(RENDER_COMMAND):
    def __init__(self,server):
        RENDER_COMMAND.__init__(self,"Cancel",server)

    def Process(self,id):
        jobid,frame=id.split("-")
        pending_file=os.path.join(self.server.jobs_directory,str(jobid),"Pending_Frames",str(frame))
        print "Removing in progress file '%s'" % (pending_file)
        os.remove(pending_file)


