import threading
import os
import time
import PDR_UTILITY


class PDR_STARTER(threading.Thread):
    def __init__(self,server,job,frame,slave):
        threading.Thread.__init__(self)
        self.server=server
        self.job=job
        self.frame=frame
        self.slave=slave
        
    def run(self):
        platform,hostname,instance=self.slave.split("__")

        remote_directory=("/local/%s-%s-%s/"%(self.server.my_hostname,self.job,self.frame))
        job_directory=os.path.join(self.server.jobs_directory,self.job)
        # Record start time
        start_filename=os.path.join(job_directory,"Status",self.frame+".start")
        time_str=str(int(time.time()))
        #print "start='%s' time_str='%s'" % (start_filename,time_str)
        open(start_filename,"w").write(time_str)

        try:
            PDR_UTILITY.Remote_Command(hostname,"mkdir -p %s"%remote_directory)
            PDR_UTILITY.Remote_Copy(hostname,os.path.join(job_directory,"Slave_Script"),remote_directory)
            if platform=="Windows":
                PDR_UTILITY.Remote_Copy(hostname,os.path.join(job_directory,"Windows_Executable.exe"),remote_directory)
            elif platform=="Linux":
                PDR_UTILITY.Remote_Copy(hostname,os.path.join(job_directory,"Linux_Executable"),remote_directory)
            PDR_UTILITY.Remote_Command(hostname,("cd %s;nohup python ./Slave_Script %s %s %s &"%(remote_directory,platform,self.job,self.frame)))
        except:
            os.system("touch %s"%os.path.join(self.server.server_root,"Commands","Failed","%s-%s"%(self.job,self.frame)))
            os.system("touch %s"%os.path.join(self.server.server_root,"Commands","Kill","%s-%s"%(self.job,self.frame)))
            
