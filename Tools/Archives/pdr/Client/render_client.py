#!/usr/bin/python
#
# This script assumes we have file system access to the render server

import os
import os.path
import sys
import shutil

class RENDER_CLIENT:
    def __init__(self,server_host="advection.stanford.edu",server_directory="/usr/local/pdr"):
        self.server_host=server_host
        self.server_directory=server_directory
        self.server_path="/n/"+self.server_host+server_directory
        self.command_path=os.path.join(self.server_path,"Commands")

    def dispatch(self,job_name,executable,command_line,frames,file_server,file_directory,memory_requirement=1,priority=2):
        old_umask=os.umask(0)
        username=os.environ['USER']
        job_directory=os.path.join(self.command_path,"Create_Data",username+"-"+job_name)
        os.mkdir(job_directory)
        # write job info file
        job_info_text="\n".join([job_name,username,str(priority),str(memory_requirement)," ".join(map(str,frames))])+"\n";
        open(os.path.join(job_directory,"Job_Info"),"w").write(job_info_text)
        # copy executable
        shutil.copy(executable,os.path.join(job_directory,"Linux_Executable"))
        # Do the script header
        open(os.path.join(job_directory,"Script_Header"),"w").writelines(
            ["file_server_input=%s:%s/Input\n" % (file_server,file_directory),
             "file_server_common=%s:%s/Common\n" % (file_server,file_directory),
             "file_server_output=%s:%s/Output\n" % (file_server,file_directory),
             "server_host=%s\n" % self.server_host,
             "server_directory=%s\n" % self.server_directory,
             "command_line_arguments=%s\n" % command_line])
        # tell server we're done setting up the info by touching the file
        open(os.path.join(self.command_path,"Create",username+"-"+job_name),"w")
        os.umask(old_umask)

    def retry(self,job_id,frames):
        old_umask=os.umask(0)
        for frame in frames: open(os.path.join(self.command_path,"Retry","%d-%d"%(job_id,frame)),"w")
        os.umask(old_umask)

    def add_slave(self,slave):
        old_umask=os.umask(0)
        open(os.path.join(self.command_path,"Add_Slave",slave),"w")
        os.umask(old_umask)

    def remove_slave(self,slave):
        old_umask=os.umask(0)
        open(os.path.join(self.command_path,"Remove_Slave",slave),"w")
        os.umask(old_umask)

    def add_chromiums_single(self):
        for i in range(1,33):
          self.add_slave("Linux__2__512__chromium" + str(i) + "__1")

    def add_spire_first(self):
        for i in range(1,17):
          self.add_slave("Linux__2__512__spire-" + str(i) + "__1")

    def add_spire_second(self):
        for i in range(1,17):
          self.add_slave("Linux__2__512__spire-" + str(i) + "__2")

    def remove_chromiums_single(self):
       for i in range(1,33):
          self.remove_slave("chromium" + str(i) + "__1")

    def remove_spire_first(self):
       for i in range(1,17):
          self.remove_slave("spire-" + str(i) + "__1")

    def remove_spire_second(self):
       for i in range(1,17):
          self.remove_slave("spire-" + str(i) + "__2")

    def framestr_to_framelist(self,framestr):
        x=framestr.strip()
        framehash={}
        for guy in x.split(","):
            guys=guy.split("-")
            if len(guys)==2:
                begin,end=guys
                begin,end=int(begin),int(end)
                for i in range(begin,end+1):
                    framehash[i]=1
            else:
                print "guy=%s"%guy
                framehash[int(guy)]=1
        frames=framehash.keys()
        frames.sort()
        return frames

def print_usage():
    print "Usage is: "
    print "   render_client.py submit <jobname> <executable> <commandline> <frame string> <file server> <file directory>"
    print "   render_client.py retry <jobid> <frame string>"
    print "   render_client.py add <machine1> <machine2> <machin3> ..."
    print "   render_client.py add_chromiums_single"
    print "   render_client.py add_spire_first"
    print "   render_client.py add_spire_second"
    print "   render_client.py remove <machine1> <machine2> <machine3> ..."
    print "   render_client.py remove_chromiums_single"
    print "   render_client.py remove_spire_first"
    print "   render_client.py remove_spire_second"
    sys.exit(1)

if __name__=="__main__":
    client=RENDER_CLIENT()

    if len(sys.argv) < 2: print_usage()
    if sys.argv[1]=="submit":
        jobname,executable,commandline,framestr,file_server,file_directory=sys.argv[2:]
        frames=client.framestr_to_framelist(framestr)
        client.dispatch(jobname,executable,commandline,frames,file_server,file_directory)
    elif sys.argv[1]=="retry":
        jobid,framestr=sys.argv[2:]
        frames=client.framestr_to_framelist(framestr)
        print "Retry with jobid=%d and frames=%s" % (int(jobid),repr(frames))
        client.retry(int(jobid),frames)
    elif sys.argv[1]=="add":
        for i in range(2,len(sys.argv)):
          print "Adding " + sys.argv[i]
          client.add_slave(sys.argv[i]) 
    elif sys.argv[1]=="add_chromiums_single":
        print "Adding all 32 chromium processors(single proc)"
        client.add_chromiums_single() 
    elif sys.argv[1]=="add_spire_single":
        print "Adding all first 16 spire processors(single proc)"
        client.add_spire_first() 
    elif sys.argv[1]=="add_spire_second":
        print "Adding all second 16 spire processors(single proc)"
        client.add_spire_second() 
    elif sys.argv[1]=="remove":
        for i in range(2,len(sys.argv)):
          print "Removing " + sys.argv[i]
          client.remove_slave(sys.argv[i]) 
    elif sys.argv[1]=="remove_chromiums_single":
        print "Removing all 32 chromium processors(single proc)"
        client.remove_chromiums_single() 
    elif sys.argv[1]=="remove_spire_first":
        print "Removing first 16 spire processors(single proc)"
        client.remove_spire_first() 
    elif sys.argv[1]=="remove_spire_second":
        print "Removing second 16 spire processors(single proc)"
        client.remove_spire_second() 
    else:
        print_usage()
