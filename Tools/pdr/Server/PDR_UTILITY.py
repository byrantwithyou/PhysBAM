import os,sys

def Make_Directory_If_Does_Not_Exist(directory_path):
    if os.path.isdir(directory_path):return
    os.mkdir(directory_path)
    os.chmod(directory_path,0777);
    print "Made directory %s" % directory_path
        
def Remote_Command(hostname,command):
    cmd=("ssh -o StrictHostKeyChecking=no -o NumberOfPasswordPrompts=0 -xnq %s \"%s\""%(hostname,command))
    if os.system(cmd)!=0:
        raise NameError,"SSH command failed '%s'"%cmd

def Remote_Copy(hostname,source_file,destination_file):
    cmd=("scp -o StrictHostKeyChecking=no -o NumberOfPasswordPrompts=0 -qB %s %s:%s"%(source_file,hostname,destination_file))
    if os.system(cmd)!=0:
        raise NameError,"scp '%s' failed"%cmd
