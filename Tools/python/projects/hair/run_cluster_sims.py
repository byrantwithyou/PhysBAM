#!/usr/bin/python
from optparse import OptionParser
import os
import re
import time
import shutil
import sys

usage="usage: %prog [options]"
parser=OptionParser(usage);
parser.add_option("-m","--model_name",action="store",type="string",dest="model_name",default=None);
parser.add_option("-n","--mon_name",action="store",type="string",dest="mon_name",default=None);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-d","--directory",action="store",type="string",dest="data_directory",default=None);
parser.add_option("-p","--parameter",action="store",type="string",dest="params",default=None);
parser.add_option("-b","--binary",action="store",type="string",dest="binary",default=None);
parser.add_option("-u","--num_procs",action="store",type="int",dest="num_procs",default=1);
parser.add_option("-t","--test",action="store",type="int",dest="test_number",default=1);
parser.add_option("-o","--output_base",action="store",type="string",dest="output_base",default="/solver/vol3/hair3/sims");
parser.add_option("-r","--recursive",action="store_true",dest="recursive",default=False);
(options,args)=parser.parse_args();
if len(args)!=0:
    parser.error("incorrect number of arguments");
    sys.exit(1);

if options.data_directory is None:
    parser.error("You must specify a data directory");
    sys.exit(1);

if options.binary is None:
    parser.error("You must specify a binary");
    sys.exit(1);

if options.recursive:
    for dir in os.listdir(options.data_directory):
        if re.match(".*10k.*",dir) and os.path.isdir(options.data_directory+"/"+dir):
            cmd="python "+" ".join(sys.argv);
            cmd=re.sub("(-d\s+[^\s]+)","\\1/"+dir,cmd);
            cmd=re.sub("-r ","",cmd);
            os.system(cmd);
    sys.exit(0);

if options.sim_name is None:
    folders=options.data_directory.split("/");
    options.sim_name=folders.pop();
    options.data_directory="/".join(folders);
    print "WARNING: sim name is not specified, guessing %s"%options.sim_name;

if options.model_name is None:
    folders = options.data_directory.split("/");
    options.model_name=folders[len(folders)-1];
    print "WARNING: model name is not specified, guessing %s"%options.model_name;

if options.mon_name is None:
    options.mon_name=options.sim_name;

datestr="%04d%02d%02d-%02d%02d"%(time.localtime()[:5])



rundirs=[]
for file in os.listdir(options.data_directory+"/"+options.sim_name):
    if file==options.params or (options.params is None and file.startswith("parameters")):
        parameters=file.replace("parameters.","");
        sim_directory_base=os.path.join(options.output_base,datestr+"-"+options.mon_name)
        if not os.path.exists(sim_directory_base): # "/solver/vol3/hair3/sims/%s-%s"%(datestr,options.mon_name)):
            os.mkdir(sim_directory_base)
        rundir=os.path.join(sim_directory_base,parameters)
        print rundir
        rundirs.append(rundir)
        mon_name=options.mon_name+"-"+parameters;

        if os.path.exists(rundir):
            print "Path exists already try new sim name"
        else:
            os.mkdir(rundir)
            os.system("%s --copy %s"%(options.binary,rundir))
            shutil.copy("%s/Scripts/linux_sim_scripts/mpi_cluster_sim.py"%os.environ["PHYSBAM"],rundir)
            os.chdir(rundir)
            cmd=("python mpi_cluster_sim.py %d %s ./solids_3d -example HAIR_SIM_TESTS -d %s -modelname %s -hairsim %s -params %s %d "%(options.num_procs,mon_name,options.data_directory,options.model_name,options.sim_name,file,options.test_number))
            print cmd
            os.system(cmd)
            if options.num_procs==1: proc_fold=""
            else: proc_fold="1/"
            restart_script="""import sys,os

if len(sys.argv)>1:
    last_frame=int(sys.argv[1])
else:
    dirs=os.listdir("Hair_Sim_Tests")
    if len(dirs)>1:
        print "FATAL: Output Directory is ambiguous, please provide a restart frame";
    last_frame=int(open("Hair_Sim_Tests/%%s/%slast_frame"%%dirs[0]).read())
print "Starting from Frame %%d"%%last_frame
resp=raw_input("Do you wish to continue? [y/n]")
if resp != 'y': sys.exit(0)
os.system("python mpi_cluster_sim.py %d %s ./solids_3d -example HAIR_SIM_TESTS -d %s -modelname %s -hairsim %s -params %s -restart %%d %d"%%last_frame)"""%(proc_fold,options.num_procs,mon_name,options.data_directory,options.model_name,options.sim_name,file,options.test_number)
            try:
                restart_file=open("%s/restart.py"%rundir,"w");
                restart_file.write(restart_script);
                restart_file.close();
            except IOError:
                print "WARNING: Couldn't write restart script"
            

print "All Run directories"
print "\n".join(rundirs)

