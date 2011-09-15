#!/usr/bin/python
from optparse import OptionParser
import sys
import os
import ProgressBar
import physbam

usage="usage: %prog [options]"
parser=OptionParser(usage);
parser.add_option("-m","--model_name",action="store",type="string",dest="model_name",default=None);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-d","--directory",action="store",type="string",dest="data_directory",default=None);
parser.add_option("-t","--proportional",action="store",type="string",dest="prop",default=None);
parser.add_option("-k","--keyframe",action="store",type="string",dest="key_frame",default=None);
parser.add_option("-l","--length",action="store",type="float",dest="default_length",default=.3);
parser.add_option("-n","--hair",action="store",type="int",dest="n",default=1000);
parser.add_option("-e","--segments",action="store",type="int",dest="n_segments",default=50);
parser.add_option("-v","--curve",action="store",type="float",dest="curve",default=0);
parser.add_option("-c","--curly",action="store_true",dest="curly",default=False);
(options,args)=parser.parse_args();
if len(args)!=0:
    parser.error("incorrect number of arguments");
    sys.exit(1);

if options.data_directory is None:
    parser.error("You must specify a data directory");
    sys.exit(1);

if options.sim_name is None:
    folders=options.data_directory.split("/");
    options.sim_name=folders.pop();
    options.data_directory="/".join(folders);
    print "WARNING: sim name is not specified, guessing %s"%options.sim_name;

if options.model_name is None:
    folders = options.data_directory.split("/");
    options.model_name=folders[len(folders)-1];
    print "WARNING: model name is not specified, guessing %s"%options.model_name;

cmd="python grow.py -s "+options.sim_name+" -l "+str(options.default_length)+" -n "+str(options.n)+" -e "+str(options.n_segments)+" -m "+options.model_name+" -v "+str(options.curve)+" -i";
if options.prop is not None:
    cmd+=" -t "+options.prop;
if options.curly:
    cmd+=" -c";
cmd+=" "+options.data_directory;
os.system(cmd);
os.system("python add_fixed_points.py -s "+options.sim_name+" -m "+options.model_name+" -d "+options.data_directory);
os.system("python graph.py -s "+options.sim_name+" -d "+options.data_directory);
if options.key_frame is not None:
    os.system("python generate_fixed_points.py -d -s "+options.sim_name+" -k "+options.key_frame+" -m "+options.model_name+" "+options.data_directory);
