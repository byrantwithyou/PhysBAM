#!/usr/bin/python
from optparse import OptionParser
import re
import sys
import physbam

usage="usage: %prog [-m mtn_file] [-r triangulated_surface] [-t tetrahedralized_volume] scale_factor";
parser=OptionParser(usage);
parser.add_option('-m','--mtn',default="",type='string',help='motion file');
parser.add_option('-r','--tri',default="",type='string',help='triangulated surface');
parser.add_option('-t','--tet',default="",type='string',help='tetrahedralized volume');
(options,args)=parser.parse_args();
if len(args)==0 or len(args)>2:
    parser.error("incorrect number of arguments");
    sys.exit(1);

scale_factor=float(args[0]);
if options.mtn!="":
    options.mtn=re.sub(".gz","",options.mtn);
    new_file=re.sub("\.mtn","_scaled_%d.mtn"%scale_factor,options.mtn);
    motion_data=physbam.BODY_MOTION_SEQUENCE_f();
    physbam.Read_From_File("float",options.mtn,motion_data);
    motion_data.Rescale(scale_factor);
    physbam.Write_To_File("float",new_file,motion_data);
    print "Wrote %s"%new_file;

if options.tri!="":
    options.tri=re.sub(".gz","",options.tri);
    new_file=re.sub("\.tri","_scaled_%d.tri"%scale_factor,options.tri);
    tri_surface=physbam.TRIANGULATED_SURFACE_f.Create();
    physbam.Read_From_File("float",options.tri,tri_surface);
    tri_surface.Rescale(scale_factor);
    physbam.Write_To_File("float",new_file,tri_surface);
    print "Wrote %s"%new_file;

if options.tet!="":
    options.tet=re.sub(".gz","",options.tet);
    new_file=re.sub("\.tet","_scaled_%d.tet"%scale_factor,options.tet);
    tet_surface=physbam.TETRAHEDRALIZED_VOLUME_f.Create();
    physbam.Read_From_File("float",options.tet,tet_surface);
    tet_surface.Rescale(scale_factor);
    physbam.Write_To_File("float",new_file,tet_surface);
    print "Wrote %s"%new_file;
