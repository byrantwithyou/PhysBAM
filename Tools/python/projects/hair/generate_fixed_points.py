from optparse import OptionParser
import physbam
import os
import sys
import ProgressBar

usage="usage: %prog [options] <data directory>"
parser=OptionParser(usage);
parser.add_option("-m","--model_name",action="store",type="string",dest="model_name",default=None);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-d","--deformable",action="store_true",dest="deformable",default=False);
parser.add_option("-k","--keyframe_motion",action="store",type="string",dest="keyframe_motion",default=None);
(options,args)=parser.parse_args();
if len(args)!=1:
    parser.error("incorrect number of arguments");
    sys.exit(1);

if options.sim_name is None:
    parser.error("You must specify a sim name to save to")
    sys.exit(1);

if options.model_name is None:
    folders = args[0].split("/");
    options.model_name=folders[len(folders)-1];
    print "WARNING: model name is not specified, guessing %s"%options.model_name;

fixed_points=physbam.LA_i()
fixed_tets=physbam.LA_i()
fixed_weights=physbam.LA_Vf3()
particles=physbam.SOLIDS_PARTICLE_Vf3()
tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create()
tri=physbam.TRIANGULATED_SURFACE_f.Create(particles)

physbam.Read_From_File("float",os.path.join(args[0],options.sim_name,"fixed_tets.gz"),fixed_tets)
physbam.Read_From_File("float",os.path.join(args[0],options.sim_name,"fixed_nodes.gz"),fixed_points)
physbam.Read_From_File("float",os.path.join(args[0],options.sim_name,"fixed_weights.gz"),fixed_weights)

last_frame=1
if options.deformable:
    motion_symlink_target=os.path.join(args[0],options.sim_name,"motion")
    if os.path.exists(motion_symlink_target) and os.path.islink(motion_symlink_target):
        os.unlink(motion_symlink_target)
    os.symlink(options.keyframe_motion,os.path.join(args[0],options.sim_name,"motion"))
    last_frame=int(open(os.path.join(args[0],options.sim_name,"motion","last_frame")).read())
else:
    tri_file=options.model_name;
    physbam.Read_From_File("float",os.path.join(args[0],options.sim_name,"particles"),particles)
    physbam.Read_From_File("float",os.path.join(args[0],"Rigid_Bodies",tri_file+".tri"),tri.mesh);

print "Calculating fixed points"
bar=ProgressBar.progressBar(0,last_frame,80)
for frame in range(0,last_frame):
    bar(frame)
    #print "%d/%d"%(frame,last_frame)
    points=physbam.LA_Vf3()
    if options.deformable:
        physbam.Read_From_File("float",os.path.join(args[0],options.sim_name,"motion","out.%d.tet"%frame),tet);
        for j in range(1,len(fixed_points)+1):
            tet_index=fixed_tets[j]
            cur_tet=tet.mesh.elements[tet_index];
            t=physbam.TETRAHEDRON_f(*map(lambda i: tet.particles.X[i],cur_tet))
            points.Append(t.Point_From_Barycentric_Coordinates(fixed_weights[j]));
    else:
        for j in range(1,len(fixed_points)+1):
            points.Append(tri.particles.X[fixed_points[j]]);
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/fixed_positions.%d"%frame,points);
print ""
