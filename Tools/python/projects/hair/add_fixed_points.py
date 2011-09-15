#!/usr/bin/python
from optparse import OptionParser
import sys
import os
import ProgressBar
import physbam

def subset(X,nodes):
    return map(lambda x: X[x],nodes)

usage="usage: %prog [options]"
parser=OptionParser(usage)
parser.add_option("-m","--model_name",action="store",type="string",dest="model_name",default=None)
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None)
parser.add_option("-d","--directory",action="store",type="string",dest="data_directory",default=None)
(options,args)=parser.parse_args()
if len(args)!=0:
    parser.error("incorrect number of arguments")
    sys.exit(1)

if options.data_directory is None:
    parser.error("You must specify a data directory")
    sys.exit(1)

if options.sim_name is None:
    folders=options.data_directory.split("/")
    options.sim_name=folders.pop()
    options.data_directory="/".join(folders)
    print "WARNING: sim name is not specified, guessing %s"%options.sim_name

if options.model_name is None:
    folders = options.data_directory.split("/")
    options.model_name=folders[len(folders)-1]
    print "WARNING: model name is not specified, guessing %s"%options.model_name

particles=physbam.SOLIDS_PARTICLE_Vf3()
tet_particles=physbam.SOLIDS_PARTICLE_Vf3()
edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(tet_particles)
fixed_nodes=physbam.LA_i()
fixed_nodes2=physbam.LA_i()
fixed_tets=physbam.LA_i()
fixed_weights=physbam.LA_Vf3()
grid=physbam.GRID_Vf3()
phi=physbam.ARRAYS_3D_f()
levelset=physbam.LEVELSET_UNIFORM_Vf3(grid,phi)

physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/particles",particles)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/edges.curve",edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes",fixed_nodes)
#physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_tets",fixed_tets)
#physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_weights",fixed_weights)

physbam.Read_From_File("float",os.path.join(options.data_directory,"Rigid_Bodies",options.model_name+".tet"),tet)
physbam.Read_From_File("float",os.path.join(options.data_directory,"Rigid_Bodies",options.model_name+".0.phi"),levelset)
levelset_implicit_object=physbam.LEVELSET_IMPLICIT_OBJECT_Vf3(levelset.grid,levelset.phi)
levelset_implicit_object.Compute_Cell_Minimum_And_Maximum(False)
tet.Initialize_Hierarchy(True)

print "Calculating Fixed Nodes"
bar=ProgressBar.progressBar(1,len(edges.mesh.elements),80)
for i in range(1,len(edges.mesh.elements)+1):
    bar(i)
    edge=edges.mesh.elements[i]
    if levelset_implicit_object.Lazy_Inside(particles.X[edge[1]],0.) or levelset_implicit_object.Lazy_Inside(particles.X[edge[2]],0.):
        fixed_nodes.Append(edge[1])
        fixed_nodes.Append(edge[2])
physbam.Sort(fixed_nodes)
for i in range(len(fixed_nodes),1,-1):
    if (fixed_nodes[i]==fixed_nodes[i-1]): fixed_nodes.Remove_Index_Lazy(i)
physbam.Sort(fixed_nodes)
print ""

print "Calculating Embedded Tets"
bar=ProgressBar.progressBar(1,len(fixed_nodes),80)
tets=physbam.LA_i()
tolerance=1e-5
failure=False
for i in range(1,len(fixed_nodes)+1):
    bar(i)
    tets.Remove_All()
    tet.hierarchy.Intersection_List(particles.X[fixed_nodes[i]],tets,tolerance)
    for tet_num in tets:
        t=physbam.TETRAHEDRON_f(*subset(tet.particles.X,tet.mesh.elements[tet_num]))
        bary=t.Barycentric_Coordinates(particles.X[fixed_nodes[i]])
        if (bary.x>-tolerance and bary.y>-tolerance and bary.z>-tolerance and bary.x+bary.y+bary.z-1<tolerance):
            fixed_tets.Append(tet_num)
            fixed_weights.Append(bary)
            break
print ""
        
physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes",fixed_nodes)
physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_tets",fixed_tets)
physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_weights",fixed_weights)
