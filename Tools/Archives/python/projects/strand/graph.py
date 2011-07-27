#!/usr/bin/python
from optparse import OptionParser
import sys

import physbam

usage="usage: %prog [options] <test number>"
parser=OptionParser(usage);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-l","--length",action="store",type="float",dest="default_length",default=.05);
parser.add_option("-d","--directory",action="store",type="string",dest="data_directory",default=None);
parser.add_option("-e","--segments",action="store",type="int",dest="n_segments",default=10);
parser.add_option("-a","--mass",action="store",type="int",dest="mass_per_unit",default=100);
parser.add_option("-c","--curly",action="store_true",dest="curly",default=False);
parser.add_option("-p","--pretend",action="store_true",dest="pretend",default=False);
(options,args)=parser.parse_args();
if len(args)!=1:
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


particles=physbam.SOLIDS_PARTICLE_Vf3()
particles.Store_Velocity(True);
particles.Store_Mass(True);
edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
fixed_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
bending_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
torsion_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
triangles=physbam.TRIANGULATED_SURFACE_f.Create(particles)
volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)

fixed_nodes_start=physbam.LA_i()
fixed_nodes_end=physbam.LA_i()
fixed_positions=physbam.LA_Vf3()
start_points=physbam.LA_Vf3()

physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/particles",particles)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/edges.curve",edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_edges.curve",fixed_edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/extra_edges.curve",extra_edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/bending_edges.curve",bending_edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/torsion_edges.curve",torsion_edges.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/triangles.tri",triangles.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/tets",volume.mesh)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes_start",fixed_nodes_start)
physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes_end",fixed_nodes_end)
#physbam.Read_From_File("float",options.data_directory+"/"+options.sim_name+"/masses",masses)


epoch=physbam.A_i();epoch.Resize(len(particles),True,True)

epoch_members=physbam.LA_LA_i();
depends=physbam.LA_LA_i();depends.Resize(len(particles),True,True)
depends_priority=physbam.LA_LA_i();depends_priority.Resize(len(particles),True,True)
depends_lengths=physbam.LA_LA_f();depends_lengths.Resize(len(particles),True,True)

print "edge"
print edges.mesh.elements
print "extra_edges"
print extra_edges.mesh.elements
print "fixed"
print fixed_nodes_start
print fixed_nodes_end


project_curve=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
project_mesh=project_curve.mesh

meshes=[edges.mesh,extra_edges.mesh]
priority=[0,1,2,3]
for i in meshes: i.Initialize_Neighbor_Nodes()

max_epoch=1
for mynode in fixed_nodes_start:
    epoch[mynode]=1
#for mynode in fixed_nodes_end:
#    epoch[mynode]=1

print "Calculating dependencies"
for mynode in xrange(1,len(particles)+1):
    if epoch[mynode]==1:
        continue # was a fixed node
    depends=[]
    for mesh_index in xrange(len(meshes)):
        mesh=meshes[mesh_index]
        for other_node in mesh.neighbor_nodes[mynode]:
            if epoch[other_node]!=0:
                #print "pri %s "%repr(priority)
                if priority[mesh_index]!=0 and priority[mesh_index]!=1:
                    print "Expected only priority 0 and 1 meshes"
                    sys.exit(1)
                depends.append((mynode,other_node,priority[mesh_index]))
    epoch[mynode]=1
    #if len(depends)==0:
    #    print "Expected dependency for node %d"%mynode
    #    sys.exit(1)
    def minimum_priority(x,y):
        if x[2]<y[2]: return x
        else: return y
    if len(depends)!=0:
        mynode,other_node,element_priority=reduce(minimum_priority,depends)
        project_mesh.elements.Append(mynode,other_node)
project_curve.Update_Number_Nodes()
print "";

print project_curve.mesh.elements

physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/project_mesh",project_mesh)
