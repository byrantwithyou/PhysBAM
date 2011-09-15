#!/usr/bin/python
from optparse import OptionParser
import physbam
import os
import math
import sys

#TEST NUMBERS
#1 - Twisting
#2 - Pushing
#3 - Buckling
#4 - Pulling
#5 - Adhesion
#6 - Projected Spring Test

def Curve(r,t,b,z):
    return physbam.Vf3(r*math.cos(b*t),r*math.sin(b*t),z*t);

def Distance(a,b):
    return math.pow(a.x-b.x,2)+math.pow(a.y-b.y,2)+math.pow(a.z-b.z,2);

usage="usage: %prog [options] <test number>"
parser=OptionParser(usage);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-l","--length",action="store",type="float",dest="default_length",default=.05);
parser.add_option("-d","--directory",action="store",type="string",dest="data_directory",default=None);
parser.add_option("-e","--segments",action="store",type="int",dest="n_segments",default=10);
parser.add_option("-r","--perturb",action="store",type="int",dest="perturb",default=1);
parser.add_option("-f","--fixed",action="store",type="int",dest="fixed",default=1);
parser.add_option("-a","--mass",action="store",type="float",dest="mass_per_unit",default=100);
parser.add_option("-c","--curly",action="store_true",dest="curly",default=False);
parser.add_option("-v","--curve",action="store",type="float",dest="curve",default=0);
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

mass_per_unit=options.mass_per_unit;
test_number=int(args[0]);

perturb_threshold=math.fabs(1-math.cos(25*math.pi/180));

# make meshes for each type of hair
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

import random
n_segments=options.n_segments;
default_length=options.default_length;
mass_per_hair=mass_per_unit*default_length;
mass_per_segment=mass_per_hair/n_segments;

fixed_nodes_start=physbam.LA_i()
fixed_nodes_end=physbam.LA_i()
fixed_positions=physbam.LA_Vf3()
start_points=physbam.LA_Vf3()
normals=physbam.LA_Vf3()
curve=physbam.LA_Vf3()
masses=physbam.LA_f()

particles.Increase_Array_Size(2*(n_segments*3))
if(test_number==1 or test_number==2 or test_number==6):
    start_points.Append(physbam.Vf3(0,0,0));
    normals.Append(physbam.Vf3(1,0,0));
    curve.Append(physbam.Vf3(0,options.curve*default_length,0));
if(test_number==3):
    start_points.Append(physbam.Vf3(0,1,0));
    normals.Append(physbam.Vf3(1,0,0));
    curve.Append(physbam.Vf3(0,options.curve*default_length,0));
if(test_number==4):
    start_points.Append(physbam.Vf3(0,1,0));
    normals.Append(physbam.Vf3(0,-1,0));
    curve.Append(physbam.Vf3(options.curve*default_length,0,0));
if(test_number==5): 
    start_points.Append(physbam.Vf3(0,0,0));
    #start_points.Append(physbam.Vf3(default_length/5,default_length/5,-default_length/20));
    start_points.Append(physbam.Vf3(default_length/5,default_length/5,0));
    normals.Append(physbam.Vf3(1,0,0));
    #normals.Append(physbam.Vf3(0,-2,1).Normalized());
    normals.Append(physbam.Vf3(0,-1,0).Normalized());
    curve.Append(physbam.Vf3(0,options.curve*default_length,0));
    curve.Append(physbam.Vf3(options.curve*default_length,0,0));
if(test_number==3 or test_number==6): control_end=False;
else: control_end=True;
threshold=options.fixed;

for j in range(1,len(start_points)+1):
    if (j>1): threshold=1;
    #perturb_amount=.05*default_length
    perturb_amount=.02*default_length
    points=[]
    if options.curly: 
        frame=physbam.FRAME_Vf3();
        frame.r=physbam.ROTATION_Vf3.From_Rotated_Vector(physbam.Vf3(0,0,1),normals[j]);
        frame.t=start_points[j];
        for i in range(1,n_segments+1): points.append(frame*Curve(perturb_amount,i-1,math.pi/3,default_length/n_segments))
    else:
        points=map(lambda t: start_points[j]+normals[j]*default_length*float(t-1)/n_segments, range(1,n_segments+1))
    if options.curve:
        for i in range(len(points)):
            if(i<threshold or (control_end and i>len(points)-threshold-1)): continue;
            if(control_end):
                x=float(i-threshold+1)/float(len(points)-threshold*2+1);
            else:
                x=float(i-threshold+1)/float(len(points)-threshold);
            points[i]+=curve[j]*4*x*(1-x)
    # find segments that need to be perturbed
    segment_perturbed=[False]
    for i in range(len(points)-1):
        v1,v2=points[i]-points[i-1],points[i+1]-points[i-1]
        v1.Normalize();v2.Normalize()
        if i>0 and (v1.Magnitude()*v2.Magnitude()-physbam.Vf3.Dot_Product(v1,v2))<perturb_threshold: segment_perturbed.append(True)
        elif i<len(points)-2 and (v1.Magnitude()*v2.Magnitude()-physbam.Vf3.Dot_Product(v1,v2))<perturb_threshold: segment_perturbed.append(True)
        else: segment_perturbed.append(False)

    # Incrementally add original segments 
    i=1
    last_particle=particles.Add_Particle()
    particles.X[last_particle]=points[0]
    if(threshold>0): fixed_nodes_start.Append(last_particle);
    masses.Append(0);
    previous=[()]
    last=None
    perturb=(points[0]-points[1]).Orthogonal_Vector()
    perturb_amount=options.perturb
    while i<len(points):
        if not options.curly or i<=threshold:
        #if segment_perturbed[i]: ## segment has been subdivided
            # find perturb direction
            segment_direction=(points[i-1]-points[i]).Normalized()
            #perturb=physbam.ROTATION_Vf3(random.uniform(45.*math.pi/180.,135.*math.pi/180.),segment_direction).Rotate(perturb)
            perturb=physbam.ROTATION_Vf3(math.pi/2,segment_direction).Rotate(perturb)
            perturb.Project_Orthogonal_To_Unit_Direction(segment_direction)

            # make perturbed particle and real particle
            e0,r0=particles.Add_Particle(),particles.Add_Particle()
            masses[len(masses)]+=mass_per_segment/3.;
            masses.Append(mass_per_segment/3.);
            masses.Append(mass_per_segment/3.);
            particles.X[e0]=.5*(points[i-1]+points[i])+perturb_amount*perturb
            particles.X[r0]=points[i]
            # make structural edges
            if (i<threshold):
                fixed_nodes_start.Append(e0); 
                fixed_nodes_start.Append(r0);
                fixed_edges.mesh.elements.Append(last_particle,r0)
            elif (i>len(points)-threshold and control_end):
                fixed_nodes_end.Append(e0); 
                fixed_nodes_end.Append(last_particle);
                fixed_edges.mesh.elements.Append(last_particle,r0)
                
                if len(previous[i-1])==3 and i==len(points)-threshold+1: # Case 2
                    r2,e1,r1=previous[i-1]
                    torsion_edges.mesh.elements.Append(r2,e0)
                    volume.mesh.elements.Append(r2,e1,r1,e0)
                    torsion_edges.mesh.elements.Append(e1,r0)
                    volume.mesh.elements.Append(e1,r1,r0,e0)
                    triangles.mesh.elements.Append(r1,e0,r0)

                    #bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
                    bending_edges.mesh.elements.Append(e1,e0) # opposite spring
                elif len(previous[i-1])==2: # Case 1
                    r2,r1=previous[i-1]
                    if i==len(points)-theshold+1:
                        torsion_edges.mesh.elements.Append(e0,r2)
                        volume.mesh.elements.Append(r2,r1,e0,r0)
                        bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
                    if i-2>=0 and len(previous[i-2])==2 and (i==len(points)-threshold+1 or i==len(points)-threshold+2): # Case 0
                        r3,foo=previous[i-2]
                        torsion_edges.mesh.elements.Append(r3,r0)
                        volume.mesh.elements.Append(r3,r2,r1,r0)
                    if i==len(points)-theshold+1:
                        triangles.mesh.elements.Append(r1,e0,r0)
                        triangles.mesh.elements.Append(r2,r1,r0)
            else:
                edges.mesh.elements.Append(last_particle,r0)
                extra_edges.mesh.elements.Append(last_particle,e0)
                extra_edges.mesh.elements.Append(e0,r0)

                if len(previous[i-1])==3: # Case 2
                    r2,e1,r1=previous[i-1]
                    torsion_edges.mesh.elements.Append(r2,e0)
                    volume.mesh.elements.Append(r2,e1,r1,e0)
                    torsion_edges.mesh.elements.Append(e1,r0)
                    volume.mesh.elements.Append(e1,r1,r0,e0)
                    triangles.mesh.elements.Append(r1,e0,r0)

                    #bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
                    bending_edges.mesh.elements.Append(e1,e0) # opposite spring
                elif len(previous[i-1])==2: # Case 1
                    r2,r1=previous[i-1]
                    torsion_edges.mesh.elements.Append(e0,r2)
                    volume.mesh.elements.Append(r2,r1,e0,r0)
                    bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
                    if i-2>=0 and len(previous[i-2])==2: # Case 0
                        r3,foo=previous[i-2]
                        torsion_edges.mesh.elements.Append(r3,r0)
                        volume.mesh.elements.Append(r3,r2,r1,r0)
                    triangles.mesh.elements.Append(r1,e0,r0)
                    triangles.mesh.elements.Append(r2,r1,r0)

            previous.append((last_particle,e0,r0))
            
        else: ## segment not subdivided
            # make only real particle 
            r0=particles.Add_Particle()
            masses[len(masses)]+=mass_per_segment/2.;
            masses.Append(mass_per_segment/2.);
            particles.X[r0]=points[i]
            # only one real edge
            if (i<threshold):
                fixed_nodes_start.Append(r0); 
                fixed_edges.mesh.elements.Append(last_particle,r0)
            else:
                if (i>len(points)-threshold and control_end):
                    fixed_nodes_end.Append(last_particle);
                    fixed_edges.mesh.elements.Append(last_particle,r0)
                
                edges.mesh.elements.Append(last_particle,r0)

                if len(previous[i-1])==2: # Case 3 & 5
                    r2,r1=previous[i-1]
                    triangles.mesh.elements.Append(r2,r1,r0)
                    bending_edges.mesh.elements.Append(r2,r0)
                    if i-2>=0:
                        if len(previous[i-2])==2: # Case 3
                            r3,dummy=previous[i-2]
                            torsion_edges.mesh.elements.Append(r3,r0)
                            volume.mesh.elements.Append(r3,r2,r1,r0)
                        elif len(previous[i-2])==3: # Case 5
                            r3,e2,dummy=previous[i-2]
                            torsion_edges.mesh.elements.Append(r3,r0)
                            volume.mesh.elements.Append(r3,r2,r1,r0)
                elif len(previous[i-1])==3: # Case 4
                    r2,e1,r1=previous[i-1]
                    triangles.mesh.elements.Append(r2,r1,r0)
                    bending_edges.mesh.elements.Append(r2,r0)
                    torsion_edges.mesh.elements.Append(e1,r0)
                    volume.mesh.elements.Append(r2,e1,r1,r0)
            
            previous.append((last_particle,r0))

        last_particle=r0
        i+=1
        if(i==len(points) and threshold>0 and control_end): fixed_nodes_end.Append(last_particle)

# Update number nodes and write meshes
edges.Update_Number_Nodes()
fixed_edges.Update_Number_Nodes()
extra_edges.Update_Number_Nodes()
bending_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
volume.Update_Number_Nodes()
triangles.Update_Number_Nodes()
os.system("mkdir -p "+options.data_directory+"/"+options.sim_name);
if not options.pretend:
    open(options.data_directory+"/"+options.sim_name+"/cmdline.txt","w").write(" ".join(sys.argv)+"\n")
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/particles",particles)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/edges.curve",edges.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_edges.curve",fixed_edges.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/extra_edges.curve",extra_edges.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/bending_edges.curve",bending_edges.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/torsion_edges.curve",torsion_edges.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/triangles.tri",triangles.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/tets",volume.mesh)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes_start",fixed_nodes_start)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/fixed_nodes_end",fixed_nodes_end)
    physbam.Write_To_File("float",options.data_directory+"/"+options.sim_name+"/masses",masses)
    print "Files written to "+options.data_directory+"/"+options.sim_name+"/";
else:
    print "Files would have been written to "+options.data_directory+"/"+options.sim_name+"/";

