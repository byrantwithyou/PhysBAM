#!/usr/bin/python
from optparse import OptionParser
import physbam
import os
import math
import ProgressBar
import sys
import random

def curve(r,t,b,z):
    return physbam.Vf3(r*math.cos(b*t),r*math.sin(b*t),z*t);

def Distance(a,b):
    return math.pow(a.x-b.x,2)+math.pow(a.y-b.y,2)+math.pow(a.z-b.z,2);

usage="usage: %prog [options] <data directory>"
parser=OptionParser(usage);
parser.add_option("-m","--model_name",action="store",type="string",dest="model_name",default=None);
parser.add_option("-s","--sim_name",action="store",type="string",dest="sim_name",default=None);
parser.add_option("-f","--guide_sim",action="store",type="string",dest="guide_name",default=None);
parser.add_option("-t","--proportional",action="store",type="string",dest="prop",default=None);
parser.add_option("-l","--length",action="store",type="float",dest="default_length",default=.05);
parser.add_option("-n","--hair",action="store",type="int",dest="n",default=100);
parser.add_option("-e","--segments",action="store",type="int",dest="n_segments",default=10);
parser.add_option("-a","--mass",action="store",type="int",dest="mass_per_unit",default=1000);
parser.add_option("-r","--samples",action="store",type="int",dest="num_interp",default=10);
parser.add_option("-v","--curve",action="store",type="float",dest="curve",default=0);
parser.add_option("-c","--curly",action="store_true",dest="curly",default=False);
parser.add_option("-p","--pretend",action="store_true",dest="pretend",default=False);
parser.add_option("-i","--ignore",action="store_true",dest="ignore_len",default=False);
parser.add_option("-g","--guide",action="store_true",dest="guide",default=False);
parser.add_option("-x","--static",action="store_true",dest="static",default=False)
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

if options.guide and options.guide_name is not None:
    parser.error("You cannot create guide hairs from a guide hair file")
    sys.exit(1);

if options.guide_name is not None and not os.path.isfile(args[0]+"/"+options.guide_name+"/surface_points.gz"):
    parser.error("%s was not properly generated"%options.guide_name)
    sys.exit(1);

if options.curly and options.guide:
    parser.error("Guide hairs must be straight");
    sys.exit(1);

mass_per_unit=options.mass_per_unit;
if options.guide: iterations=20
else: iterations=2;
if options.n<100: iterations=0;
vertex_densities=physbam.A_f()
physbam.Read_From_File("float",os.path.join(args[0],"maps","densities.gz"),vertex_densities)
#vertex_colors=physbam.A_i()
#physbam.Read_From_File("float",os.path.join(MODELDIR,"maps","particle_colors.gz"),vertex_colors)
hair_densities=physbam.A_f()
physbam.Read_From_File("float",os.path.join(args[0],"maps","hair_densities.gz"),hair_densities)
comb_densities=physbam.A_f()
physbam.Read_From_File("float",os.path.join(args[0],"maps","comb_densities.gz"),comb_densities)
tri=physbam.TRIANGULATED_SURFACE_f.Create()
physbam.Read_From_File("float",os.path.join(args[0],"Rigid_Bodies",options.model_name+".tri"),tri)
base_vertex_densities=physbam.A_f()
if options.prop is not None:
    physbam.Read_From_File("float",os.path.join(options.prop,"maps","densities.gz"),base_vertex_densities)

surface_points=physbam.A_Vf3();
surface_triangles=physbam.A_i();
guide_points=physbam.A_Vf3();
guide_tris=physbam.A_i();
# cell centered arrays
comb_elements=physbam.A_f();
comb_elements.Resize(len(tri.mesh.elements),False,False)
ordering=physbam.A_f();
ordering.Resize(len(tri.mesh.elements),False,False)
element_densities=physbam.A_f()
element_densities.Resize(len(tri.mesh.elements),False,False)
base_element_densities=physbam.A_f()
base_element_densities.Resize(len(tri.mesh.elements),False,False)
element_lengths=physbam.A_f()
element_lengths.Resize(len(tri.mesh.elements),False,False)

print "Computing interpolation points"
bar=ProgressBar.progressBar(1,len(tri.mesh.elements),80)
hash=[[] for i in range(options.num_interp)]
interpolation_points=physbam.LA_i();
for element in range(1,len(tri.mesh.elements)+1):
    bar(element)
    nodes=tri.mesh.elements[element]
    comb_elements[element]=min(map(lambda x: comb_densities[x],nodes))
    if comb_elements[element]>0:
        comb_elements[element]=sum(map(lambda x: comb_densities[x],nodes))/3.
        if comb_elements[element]==1: hash[options.num_interp-1].append(element);
        else: hash[int(comb_elements[element]*options.num_interp)].append(element);

for triangles in hash:
    center=physbam.Vf3();
    for index in triangles:
        triangle=physbam.TRIANGLE_3D_f(*map(lambda node: tri.particles.X[node],tri.mesh.elements[index]))
        center+=triangle.Center();
    center/=len(triangles);
    min_tri=-1;
    min_distance=-1;
    for index in triangles:
        triangle=physbam.TRIANGLE_3D_f(*map(lambda node: tri.particles.X[node],tri.mesh.elements[index]))
        distance=math.pow(triangle.Center().x-center.x,2)+math.pow(triangle.Center().y-center.y,2)+math.pow(triangle.Center().z-center.z,2);
        if (distance < min_distance or min_distance==-1):
            min_distance=distance;
            min_tri=index;
    interpolation_points.Append(min_tri);
print "";
    
print "Computing Cell Centered Data"
bar=ProgressBar.progressBar(1,len(tri.mesh.elements),80)
length=len(tri.mesh.elements);
active=0;
base_active=0;
for element in range(len(tri.mesh.elements),0,-1):
    bar(length-element)
    nodes=tri.mesh.elements[element]
    element_densities[element]=sum(map(lambda x: vertex_densities[x],nodes))/3.
    if options.prop is not None:
        base_element_densities[element]=sum(map(lambda x: base_vertex_densities[x],nodes))/3.
        if base_element_densities[element]!=0: base_active+=1;
    if element_densities[element]==0: tri.mesh.elements.Remove_Index_Lazy(element)
    else: active+=1;
old_n=options.n
if options.prop is not None: 
    options.n=options.n*active/base_active;

for element in range(1,len(tri.mesh.elements)):
    nodes=tri.mesh.elements[element]
    element_densities[element]=sum(map(lambda x: vertex_densities[x],nodes))/3.
    element_lengths[element]=sum(map(lambda x: hair_densities[x],nodes))/3.
    triangle=physbam.TRIANGLE_3D_f(*map(lambda node: tri.particles.X[node],tri.mesh.elements[element]))
    element_densities[element]*=triangle.Area()

print " "

print "Creating point distibution"
# TODO: make an option for not switching to the repulsion way
if old_n>100:
    repulsions=physbam.POINT_REPULSION_f(tri,options.n)
    repulsions.Seed_Points()
    for i in range(iterations):
        repulsions.Move_Points()
    repulsions.Get_Points(surface_points);
    repulsions.Get_Triangles(surface_triangles);
else:
    tri.Update_Vertex_Normals();
    bar=ProgressBar.progressBar(1,len(tri.mesh.elements),80)
    surface_points.Resize(options.n,True,True)
    surface_triangles.Resize(options.n,True,True)
    for i in range(1,options.n+1):
        bar(i);
        random_number=random.randint(1,len(tri.mesh.elements))
        surface_triangles[i]=random_number;
        triangle=physbam.TRIANGLE_3D_f(*map(lambda node: tri.particles.X[node],tri.mesh.elements[random_number]))
        x1=random.uniform(0,1)
        x2=random.uniform(0,x1)
        surface_points[i]=triangle.Point_From_Barycentric_Coordinates(physbam.Vf3(1-x1-x2,x1,x2));
if options.guide_name is not None:
    dist_thresh=-1;
    for i in range(1,len(surface_points)+1):
        for j in range(i+1,len(surface_points)+1):
            if Distance(surface_points[i],surface_points[j])<dist_thresh or dist_thresh==-1:
                dist_thresh=Distance(surface_points[i],surface_points[j]);
print " "

print "Creating hair"
if options.curly: perturb_threshold=math.fabs(1-math.cos(1*math.pi/180));
else: perturb_threshold=math.fabs(1-math.cos(45*math.pi/180));

# make meshes for each type of hair
guide_particles=physbam.SOLIDS_PARTICLE_Vf3()
particles=physbam.SOLIDS_PARTICLE_Vf3()
tet_particles=physbam.SOLIDS_PARTICLE_Vf3()
particles.Store_Velocity(True);
particles.Store_Mass(True);
edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
fixed_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
bending_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
torsion_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
triangles=physbam.TRIANGULATED_SURFACE_f.Create(particles)
volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)
tet=physbam.TETRAHEDRALIZED_VOLUME_f.Create(tet_particles)
guide_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
guide_extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
guide_volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)

if options.guide_name is not None:
    physbam.Read_From_File("float",args[0]+"/"+options.guide_name+"/extra_edges.curve",guide_edges.mesh)
    physbam.Read_From_File("float",args[0]+"/"+options.guide_name+"/particles",guide_particles)
    physbam.Read_From_File("float",args[0]+"/"+options.guide_name+"/surface_points",guide_points)
    physbam.Read_From_File("float",args[0]+"/"+options.guide_name+"/surface_triangles",guide_tris)

import random
n=options.n;
n_segments=options.n_segments;
length=n_segments;
default_length=options.default_length;
mass_per_hair=mass_per_unit*default_length;
mass_per_segment=mass_per_hair/n_segments;
r=random.Random()
bar=ProgressBar.progressBar(1,n+1,80)
fixed_nodes=physbam.LA_i()
fixed_tets=physbam.LA_i()
fixed_weights=physbam.LA_Vf3()
masses=physbam.LA_f()
#particles.Increase_Array_Size(n*(n_segments*3))
tri.Use_Vertex_Normals();
skipped_hairs=0
for j in range(1,len(surface_points)+1):
    bar(j)
    threshold=2
    skip=False
    for point in guide_points:
        if Distance(point,surface_points[j])<dist_thresh: 
            skip=True
            continue
    if skip:
        skipped_hairs+=1
        continue
    point=surface_points[j];
    normal=tri.Normal(point,surface_triangles[j]);
    root=None
    guide_index=0;
    min_dist=-1;
    for i in range(1,len(guide_points)+1):
        dist=Distance(point,guide_points[i])
        if dist<min_dist or min_dist==-1:
            min_dist=dist;
            normal=tri.Normal(guide_points[i],guide_tris[i]);
            root=guide_points[i];
    for i in range(1,len(guide_edges.mesh.elements)+1):
        if guide_particles.X[guide_edges.mesh.elements[i][1]]==root: guide_index=i-1;
    if not options.ignore_len: length=int(n_segments*element_lengths[j]);
    random_nums=random.Random();
    if physbam.Vf3.Dot_Product(physbam.Vf3(0,0,1),normal)<0:
        offset_normal=physbam.Vf3(0,0,0);
    else:
        offset_normal=physbam.Vf3(1,0,0);
        if physbam.Vf3.Dot_Product(offset_normal,normal)<0:
            offset_normal=-offset_normal;
        projected=physbam.Vf3(normal.x,0,normal.z);
        projected.Normalize();
        if(math.cos(15*math.pi/180)<physbam.Vf3.Dot_Product(physbam.Vf3(0,0,1),projected) and random_nums.uniform(0,1)>.5):
            offset_normal=-offset_normal;
    if physbam.Vf3.Dot_Product(physbam.Vf3(0,1,0),normal)>0:
        normal+=offset_normal;
        normal.Normalize();
    points=map(lambda t: point+normal*default_length*float(t-1)/n_segments, range(-1,length+2))
    perturb_amount=.05*default_length
    mod=threshold%2;
    if options.curly: 
        points=map(lambda t: point+normal*default_length*float(t-1)/n_segments, range(-1,threshold))
        frame=physbam.FRAME_Vf3();
        frame.r=physbam.ROTATION_Vf3.From_Rotated_Vector(physbam.Vf3(0,0,1),normal);
        frame.t=points[threshold-1];
        for i in range(threshold,length+2): points.append(frame*curve(perturb_amount,i-threshold+1,math.pi/3,default_length/n_segments))
    else:
        points=map(lambda t: point+normal*default_length*float(t-1)/n_segments, range(-1,length+2))
    # find segments that need to be perturbed
    segment_perturbed=[False]
    for i in range(len(points)-1):
        v1,v2=points[i]-points[i-1],points[i+1]-points[i-1]
        v1.Normalize();v2.Normalize()
        if i>0 and (v1.Magnitude()*v2.Magnitude()-physbam.Vf3.Dot_Product(v1,v2))<perturb_threshold: segment_perturbed.append(True)
        elif i<len(points)-2 and (v1.Magnitude()*v2.Magnitude()-physbam.Vf3.Dot_Product(v1,v2))<perturb_threshold: segment_perturbed.append(True)
        else: segment_perturbed.append(False)

    if(len(points)==0): continue
    # Incrementally add original segments 
    i=1
    last_particle=particles.Add_Particle()
    print "\n\n\n%s"%(repr(points)) # %d of %d or %d"%(last_particle,len(particles),len(particles.X),len(points))
    #print particles.X[1]
    particles.X[last_particle]=points[0]
    fixed_nodes.Append(last_particle);
    masses.Append(0);
    previous=[()]
    last=None
    perturb=(points[0]-points[1]).Orthogonal_Vector()
    threshold=2;
    perturb_amount=1
    guide_skip=1;
    normal=tri.Normal(surface_points[j],surface_triangles[j]);
    projection=normal;
    projection.y=0;
    direction=physbam.Vf3.Cross_Product(projection,physbam.Vf3(0,1,0));
    direction=physbam.Vf3.Cross_Product(direction,normal);
    direction.Normalize();
    if(physbam.Vf3.Dot_Product(direction,physbam.Vf3(0,-1,0))>0):
        direction=-direction;
    while i<len(points):
        if(i>threshold):
            x=float(i-threshold)/float(len(points)-threshold-1);
            points[i]+=direction*4*x*(1-x)*options.curve;
        
        if not options.curly or i<=threshold:
        #if segment_perturbed[i]: ## segment has been subdivided
            # find perturb direction
            segment_direction=(points[i-1]-points[i]).Normalized()
            #perturb=physbam.ROTATION_Vf3(math.pi/2,segment_direction).Rotate(perturb)
            perturb=physbam.ROTATION_Vf3(random.uniform(45.*math.pi/180.,135.*math.pi/180.),segment_direction).Rotate(perturb)
            perturb.Project_Orthogonal_To_Unit_Direction(segment_direction)

            # make perturbed particle and real particle
            if(i<threshold):
                r0=particles.Add_Particle();
                masses[len(masses)]+=mass_per_segment/2.;
                masses.Append(mass_per_segment/2.);
                particles.X[r0]=points[i]
            else:
                e0,r0=particles.Add_Particle(),particles.Add_Particle()
                masses[len(masses)]+=mass_per_segment/3.;
                masses.Append(mass_per_segment/3.);
                masses.Append(mass_per_segment/3.);
                particles.X[e0]=.5*(points[i-1]+points[i])+perturb_amount*perturb
                particles.X[r0]=points[i]
            # make structural edges
            if (i<=threshold):
                if(i==threshold): fixed_nodes.Append(e0); 
                fixed_nodes.Append(r0);
                fixed_edges.mesh.elements.Append((last_particle,r0))
            else:
                edges.mesh.elements.Append((last_particle,r0))
                extra_edges.mesh.elements.Append((last_particle,e0))
                extra_edges.mesh.elements.Append((e0,r0))

                if len(previous[i-1])==3: # Case 2
                    r2,e1,r1=previous[i-1]
                    torsion_edges.mesh.elements.Append((r2,e0))
                    volume.mesh.elements.Append((r2,e1,r1,e0))
                    torsion_edges.mesh.elements.Append((e1,r0))
                    volume.mesh.elements.Append((e1,r1,r0,e0))
                    triangles.mesh.elements.Append((r1,e0,r0))

                    #bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
                    bending_edges.mesh.elements.Append((e1,e0)) # opposite spring
                elif len(previous[i-1])==2: # Case 1
                    r2,r1=previous[i-1]
                    torsion_edges.mesh.elements.Append((e0,r2))
                    volume.mesh.elements.Append((r2,r1,e0,r0))
                    bending_edges.mesh.elements.Append((last_particle,r0)) # guessed spring
                    if i-2>=0 and len(previous[i-2])==2: # Case 0
                        r3,foo=previous[i-2]
                        torsion_edges.mesh.elements.Append((r3,r0))
                        volume.mesh.elements.Append((r3,r2,r1,r0))
                    triangles.mesh.elements.Append((r1,e0,r0))
                    triangles.mesh.elements.Append((r2,r1,r0))

            if (i<threshold): previous.append((last_particle,r0))
            else: previous.append((last_particle,e0,r0))
            
        else: ## segment not subdivided
            # make only real particle 
            r0=particles.Add_Particle()
            masses[len(masses)]+=mass_per_segment/2.;
            masses.Append(mass_per_segment/2.);
            particles.X[r0]=points[i]
            # only one real edge
            if (i<=threshold):
                fixed_nodes.Append(r0); 
                fixed_edges.mesh.elements.Append((last_particle,r0))
            else:
                edges.mesh.elements.Append((last_particle,r0))

                if len(previous[i-1])==2: # Case 3 & 5
                    r2,r1=previous[i-1]
                    triangles.mesh.elements.Append((r2,r1,r0))
                    bending_edges.mesh.elements.Append((r2,r0))
                    if i-2>=0:
                        if len(previous[i-2])==2: # Case 3
                            r3,dummy=previous[i-2]
                            torsion_edges.mesh.elements.Append((r3,r0))
                            volume.mesh.elements.Append((r3,r2,r1,r0))
                        elif len(previous[i-2])==3: # Case 5
                            r3,e2,dummy=previous[i-2]
                            torsion_edges.mesh.elements.Append((r3,r0))
                            volume.mesh.elements.Append((r3,r2,r1,r0))
                elif len(previous[i-1])==3: # Case 4
                    r2,e1,r1=previous[i-1]
                    triangles.mesh.elements.Append((r2,r1,r0))
                    bending_edges.mesh.elements.Append((r2,r0))
                    torsion_edges.mesh.elements.Append((e1,r0))
                    volume.mesh.elements.Append((r2,e1,r1,r0))
            
            previous.append((last_particle,r0))

        if options.guide_name is not None and i%guide_skip==0 and i>threshold and i+1<len(points):
            guide_1=guide_edges.mesh.elements[guide_index+2*(i-threshold)][1]
            guide_2=guide_edges.mesh.elements[guide_index+2*(i-threshold)][2]
            guide_3=guide_edges.mesh.elements[guide_index+2*(i-threshold)+1][2]
            guide_extra_edges.mesh.elements.Append((r0,guide_1));
            guide_extra_edges.mesh.elements.Append((r0,guide_3));
            guide_extra_edges.mesh.elements.Append((r0,guide_2));
            guide_extra_edges.mesh.elements.Append((guide_1,guide_2));
            guide_extra_edges.mesh.elements.Append((guide_2,guide_3));
            guide_extra_edges.mesh.elements.Append((guide_1,guide_3));
            guide_volume.mesh.elements.Append((r0,guide_1,guide_2,guide_3));

        last_particle=r0
        i+=1
bar(n+1);
print " "

if not options.static:
    print "Calculating Embedded Tets"
    physbam.Read_From_File("float",os.path.join(args[0],"Rigid_Bodies",options.model_name+".tet"),tet)
    tet.Initialize_Hierarchy(True);
    tets=physbam.LA_i();
    tolerance=1e-5;
    bar=ProgressBar.progressBar(1,len(fixed_nodes),80)
    
    def subset(X,nodes):
        return map(lambda x: X[x],nodes)
    
    failure=False;
    for i in range(1,len(fixed_nodes)+1):
        bar(i)
        tets.Remove_All();
        tet.hierarchy.Intersection_List(particles.X[fixed_nodes[i]],tets,tolerance);
        for tet_num in tets:
            t=physbam.TETRAHEDRON_f(*subset(tet.particles.X,tet.mesh.elements[tet_num]))
            bary=t.Barycentric_Coordinates(particles.X[fixed_nodes[i]])
            if (bary.x>-tolerance and bary.y>-tolerance and bary.z>-tolerance and bary.x+bary.y+bary.z-1<tolerance):
                fixed_tets.Append(tet_num);
                fixed_weights.Append(bary);
                break
        if i!=len(fixed_tets): failure=True;
    if failure: print "ERROR: Tet intersections failed";

print " "
# Update number nodes and write meshes
edges.Update_Number_Nodes()
fixed_edges.Update_Number_Nodes()
extra_edges.Update_Number_Nodes()
bending_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
volume.Update_Number_Nodes()
triangles.Update_Number_Nodes()
guide_extra_edges.Update_Number_Nodes()
guide_volume.Update_Number_Nodes()
os.system("mkdir -p "+args[0]+"/"+options.sim_name);
if options.prop: print "Real n is "+str(n);
print "Skipped "+str(skipped_hairs)+" hairs";
if not options.pretend:
    open(args[0]+"/"+options.sim_name+"/cmdline.txt","w").write(" ".join(sys.argv)+"\n")
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/particles",particles)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/edges.curve",edges.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/fixed_edges.curve",fixed_edges.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/extra_edges.curve",extra_edges.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/bending_edges.curve",bending_edges.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/torsion_edges.curve",torsion_edges.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/triangles.tri",triangles.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/tets",volume.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/fixed_nodes",fixed_nodes)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/fixed_tets",fixed_tets)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/fixed_weights",fixed_weights)
    if options.guide: 
        physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/surface_points",surface_points)
        physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/surface_triangles",surface_triangles)
    if options.guide_name is not None:
        physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/guide_edges.curve",guide_extra_edges.mesh)
        physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/guide_tets",guide_volume.mesh)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/masses",masses)
    physbam.Write_To_File("float",args[0]+"/"+options.sim_name+"/interpolation",interpolation_points)
    print "Files written to "+args[0]+"/"+options.sim_name+"/";
else:
    print "Files would have been written to "+args[0]+"/"+options.sim_name+"/";

