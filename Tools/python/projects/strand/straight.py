#!/usr/bin/python
import math
import physbam

n_segments=100

#def X(s):
#    if s<.5:
#       return physbam.Vf3(0.,float(-s),0.)
#    else:
#       return physbam.Vf3(s-.5,float(-s),0.)

#def X(s):
#    return physbam.Vf3(math.cos(s*10)/10.,-s/5,math.sin(s*10)/10.)
def X(s):
    return physbam.Vf3(s,0,0)

def s(i):
    return float(i-1)/(n_segments)

perturb_threshold=1e-4
points=map(lambda t: X(s(t)),range(1,n_segments+2))
#points=[physbam.Vf3(0,0,0),
#        physbam.Vf3(.2,0,0),
#        physbam.Vf3(.4,.2,0),
#        physbam.Vf3(.6,.2,0),
#        physbam.Vf3(.8,.2,0)]

# find segments that need to be perturbed
segment_perturbed=[False]
for i in range(len(points)-1):
    if i>0 and physbam.TRIANGLE_3D_f(points[i-1],points[i],points[i+1]).Area()<perturb_threshold: segment_perturbed.append(True)
    elif i<len(points)-2 and physbam.TRIANGLE_3D_f(points[i],points[i+1],points[i+2]).Area()<perturb_threshold: segment_perturbed.append(True)
    else: segment_perturbed.append(False)
#print segment_perturbed        

# make meshes for each type of hair
particles=physbam.SOLIDS_PARTICLE_Vf3()
edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
extra_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
bending_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
torsion_edges=physbam.SEGMENTED_CURVE_Vf3.Create(particles)
triangles=physbam.TRIANGULATED_SURFACE_f.Create(particles)
volume=physbam.TETRAHEDRALIZED_VOLUME_f.Create(particles)

# Incrementally add original segments 
i=1
last_particle=particles.Add_Particle()
particles.X[last_particle]=points[0]
previous=[()]
last=None
perturb=(points[0]-points[1]).Orthogonal_Vector()

perturb_amount=1
while i<len(points):
    #print i
    if segment_perturbed[i]: ## segment has been subdivided
        # find perturb direction
        segment_direction=(points[i-1]-points[i]).Normalized()
        perturb=physbam.ROTATION_Vf3(math.pi/2,segment_direction).Rotate(perturb)
        perturb.Project_Orthogonal_To_Unit_Direction(segment_direction)

        # make perturbed particle and real particle
        e0,r0=particles.Add_Particle(),particles.Add_Particle()
        particles.X[e0]=.5*(points[i-1]+points[i])+perturb_amount*perturb
        particles.X[r0]=points[i]
        # make structural edges
        edges.mesh.elements.Append(physbam.Vi2(last_particle,r0))
        extra_edges.mesh.elements.Append(physbam.Vi2(last_particle,e0))
        extra_edges.mesh.elements.Append(physbam.Vi2(e0,r0))
        
        if len(previous[i-1])==3: # Case 2
            r2,e1,r1=previous[i-1]
            torsion_edges.mesh.elements.Append(physbam.Vi2(r2,e0))
            volume.mesh.elements.Append(physbam.Vi4(r2,e1,r1,e0))
            torsion_edges.mesh.elements.Append(physbam.Vi2(e1,r0))
            volume.mesh.elements.Append(physbam.Vi4(e1,r1,r0,e0))
            triangles.mesh.elements.Append(physbam.Vi3(r1,e0,r0))

            #bending_edges.mesh.elements.Append(last_particle,r0) # guessed spring
            bending_edges.mesh.elements.Append(physbam.Vi2(e1,e0)) # opposite spring
        elif len(previous[i-1])==2: # Case 1
            r2,r1=previous[i-1]
            torsion_edges.mesh.elements.Append(physbam.Vi2(e0,r2))
            volume.mesh.elements.Append(physbam.Vi4(r2,r1,e0,r0))
            bending_edges.mesh.elements.Append(physbam.Vi2(last_particle,r0)) # guessed spring
            if i-2>=0 and len(previous[i-2])==2: # Case 0
                r3,foo=previous[i-2]
                torsion_edges.mesh.elements.Append(physbam.Vi2(r3,r0))
                volume.mesh.elements.Append(physbam.Vi4(r3,r2,r1,r0))
            triangles.mesh.elements.Append(physbam.Vi3(r1,e0,r0))
            triangles.mesh.elements.Append(physbam.Vi3(r2,r1,r0))


        previous.append(physbam.Vi3(last_particle,e0,r0))
        last_particle=r0
        i+=1
        
    else: ## segment not subdivided
        perturb=(point[i-1]-points[i]).Orthogonal_Vector()
        # make only real particle 
        r0=particles.Add_Particle()
        particles.X[r0]=points[i]
        # only one real edge
        edges.mesh.elements.Append(physbam.Vi2(last_particle,r0))

        triangles.mesh.elements.Append(physbam.Vi3(r2,r1,r0))
        if len(previous[i-1])==2: # Case 3 & 5
            r2,r1=previous[i-1]
            bending_edges.mesh.elements.Append(physbam.Vi2(r2,r0))
            if i-2>=0:
                if len(previous[i-2])==2: # Case 3
                    r3,dummy=previous[i-2]
                    torsion_edges.mesh.elements.Append(physbam.Vi2(r3,r0))
                    volume.mesh.elements.Append(physbam.Vi4(r3,r2,r1,r0))
                elif len(previous[i-2])==3: # Case 5
                    r3,e2,dummy=previous[i-2]
                    torsion_edges.mesh.elements.Append(physbam.Vi2(r3,r0))
                    volume.mesh.elements.Append(physbam.Vi4(r3,r2,r1,r0))
        elif len(previous[i-1])==3: # Case 4
            r2,e1,r1=previous[i-1]
            bending_edges.mesh.elements.Append(physbam.Vi2(r2,r0))
            torsion_edges.mesh.elements.Append(physbam.Vi2(e1,r0))
            volume.mesh.elements.Append(physbam.Vi4(r2,e1,r1,r0))
        
        previous.append(physbam.Vi2(last_particle,r0))
        last_particle=r0
        i+=1

# Update number nodes and write meshes
edges.Update_Number_Nodes()
extra_edges.Update_Number_Nodes()
bending_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
torsion_edges.Update_Number_Nodes()
volume.Update_Number_Nodes()
triangles.Update_Number_Nodes()
physbam.Write_To_File("float","edges.curve",edges)
physbam.Write_To_File("float","extra_edges.curve",extra_edges)
physbam.Write_To_File("float","bending_edges.curve",bending_edges)
physbam.Write_To_File("float","torsion_edges.curve",torsion_edges)
physbam.Write_To_File("float","triangles.tri",triangles)
physbam.Write_To_File("float","tets",volume)
