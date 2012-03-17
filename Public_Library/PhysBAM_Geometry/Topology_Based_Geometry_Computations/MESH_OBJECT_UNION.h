//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MESH_OBJECT_UNION__
#define __MESH_OBJECT_UNION__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>

namespace PhysBAM{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS{
//#####################################################################
// Function Union_Mesh_Objects_Relatively
//#####################################################################
template<class TV,class T_OBJECT>
void Union_Mesh_Objects_Relatively(T_OBJECT *object,const ARRAY<T_OBJECT*>& object_list,const ARRAY<FRAME<TV> >& relative_frames)
{
    GEOMETRY_PARTICLES<TV>& particles=object->particles;
    particles.Clean_Memory();
    object->mesh.elements.Remove_All();
    // resize
    {int total_particles=0,total_elements=0;
    for(int i=0;i<object_list.m;i++){
        total_particles+=object_list(i)->particles.Size();
        total_elements+=object_list(i)->mesh.elements.m;}
    particles.Preallocate(total_particles);object->mesh.elements.Preallocate(total_elements);}
    // copy
    for(int i=0;i<object_list.m;i++){
        int particle_offset=particles.Size();
        particles.Add_Arrays(object_list(i)->particles);
        particles.Append(object_list(i)->particles);
        for(int p=0;p<object_list(i)->particles.Size();p++){int p2=p+particle_offset;
            particles.X(p2)=relative_frames(i)*particles.X(p2);
            if(particles.store_velocity) particles.V(p2)=relative_frames(i).r.Rotate(particles.V(p2));}
        object->mesh.elements.Append_Elements(object_list(i)->mesh.elements+particle_offset);}
    object->Update_Number_Nodes();
}

}
}
#endif
