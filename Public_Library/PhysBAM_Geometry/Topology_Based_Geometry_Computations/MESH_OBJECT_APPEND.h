//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __MESH_OBJECT_APPEND__
#define __MESH_OBJECT_APPEND__

#include <PhysBAM_Tools/Arrays/ARRAY.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/MESH_OBJECT.h>

namespace PhysBAM{
namespace TOPOLOGY_BASED_GEOMETRY_COMPUTATIONS
{
//#####################################################################
// Function Append_Particles_And_Create_Copy
//#####################################################################
template<class TV,class T_MESH>
STRUCTURE<TV>* Append_Particles_And_Create_Copy(const MESH_OBJECT<TV,T_MESH>& mo,GEOMETRY_PARTICLES<TV>& new_particles,ARRAY<int>* particle_indices) // number_nodes must be set elsewhere
{
    typename MESH_TO_OBJECT<TV,T_MESH>::TYPE* object=mo.Create(new_particles);
    int offset=new_particles.Size();
    new_particles.Append(mo.particles);
    if(particle_indices) for(int p=0;p<mo.particles.Size();p++) particle_indices->Append(p+offset);
    object->mesh.Initialize_Mesh_With_Particle_Offset(mo.mesh,offset);
    return object;
}
}
}
#endif
