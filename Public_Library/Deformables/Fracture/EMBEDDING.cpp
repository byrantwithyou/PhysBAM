//#####################################################################
// Copyright 2009, Jon Gretarsson.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class EMBEDDING
//#####################################################################
#include <Geometry/Registry/STRUCTURE_REGISTRY.h>
#include <Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <Deformables/Fracture/EMBEDDED_TRIANGULATED_OBJECT.h>
#include <Deformables/Fracture/EMBEDDING.h>
#include <Deformables/Fracture/TRIANGLES_OF_MATERIAL.h>
#include <Deformables/Particles/FREE_PARTICLES.h>
using namespace PhysBAM;
namespace PhysBAM{
bool Register_Solids_Structures(){
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<EMBEDDING<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<float> >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDING<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<EMBEDDING<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE<double> >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDING<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<TRIANGLES_OF_MATERIAL<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<EMBEDDED_TRIANGULATED_OBJECT<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,3> >::Register<FREE_PARTICLES<VECTOR<double,3> > >();
    STRUCTURE_REGISTRY<VECTOR<double,2> >::Register<FREE_PARTICLES<VECTOR<double,2> > >();
    STRUCTURE_REGISTRY<VECTOR<float,3> >::Register<FREE_PARTICLES<VECTOR<float,3> > >();
    STRUCTURE_REGISTRY<VECTOR<float,2> >::Register<FREE_PARTICLES<VECTOR<float,2> > >();
    return true;
}
bool registered_solids_structures_asdf=Register_Solids_Structures();
//#####################################################################
// Constructor
//#####################################################################
template<class TV> EMBEDDING<TV>::
EMBEDDING(GEOMETRY_PARTICLES<TV>& particles_input)
    :particles(particles_input),material_surface(material_surface_mesh,particles),need_destroy_particles(false)
{
    PHYSBAM_ASSERT(registered_solids_structures_asdf);
}
//#####################################################################
template class EMBEDDING<VECTOR<float,2> >;
template class EMBEDDING<VECTOR<float,3> >;
template class EMBEDDING<VECTOR<double,2> >;
template class EMBEDDING<VECTOR<double,3> >;
}
