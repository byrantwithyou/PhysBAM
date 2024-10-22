//#####################################################################
// Copyright 2004-2008, Ron Fedkiw, Nipun Kwatra, Frank Losasso, Andrew Selle
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/VECTOR.h>
#include <Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <Dynamics/Particles/PARTICLE_LEVELSET_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
PARTICLE_LEVELSET_PARTICLES()
    :quantized_collision_distance(0,0),age(0,0),radius(0,0),next(0)
{
    Add_Array("quantized_collision_distance",&quantized_collision_distance);
    Add_Array("age",&age);
    Add_Array("radius",&radius);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_PARTICLES<TV>::
~PARTICLE_LEVELSET_PARTICLES()
{}
//#####################################################################
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<float,3> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_PARTICLES<VECTOR<double,3> >;
}
