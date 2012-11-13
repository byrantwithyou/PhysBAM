//#####################################################################
// Copyright 2004, Ron Fedkiw, Frank Losasso, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Vectors/VECTOR.h>
#include <PhysBAM_Dynamics/Particles/DYNAMICS_PARTICLES_FORWARD.h>
#include <PhysBAM_Dynamics/Particles/PARTICLE_LEVELSET_REMOVED_PARTICLES.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>::
PARTICLE_LEVELSET_REMOVED_PARTICLES()
    :V(0,0),next(0)
{
    this->Add_Array(ATTRIBUTE_ID_V,&V);
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> PARTICLE_LEVELSET_REMOVED_PARTICLES<TV>::
~PARTICLE_LEVELSET_REMOVED_PARTICLES()
{}
//#####################################################################
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,1> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,2> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<float,3> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,1> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,2> >;
template class PARTICLE_LEVELSET_REMOVED_PARTICLES<VECTOR<double,3> >;
}
