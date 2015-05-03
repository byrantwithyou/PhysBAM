//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT<TV>::
~MPM_COLLISION_OBJECT()
{
}
template class MPM_COLLISION_OBJECT<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,3> >;
}
