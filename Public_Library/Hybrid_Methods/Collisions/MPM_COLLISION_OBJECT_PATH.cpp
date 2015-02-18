//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT_PATH.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT_PATH<TV>::
~MPM_COLLISION_OBJECT_PATH()
{
}
template class MPM_COLLISION_OBJECT_PATH<VECTOR<float,1> >;
template class MPM_COLLISION_OBJECT_PATH<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT_PATH<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT_PATH<VECTOR<double,1> >;
template class MPM_COLLISION_OBJECT_PATH<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT_PATH<VECTOR<double,3> >;
}
