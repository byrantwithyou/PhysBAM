//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/FRAME.h>
#include <Tools/Vectors/TWIST.h>
#include <Tools/Vectors/VECTOR.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT_UNIFORM_PATH.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
MPM_COLLISION_OBJECT_UNIFORM_PATH(const FRAME<TV>& frame,const TWIST<TV>& twist)
    :frame(frame),twist(twist)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
~MPM_COLLISION_OBJECT_UNIFORM_PATH()
{
}
//#####################################################################
// Function Orientation
//#####################################################################
template<class TV> FRAME<TV> MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Orientation(const T t) const
{
    return frame;
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TWIST<TV> MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Velocity(const T t) const
{
    return twist;
}
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<double,3> >;
}
