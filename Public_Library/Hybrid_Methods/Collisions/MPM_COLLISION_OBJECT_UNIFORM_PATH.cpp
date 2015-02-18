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
// Function Translation
//#####################################################################
template<class TV> TV MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Translation(const T t) const
{
    return frame.t;
}
//#####################################################################
// Function Translation_Velocity
//#####################################################################
template<class TV> TV MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Translation_Velocity(const T t) const
{
    return twist.linear;
}
//#####################################################################
// Function Rotation
//#####################################################################
template<class TV> ROTATION<TV> MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Rotation(const T t) const
{
    return frame.r;
}
//#####################################################################
// Function Angular_Velocity
//#####################################################################
template<class TV> typename TV::SPIN MPM_COLLISION_OBJECT_UNIFORM_PATH<TV>::
Angular_Velocity(const T t) const
{
    return twist.angular;
}
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<float,1> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<double,1> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT_UNIFORM_PATH<VECTOR<double,3> >;
}
