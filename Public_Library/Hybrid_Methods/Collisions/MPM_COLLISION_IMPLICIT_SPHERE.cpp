//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/FRAME.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Vectors/TWIST.h>
#include <Geometry/Basic_Geometry/SPHERE.h>
#include <Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_INVERT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_IMPLICIT_SPHERE.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_COLLISION_IMPLICIT_SPHERE<TV>::
MPM_COLLISION_IMPLICIT_SPHERE(COLLISION_TYPE type_input,T friction_input,
    std::function<FRAME<TV>(T)> func_frame,std::function<TWIST<TV>(T)> func_twist,TV center,T radius,T shrinking_speed)
    :MPM_COLLISION_IMPLICIT_OBJECT<TV>(0,type_input,friction_input,func_frame,func_twist),sphere(new ANALYTIC_IMPLICIT_OBJECT<SPHERE<TV> >(SPHERE<TV>(center,radius))),radius(radius),shrinking_speed(shrinking_speed)
{
    io=new IMPLICIT_OBJECT_INVERT<TV>(sphere);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_IMPLICIT_SPHERE<TV>::
~MPM_COLLISION_IMPLICIT_SPHERE()
{
}
//#####################################################################
// Function Phi
//#####################################################################
template<class TV> typename TV::SCALAR MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Phi(const TV& X,T time) const
{
    Update(time);
    return BASE::Phi(X,time);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Normal(const TV& X,T time) const
{
    Update(time);
    return BASE::Normal(X,time);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Hessian(const TV& X,T time) const
{
    Update(time);
    return BASE::Hessian(X,time);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Velocity(const TV& X,T time) const
{
    TV n=Normal(X,time);
    return shrinking_speed*n;
}
//#####################################################################
// Function Get_Implicit_Object
//#####################################################################
template<class TV> IMPLICIT_OBJECT<TV>* MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Get_Implicit_Object(T time)
{
    Update(time);
    return BASE::Get_Implicit_Object(time);
}
//#####################################################################
// Function Update
//#####################################################################
template<class TV> void MPM_COLLISION_IMPLICIT_SPHERE<TV>::
Update(T time) const
{
    sphere->analytic.radius=radius-time*shrinking_speed;
}
//#####################################################################
template class MPM_COLLISION_IMPLICIT_SPHERE<VECTOR<float,2> >;
template class MPM_COLLISION_IMPLICIT_SPHERE<VECTOR<float,3> >;
template class MPM_COLLISION_IMPLICIT_SPHERE<VECTOR<double,2> >;
template class MPM_COLLISION_IMPLICIT_SPHERE<VECTOR<double,3> >;
}
