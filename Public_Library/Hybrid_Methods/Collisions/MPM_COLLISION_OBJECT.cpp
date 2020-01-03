//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Vectors/TWIST.h>
#include <Core/Vectors/VECTOR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
namespace PhysBAM{
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT<TV>::
~MPM_COLLISION_OBJECT()
{
    delete iot_frame;
    delete iot;
    delete io;
}
//#####################################################################
// Function Phi
//#####################################################################
template<class TV> typename TV::SCALAR MPM_COLLISION_OBJECT<TV>::
Phi(const TV& X,T time) const
{
    if(func_frame){
        FRAME<TV> frame=func_frame(time);
        TV Z=frame.Inverse_Times(X);
        return io->Extended_Phi(Z);}
    return io->Extended_Phi(X);
}
//#####################################################################
// Function Normal
//#####################################################################
template<class TV> TV MPM_COLLISION_OBJECT<TV>::
Normal(const TV& X,T time) const
{
    if(func_frame){
        FRAME<TV> frame=func_frame(time);
        TV Z=frame.Inverse_Times(X);
        TV N=io->Extended_Normal(Z);
        return frame.r.Rotate(N);}
    return io->Extended_Normal(X);
}
//#####################################################################
// Function Hessian
//#####################################################################
template<class TV> SYMMETRIC_MATRIX<typename TV::SCALAR,TV::m> MPM_COLLISION_OBJECT<TV>::
Hessian(const TV& X,T time) const
{
    if(func_frame){
        FRAME<TV> frame=func_frame(time);
        TV Z=frame.Inverse_Times(X);
        SYMMETRIC_MATRIX<T,TV::m> H=io->Hessian(Z);
        return SYMMETRIC_MATRIX<T,TV::m>::Conjugate(frame.r.Rotation_Matrix(),H);}
    return io->Hessian(X);
}
//#####################################################################
// Function Velocity
//#####################################################################
template<class TV> TV MPM_COLLISION_OBJECT<TV>::
Velocity(const TV& X,T time) const
{
    if(func_frame){
        FRAME<TV> frame=func_frame(time);
        TWIST<TV> twist=func_twist(time);
        return twist.linear+twist.angular.Cross(X-frame.t);}
    return TV();
}
//#####################################################################
// Function Get_Implicit_Object
//#####################################################################
template<class TV> IMPLICIT_OBJECT<TV>* MPM_COLLISION_OBJECT<TV>::
Get_Implicit_Object(T time)
{
    if(func_frame){
        if(!iot){
            iot_frame=new FRAME<TV>;
            iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(io,false,iot_frame);}
        *iot_frame=func_frame(time);
        return iot;}
    return io;
}
//#####################################################################
template class MPM_COLLISION_OBJECT<VECTOR<float,1> >;
template class MPM_COLLISION_OBJECT<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,1> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,3> >;
}
