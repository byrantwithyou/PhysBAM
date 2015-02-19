//#####################################################################
// Copyright 2015, Chenfanfu Jiang.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/FRAME.h>
#include <Tools/Vectors/TWIST.h>
#include <Tools/Vectors/VECTOR.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT.h>
#include <Hybrid_Methods/Collisions/MPM_COLLISION_OBJECT_PATH.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT<TV>::
MPM_COLLISION_OBJECT(IMPLICIT_OBJECT<TV>* impob,MPM_COLLISION_OBJECT_PATH<TV>* path,const bool sticky,const T friction)
    :impob(impob),path(path),sticky(sticky),friction(friction)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_COLLISION_OBJECT<TV>::
~MPM_COLLISION_OBJECT()
{
    if(impob) delete impob;
    if(path) delete path;
}
//#####################################################################
// Function Detect
//#####################################################################
template<class TV> bool MPM_COLLISION_OBJECT<TV>::
Detect(const T t,const TV& x,T* phi,TV* n) const
{
    TV X=path->Orientation(t).Inverse_Times(x);
    T ep=impob->Extended_Phi(X);
    if(ep<0){
        if(phi) *phi=ep;
        if(n) *n=path->Orientation(t).r.Rotate(impob->Normal(X));
        return true;}
    return false;
}
//#####################################################################
// Function Collide
//#####################################################################
template<class TV> bool MPM_COLLISION_OBJECT<TV>::
Collide(const T t,const TV& x,TV& v,T* phi,TV* n,bool apply_friction) const
{
    TV normal;
    if(Detect(t,x,phi,&normal)){
        TV V=TV::Cross_Product(path->Velocity(t).angular,x-path->Orientation(t).t)+path->Velocity(t).linear;
        v-=V;
        Collide_Static(t,x,normal,v,apply_friction);
        v+=V;
        if(n) *n=normal;
        return true;}
    return false;
}
//#####################################################################
// Function Collide_Static
//#####################################################################
template<class TV> void MPM_COLLISION_OBJECT<TV>::
Collide_Static(const T t,const TV& x,const TV& n,TV& v,bool apply_friction) const
{
    if(sticky){v=TV();return;}
    T vn=TV::Dot_Product(v,n);
    if(vn<0){
        v-=n*vn;
        if(apply_friction){
            if(-vn*friction<v.Magnitude())
                v+=v.Normalized()*vn*friction;
            else v=TV();}}
}
template class MPM_COLLISION_OBJECT<VECTOR<float,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<float,3> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,2> >;
template class MPM_COLLISION_OBJECT<VECTOR<double,3> >;
}
