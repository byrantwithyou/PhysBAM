//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Matrices/SYMMETRIC_MATRIX.h>
#include <Rigids/Collisions/COLLISION_HELPER.h>
using namespace PhysBAM;
//#####################################################################
// Function Compute_Collision_Impulse
//#####################################################################
template<class TV,class T,int d> TV PhysBAM::Compute_Collision_Impulse(const TV& normal,const SYMMETRIC_MATRIX<T,d>& impulse_factor,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,bool* applied_sticking_impulse)
{
    if(!impulse_factor.Determinant()) return TV();
    T relative_normal_velocity=TV::Dot_Product(relative_velocity,normal);
    if(relative_normal_velocity>=0) return TV();

    // frictionless case
    if(applied_sticking_impulse) *applied_sticking_impulse=false;
    if(!coefficient_of_friction)
        return -(1+coefficient_of_restitution)*relative_normal_velocity/normal.Dot(impulse_factor*normal)*normal;

    // friction case
    // see if friction stops sliding
    TV sticking_acceleration=-coefficient_of_restitution*relative_normal_velocity*normal-relative_velocity,sticking_impulse=impulse_factor.Inverse_Times(sticking_acceleration);
    T normal_component=sticking_impulse.Dot(normal);
    if((sticking_impulse-normal_component*normal).Magnitude()<=coefficient_of_friction*normal_component){
        if(sticking_impulse.Dot(relative_velocity+(T).5*sticking_acceleration)<=0){
            if(applied_sticking_impulse) *applied_sticking_impulse=true;
            return sticking_impulse;}}

    // friction does not stop sliding
    TV relative_tangential_velocity=relative_velocity-relative_normal_velocity*normal;
    TV tangential_direction=relative_tangential_velocity.Normalized();
    TV impulse_direction=normal-coefficient_of_friction*tangential_direction;
    TV impulse=-(1+coefficient_of_restitution)*relative_normal_velocity/normal.Dot(impulse_factor*impulse_direction)*impulse_direction;
    T a=impulse.Dot(impulse_factor*impulse),b=impulse.Dot(relative_velocity);
    if(a+2*b<=0) return impulse;
    if(b>=0) return TV();
    return (-2*b/a)*impulse;
}
namespace PhysBAM{
template VECTOR<float,1> Compute_Collision_Impulse<VECTOR<float,1>,float,1>(VECTOR<float,1> const&,SYMMETRIC_MATRIX<float,1> const&,VECTOR<float,1> const&,float,float,bool*);
template VECTOR<float,2> Compute_Collision_Impulse<VECTOR<float,2>,float,2>(VECTOR<float,2> const&,SYMMETRIC_MATRIX<float,2> const&,VECTOR<float,2> const&,float,float,bool*);
template VECTOR<float,3> Compute_Collision_Impulse<VECTOR<float,3>,float,3>(VECTOR<float,3> const&,SYMMETRIC_MATRIX<float,3> const&,VECTOR<float,3> const&,float,float,bool*);
template VECTOR<double,1> Compute_Collision_Impulse<VECTOR<double,1>,double,1>(VECTOR<double,1> const&,SYMMETRIC_MATRIX<double,1> const&,VECTOR<double,1> const&,double,double,bool*);
template VECTOR<double,2> Compute_Collision_Impulse<VECTOR<double,2>,double,2>(VECTOR<double,2> const&,SYMMETRIC_MATRIX<double,2> const&,VECTOR<double,2> const&,double,double,bool*);
template VECTOR<double,3> Compute_Collision_Impulse<VECTOR<double,3>,double,3>(VECTOR<double,3> const&,SYMMETRIC_MATRIX<double,3> const&,VECTOR<double,3> const&,double,double,bool*);
}
