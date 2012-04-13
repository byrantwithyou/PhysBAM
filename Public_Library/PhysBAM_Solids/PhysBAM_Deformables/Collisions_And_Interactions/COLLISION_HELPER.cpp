//#####################################################################
// Copyright 2012.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COLLISION_HELPER.h>
namespace PhysBAM{
//#####################################################################
// Function Compute_Collision_Impulse
//#####################################################################
template<class TV,class T,int d> TV Compute_Collision_Impulse(const TV& normal,const SYMMETRIC_MATRIX<T,d>& impulse_factor,
    const TV& relative_velocity,const T coefficient_of_restitution,const T coefficient_of_friction,bool* applied_sticking_impulse)
{
    T relative_normal_velocity=TV::Dot_Product(relative_velocity,normal);
    if(relative_normal_velocity>=0) return TV();

    // frictionless case
    if(applied_sticking_impulse) *applied_sticking_impulse=false;
    if(!coefficient_of_friction)
        return -(1+coefficient_of_restitution)*relative_normal_velocity/normal.Dot(impulse_factor*normal)*normal;

    // friction case
    // see if friction stops sliding
    TV sticking_impulse=impulse_factor.Solve_Linear_System(-coefficient_of_restitution*relative_normal_velocity*normal-relative_velocity);
    T normal_component=sticking_impulse.Dot(normal);
    if((sticking_impulse-normal_component*normal).Magnitude()<=coefficient_of_friction*normal_component){
        if(applied_sticking_impulse) *applied_sticking_impulse=true;
        return sticking_impulse;}

    // friction does not stop sliding
    TV relative_tangential_velocity=relative_velocity-relative_normal_velocity*normal;
    TV tangential_direction=relative_tangential_velocity.Normalized();
    TV impulse_direction=normal-coefficient_of_friction*tangential_direction;
    return -(1+coefficient_of_restitution)*relative_normal_velocity/normal.Dot(impulse_factor*impulse_direction)*impulse_direction;
}
}
