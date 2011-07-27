//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_REPULSION_BASE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_REPULSION_BASE<TV>::
COMBINED_COLLISIONS_REPULSION_BASE(TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,bool prune_pairs_input)
    :triangle_repulsions(triangle_repulsions_input),prune_pairs(prune_pairs_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_REPULSION_BASE<TV>::
~COMBINED_COLLISIONS_REPULSION_BASE()
{
}
//#####################################################################
// Function Count
//#####################################################################
template<class TV> int COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Count() const
{
    return collisions.m;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Precompute(int e)
{
}
//#####################################################################
// Function Affected_Bodies
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const
{
    const COLLISION& col=collisions(e);
    for(int i=1;i<=col.nodes.m;i++)
        hash.Set(Deformable_Particle_To_Combined_Body_Id(col.nodes(i)));
}
//#####################################################################
// Function Compute_Impulse
//#####################################################################
template<class TV> TV COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Compute_Impulse(int e,T scale,T dt,T time,bool diff) const
{
    const COLLISION& col=collisions(e);
    T relative_normal_velocity=TV::Dot_Product(col.relative_velocity,col.normal);
    TV tangential_velocity=col.relative_velocity-relative_normal_velocity*col.normal;
    if(relative_normal_velocity>0) relative_normal_velocity=0;
    TV tangential_direction=tangential_velocity;
    T tangential_velocity_magnitude=tangential_direction.Normalize();
    T effective_mass=Effective_Mass(e);
    // see if friction stops sliding
    if(tangential_velocity_magnitude<=-triangle_repulsions.friction_coefficient*scale*relative_normal_velocity){
        if(diff) return -effective_mass*relative_normal_velocity*col.normal;
        else return -effective_mass*(tangential_velocity+scale*relative_normal_velocity*col.normal);}

    TV j=-effective_mass*relative_normal_velocity*(col.normal-triangle_repulsions.friction_coefficient*tangential_direction);
    if(!diff) j*=scale;
    return j;
}
//#####################################################################
// Function Accumulate_Impulse
//#####################################################################
template<class TV> bool COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return true;
    const COLLISION& col=collisions(e);
    TV j=Compute_Impulse(e,scale,dt,time,false);

    COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di=dynamic_cast<COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>&>(impulse);
    for(int i=1;i<=col.nodes.m;i++)
        di.impulse(col.nodes(i))+=col.weights(i)*j;
    return true;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const
{
    const COLLISION& col=collisions(e);
    ARRAY_VIEW<const TV> V=triangle_repulsions.geometry.deformable_body_collection.particles.V;

    const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di=dynamic_cast<const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>&>(impulse);
    TV rel_vel;
    for(int i=1;i<=col.nodes.m;i++)
        rel_vel+=col.weights(i)*di.Particle_Velocity(col.nodes(i));
    T orig_rel_vel_N=TV::Dot_Product(col.relative_velocity,col.normal);
    if(!orig_rel_vel_N) return 0;
    T rel_vel_N=TV::Dot_Product(rel_vel,col.normal);
    return 1-rel_vel_N/orig_rel_vel_N;
}
//#####################################################################
// Function Moved
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Moved(int e)
{
}
//#####################################################################
// Function Effective_Mass
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Effective_Mass(int e) const
{
    PARTICLES<TV>& particles=triangle_repulsions.geometry.deformable_body_collection.particles;
    const COLLISION& col=collisions(e);
    T impulse_factor=0;
    for(int i=1;i<=col.weights.m;i++)
        impulse_factor+=sqr(col.weights(i))*particles.one_over_mass(col.nodes(i));
    return 1/impulse_factor;
}
//#####################################################################
// Function Kinetic_Energy_Gradient
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_BASE<TV>::
Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return 0;
    const COLLISION& col=collisions(e);
    TV j=Compute_Impulse(e,scale,dt,time,true);
    T ke=0;
    const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di=dynamic_cast<const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>&>(impulse);
    for(int i=1;i<=col.nodes.m;i++)
        ke+=di.Kinetic_Energy_Change(col.nodes(i),col.weights(i)*j);
    return ke;
}
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_REPULSION_BASE<VECTOR<double,3> >;
#endif
