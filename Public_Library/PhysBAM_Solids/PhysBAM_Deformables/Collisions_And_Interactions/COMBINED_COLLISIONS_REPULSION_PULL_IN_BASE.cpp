//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_DEFORMABLE_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV,int d> COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE(TRIANGLE_COLLISIONS<TV>& triangle_collisions_input,TRIANGLE_REPULSIONS<TV>& triangle_repulsions_input,const bool update_swept_hierarchies_input)
    :triangle_collisions(triangle_collisions_input),triangle_repulsions(triangle_repulsions_input),update_swept_hierarchies(update_swept_hierarchies_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV,int d> COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
~COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE()
{
}
//#####################################################################
// Function Count
//#####################################################################
template<class TV,int d> int COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Count() const
{
    return collisions.m;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV,int d> void COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Precompute(int e)
{
}
//#####################################################################
// Function Affected_Bodies
//#####################################################################
template<class TV,int d> void COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const
{
    const COLLISION& col=collisions(e);
    for(int i=1;i<=col.nodes.m;i++)
        hash.Set(Deformable_Particle_To_Combined_Body_Id(col.nodes(i)));
}
//#####################################################################
// Function Compute_Impulse
//#####################################################################
template<class TV,int d> TV COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Compute_Impulse(int e,const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di,T scale,T dt,T time,bool diff) const
{
    const COLLISION& col=collisions(e);
    TV rel_vel;
    if(!diff){
        for(int i=1;i<=col.nodes.m;i++)
            rel_vel+=col.weights(i)*di.Particle_Velocity(col.nodes(i));}
    else{
        ARRAY_VIEW<const TV> V=triangle_collisions.geometry.deformable_body_collection.particles.V;
        for(int i=1;i<=col.nodes.m;i++)
            rel_vel+=col.weights(i)*V(col.nodes(i));}

    T relative_normal_velocity=TV::Dot_Product(rel_vel,col.normal);
    T effective_mass=Effective_Mass(e);
    
    TV j=-effective_mass*(relative_normal_velocity-col.target_relative_velocity)*col.normal;
    if(!diff)
        col.total_impulse+=j;
    return j;
}
//#####################################################################
// Function Accumulate_Impulse
//#####################################################################
template<class TV,int d> bool COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const
{
    const COLLISION& col=collisions(e);
    COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>& di=dynamic_cast<COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>&>(impulse);
    TV j=Compute_Impulse(e,di,scale,dt,time,false);

    for(int i=1;i<=col.nodes.m;i++)
        di.impulse(col.nodes(i))+=col.weights(i)*j;
    if(col.nodes(1)==4){
        TV rel_vel;
        for(int i=1;i<=col.nodes.m;i++)
            rel_vel+=col.weights(i)*di.Particle_Velocity(col.nodes(i));
        T relative_normal_velocity=TV::Dot_Product(rel_vel,col.normal);
        LOG::cout<<" measured relative velocity now "<<relative_normal_velocity<<std::endl;}

    return true;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV,int d> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const
{
    const COLLISION& col=collisions(e);

    if(TV::Dot_Product(col.total_impulse,col.normal)<=0){ // turn this off
        PHYSBAM_FATAL_ERROR("Review the following line, and if correct, use \"e!=0\"");
        // flagged_for_removal has type ARRAY<bool>, so I don't think the following line is correct.
        //flagged_for_removal.Append(e);
    }
    return 0;
}
//#####################################################################
// Function Moved
//#####################################################################
template<class TV,int d> void COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Moved(int e)
{

}
//#####################################################################
// Function Effective_Mass
//#####################################################################
template<class TV,int d> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Effective_Mass(int e) const
{
    PARTICLES<TV>& particles=triangle_collisions.geometry.deformable_body_collection.particles;
    const COLLISION& col=collisions(e);
    T impulse_factor=0;
    for(int i=1;i<=col.weights.m;i++)
        impulse_factor+=sqr(col.weights(i))*particles.one_over_mass(col.nodes(i));
    return 1/impulse_factor;
}
//#####################################################################
// Function Kinetic_Energy_Gradient
//#####################################################################
template<class TV,int d> typename TV::SCALAR COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<TV,d>::
Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const
{
    PHYSBAM_NOT_IMPLEMENTED();
    return 0;
}
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<float,1>,0>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<float,1>,2>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<float,2>,2>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<float,2>,3>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<float,3>,4>;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<double,1>,0>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<double,1>,2>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<double,2>,2>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<double,2>,3>;
template class COMBINED_COLLISIONS_REPULSION_PULL_IN_BASE<VECTOR<double,3>,4>;
#endif
