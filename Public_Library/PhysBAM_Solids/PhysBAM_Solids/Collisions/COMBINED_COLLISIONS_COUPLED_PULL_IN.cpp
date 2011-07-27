//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Tools/Vectors/Dot_Product.h>
#include <PhysBAM_Geometry/Collision_Detection/COLLISION_GEOMETRY_SPATIAL_PARTITION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_1D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_2D.h>
#include <PhysBAM_Geometry/Collisions/RIGID_COLLISION_GEOMETRY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/COMBINED_COLLISIONS_COUPLED_PULL_IN.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
COMBINED_COLLISIONS_COUPLED_PULL_IN(RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions_input,const T dt_input)
    :rigid_deformable_collisions(rigid_deformable_collisions_input),collision_body_thickness(rigid_deformable_collisions_input.solids_parameters.rigid_body_collision_parameters.collision_body_thickness),dt(dt_input)
{
    // clear the contact list
    rigid_deformable_collisions.particles_contacting_rigid_body.Remove_All();
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
~COMBINED_COLLISIONS_COUPLED_PULL_IN()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Discover(const T dt,const T time)
{
    for(int i=flagged_for_removal.m;i>0;i--){
        int index=flagged_for_removal(i);
        rigid_deformable_collisions.particles_contacting_rigid_body.Get_Or_Insert(collisions(index).body).Delete(collisions(index).p);
        collisions.Remove_Index_Lazy(index);}
    flagged_for_removal.Remove_All();
        
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=rigid_deformable_collisions.solid_body_collection.deformable_body_collection;
    deformable_body_collection.collisions.Compute_Candidate_Nodes_For_Collision_Body_Collisions(rigid_deformable_collisions.rigid_collision_bodies);

    for(COLLISION_GEOMETRY_ID collision_body_id(1);collision_body_id<=rigid_deformable_collisions.rigid_collision_bodies.m;collision_body_id++){
        RIGID_BODY<TV>& rigid_body=dynamic_cast<RIGID_BODY<TV>&>(dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>&>(*rigid_deformable_collisions.rigid_collision_bodies(collision_body_id)).rigid_geometry);
        for(int j=1;j<=deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id).m;j++){
            int p=deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id)(j);
            if(rigid_deformable_collisions.particles_contacting_rigid_body.Get_Or_Insert(rigid_body.particle_index).Contains(p)) continue; // already got one

            T depth=rigid_body.Implicit_Geometry_Extended_Value(particles.X(p));
            if(depth>=.99*collision_body_thickness) continue;

            TV collision_normal=rigid_body.Implicit_Geometry_Normal(particles.X(p));
            TV body_V=rigid_body.Pointwise_Object_Velocity(particles.X(p));
            T relative_speed=TV::Dot_Product(particles.V(p)-body_V,collision_normal);
            T target_relative_speed=relative_speed+(-depth+collision_body_thickness)/dt;
            //if(target_relative_speed<0)
            rigid_deformable_collisions.particles_contacting_rigid_body.Get_Or_Insert(rigid_body.particle_index).Set(p);

            COLLISION col;
            col.body=rigid_body.particle_index;
            col.p=p;
            col.normal=collision_normal;
            col.relative_speed=relative_speed;
            col.target_relative_speed=target_relative_speed;
            col.total_impulse=TV();
            collisions.Append(col);}}
    LOG::cout<<"PROCESSING: "<<collisions.m<<std::endl;
}
//#####################################################################
// Function Count
//#####################################################################
template<class TV> int COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Count() const
{
    return collisions.m;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Precompute(int e)
{
}
//#####################################################################
// Function Affected_Bodies
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const
{
    const COLLISION& col=collisions(e);
    hash.Set(Rigid_Body_To_Combined_Body_Id(col.body));
    hash.Set(Deformable_Particle_To_Combined_Body_Id(col.p));
}
//#####################################################################
// Function Compute_Impulse
//#####################################################################
template<class TV> TV COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Compute_Impulse(int e,const COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>& ci,T scale,T dt,T time,bool diff) const
{
    const COLLISION& col=collisions(e);

    T effective_mass=Effective_Mass(e);
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    TV rel_vel;
    if(!diff) rel_vel=ci.Particle_Velocity(col.p)-ci.Velocity_At_Point(col.body,particles.X(col.p));
    else{
        PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
        RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(col.body);
        rel_vel=particles.V(col.p)-rigid_body.Pointwise_Object_Velocity(particles.X(col.p));}
    T rel_vel_N=TV::Dot_Product(rel_vel,col.normal);
    TV j=-effective_mass*(rel_vel_N-col.target_relative_speed)*col.normal;
    if(!diff)
        col.total_impulse+=j;
    return j;
}
//#####################################################################
// Function Accumulate_Impulse
//#####################################################################
template<class TV> bool COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const
{
    const COLLISION& col=collisions(e);

    // TODO: dynamic cast in inner loop is slow
    COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>& ci=dynamic_cast<COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>&>(impulse);
    TV j=Compute_Impulse(e,ci,scale,dt,time,false);
   
    ci.Apply_Coupled_Impulse(col.p,col.body,j);
    return true;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const
{
    const COLLISION& col=collisions(e);

    if(TV::Dot_Product(col.total_impulse,col.normal)<=0) // take it out of the hash table
        flagged_for_removal.Append(e);
    return 0;
}
//#####################################################################
// Function Effective_Mass
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Effective_Mass(int e) const
{
    const COLLISION& col=collisions(e);
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    T nT_impulse1_n=particles.one_over_mass(col.p),nT_impulse2_n=0;
    RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(col.body);
    if(!rigid_body.Has_Infinite_Inertia()){
        T_SPIN r2xn=TV::Cross_Product(particles.X(col.p)-rigid_body.X(),col.normal);
        nT_impulse2_n=1/rigid_body.Mass()+Dot_Product(r2xn,rigid_body.World_Space_Inertia_Tensor_Inverse()*r2xn);}
    return 1/(nT_impulse1_n+nT_impulse2_n);
}
//#####################################################################
// Function Moved
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED_PULL_IN<TV>::
Moved(int e)
{
}
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_COUPLED_PULL_IN<VECTOR<double,3> >;
#endif
