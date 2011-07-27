//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_HASHTABLE.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Tools/Read_Write/Vectors/READ_WRITE_VECTOR.h>
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
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
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/COMBINED_COLLISIONS_COUPLED.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED<TV>::
COMBINED_COLLISIONS_COUPLED(RIGID_DEFORMABLE_COLLISIONS<TV>& rigid_deformable_collisions_input,ARRAY_VIEW<const TV> X_save_input,bool do_collision,bool use_saved_pairs_input)
    :rigid_deformable_collisions(rigid_deformable_collisions_input),collision_body_thickness(0),use_particles_collided_with_rigid_body(do_collision),
    use_particles_contacting_rigid_body(!do_collision),clamp_friction_magnitude(do_collision),use_collisions(do_collision),use_saved_pairs(use_saved_pairs_input),X_save(X_save_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_COUPLED<TV>::
~COMBINED_COLLISIONS_COUPLED()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Discover(const T dt,const T time)
{
    if(use_collisions) Collision_Discover();
    else Contact_Discover();
    LOG::cout<<"PROCESSING: "<<collisions.m<<std::endl;
}
//#####################################################################
// Function Contact_Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Contact_Discover()
{
    collisions.Remove_All();
    HASHTABLE<int,HASHTABLE<int> >& hash=use_saved_pairs?rigid_deformable_collisions.particles_contacting_rigid_body:rigid_deformable_collisions.particles_collided_with_rigid_body;
    for(HASHTABLE<int,HASHTABLE<int> >::ITERATOR i(hash);i.Valid();i.Next()){int b=i.Key();
        rigid_deformable_collisions.rigid_body_collisions.collision_callbacks.Swap_State(b);
        for(HASHTABLE<int>::ITERATOR j(i.Data());j.Valid();j.Next())
            Discover_Helper(b,j.Key(),X_save,COLLISION_GEOMETRY_ID());
        rigid_deformable_collisions.rigid_body_collisions.collision_callbacks.Swap_State(b);}
}
//#####################################################################
// Function Collision_Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Collision_Discover()
{
    collisions.Remove_All();
    rigid_deformable_collisions.solid_body_collection.collision_body_list.spatial_partition->Reinitialize();
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=rigid_deformable_collisions.solid_body_collection.deformable_body_collection;
    collisions.Remove_All();
    deformable_body_collection.collisions.Compute_Candidate_Nodes_For_Collision_Body_Collisions(rigid_deformable_collisions.rigid_collision_bodies);
    for(COLLISION_GEOMETRY_ID collision_body_id(1);collision_body_id<=rigid_deformable_collisions.rigid_collision_bodies.m;collision_body_id++){
        RIGID_BODY<TV>& rigid_body=dynamic_cast<RIGID_BODY<TV>&>(dynamic_cast<RIGID_COLLISION_GEOMETRY<TV>&>(*rigid_deformable_collisions.rigid_collision_bodies(collision_body_id)).rigid_geometry);
        for(int j=1;j<=deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id).m;j++)
            Discover_Helper(rigid_body.particle_index,deformable_body_collection.collisions.collision_body_candidate_nodes(collision_body_id)(j),
                rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles.X,collision_body_id);}
}
//#####################################################################
// Function Discover_Helper
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Discover_Helper(int b,int p,ARRAY_VIEW<const TV> X,COLLISION_GEOMETRY_ID collision_body_id)
{
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(b);
    COLLISION col;
    col.body=b;
    col.p=p;
    if(!use_saved_pairs && !rigid_body.Implicit_Geometry_Lazy_Inside_And_Value(particles.X(p),col.depth,collision_body_thickness)) return;
    col.depth-=collision_body_thickness;
    if(use_particles_collided_with_rigid_body)
        rigid_deformable_collisions.particles_collided_with_rigid_body.Get_Or_Insert(b).Set(col.p);
    if(!use_saved_pairs && use_particles_contacting_rigid_body)
        rigid_deformable_collisions.particles_contacting_rigid_body.Get_Or_Insert(b).Set(col.p);
    col.normal=rigid_body.Implicit_Geometry_Normal(X(col.p));
    col.relative_velocity=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles.V(col.p)-rigid_body.Pointwise_Object_Velocity(X(col.p));
    if(!use_saved_pairs && TV::Dot_Product(col.relative_velocity,col.normal)>=0) return;
    col.collision_body_id=collision_body_id;
    collisions.Append(col);
}
//#####################################################################
// Function Count
//#####################################################################
template<class TV> int COMBINED_COLLISIONS_COUPLED<TV>::
Count() const
{
    return collisions.m;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Precompute(int e)
{
}
//#####################################################################
// Function Affected_Bodies
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const
{
    const COLLISION& col=collisions(e);
    hash.Set(Rigid_Body_To_Combined_Body_Id(col.body));
    hash.Set(Deformable_Particle_To_Combined_Body_Id(col.p));
}
//#####################################################################
// Function Accumulate_Impulse
//#####################################################################
template<class TV> bool COMBINED_COLLISIONS_COUPLED<TV>::
Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return true;
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    const COLLISION& col=collisions(e);
    RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(col.body);
    TV collision_impulse;
    rigid_deformable_collisions.Apply_Rigid_Deformable_Collision_Impulse(rigid_body,col.p,use_collisions?particles.X(col.p):X_save(col.p),col.normal,col.relative_velocity,
        scale-1,rigid_body.coefficient_of_friction,clamp_friction_magnitude,collision_impulse,use_saved_pairs,false);

    dynamic_cast<COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>&>(impulse).Apply_Coupled_Impulse(col.p,col.body,collision_impulse);
    return true;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_COUPLED<TV>::
Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const
{
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    const COLLISION& col=collisions(e);
    const COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>& coupled_impulse=dynamic_cast<const COMBINED_COLLISIONS_COUPLED_IMPULSE<TV>&>(impulse);
    TV rel_vel=coupled_impulse.Particle_Velocity(col.p)-coupled_impulse.Velocity_At_Point(col.body,particles.X(col.p));
    T orig_rel_vel_N=TV::Dot_Product(col.relative_velocity,col.normal);
    if(!orig_rel_vel_N) return 0;
    T rel_vel_N=TV::Dot_Product(rel_vel,col.normal);
    return 1-rel_vel_N/orig_rel_vel_N;
}
//#####################################################################
// Function Moved
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_COUPLED<TV>::
Moved(int e)
{
    const COLLISION& col=collisions(e);
    if(rigid_deformable_collisions.solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects && col.collision_body_id)
        rigid_deformable_collisions.solid_body_collection.collision_body_list.spatial_partition->Update_Body(col.collision_body_id);
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_COUPLED<TV>::
Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return 0;
    PARTICLES<TV>& particles=rigid_deformable_collisions.solid_body_collection.deformable_body_collection.particles;
    const COLLISION& col=collisions(e);
    RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(col.body);
    TV j=Compute_Impulse_Derivative(e,scale,0,rigid_body.coefficient_of_friction,clamp_friction_magnitude);
    T rj=dynamic_cast<const COMBINED_COLLISIONS_RIGID_IMPULSE<TV>&>(impulse).Kinetic_Energy_Change(col.body,particles.X(col.p),TWIST<TV>(-j,typename TV::SPIN()));
    T dj=dynamic_cast<const COMBINED_COLLISIONS_DEFORMABLE_IMPULSE<TV>&>(impulse).Kinetic_Energy_Change(col.p,j);
    return rj+dj;
}
//#####################################################################
// Function Compute_Impulse_Derivative
//#####################################################################
template<class TV> TV COMBINED_COLLISIONS_COUPLED<TV>::
Compute_Impulse_Derivative(int e,T scale,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude) const 
{
    const COLLISION& col=collisions(e);
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=rigid_deformable_collisions.solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY<TV>& rigid_body=rigid_deformable_collisions.solid_body_collection.rigid_body_collection.Rigid_Body(col.body);

    T adjusted_cor=(coefficient_of_restitution+1)*scale-1;
    T d_adjusted_cor=coefficient_of_restitution+1;
    TV impulse;
    TV d_impulse;
    bool allow_pull=use_saved_pairs;

    T relative_normal_velocity=TV::Dot_Product(col.relative_velocity,col.normal);
    if(relative_normal_velocity>0 && !allow_pull) relative_normal_velocity=0;

    if(!coefficient_of_friction || (allow_pull && relative_normal_velocity>=0)){ // frictionless case
        T nT_impulse1_n=particles.one_over_mass(col.p),nT_impulse2_n=0;
        if(!rigid_body.Has_Infinite_Inertia()){
            T_SPIN r2xn=TV::Cross_Product(particles.X(col.p)-rigid_body.X(),col.normal);
            nT_impulse2_n=1/rigid_body.Mass()+TV::SPIN::Dot_Product(r2xn,rigid_body.World_Space_Inertia_Tensor_Inverse()*r2xn);}
        impulse=-(1+adjusted_cor)*relative_normal_velocity/(nT_impulse1_n+nT_impulse2_n)*col.normal;
        d_impulse=-d_adjusted_cor*relative_normal_velocity/(nT_impulse1_n+nT_impulse2_n)*col.normal;}
    else{    
        typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX impulse_factor=rigid_body.Impulse_Factor(particles.X(col.p))+particles.one_over_mass(col.p),impulse_factor_inverse=impulse_factor.Inverse();
        // see if friction stops sliding
        impulse=-impulse_factor_inverse*(adjusted_cor*relative_normal_velocity*col.normal+col.relative_velocity); // sticking impulse
        d_impulse=-impulse_factor_inverse*(d_adjusted_cor*relative_normal_velocity*col.normal); // sticking impulse
        T normal_component=TV::Dot_Product(impulse,col.normal);
        if((impulse-normal_component*col.normal).Magnitude()>coefficient_of_friction*normal_component){ // sticking impulse was not admissible
            // friction does not stop sliding
            TV relative_tangential_velocity=col.relative_velocity.Projected_Orthogonal_To_Unit_Direction(col.normal);
            TV tangential_direction=relative_tangential_velocity;
            T relative_tangential_velocity_magnitude=tangential_direction.Normalize();
            TV impulse_factor_times_direction=impulse_factor*(col.normal-coefficient_of_friction*tangential_direction);
            assert(TV::Dot_Product(impulse_factor_times_direction,col.normal));
            TV delta=-(1+adjusted_cor)*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,col.normal)*impulse_factor_times_direction;
            TV d_delta=-d_adjusted_cor*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,col.normal)*impulse_factor_times_direction;
            // clamp friction magnitude: should be (true) in the elastic collision case, (false) in the inelastic collision case
            if(clamp_friction_magnitude){ // should only clamp friction magnitude in the elastic case!
                TV new_relative_velocity=col.relative_velocity+delta;
                TV d_new_relative_velocity=col.relative_velocity+d_delta;
                TV new_normal_velocity=new_relative_velocity.Projected_On_Unit_Direction(col.normal);
                TV d_new_normal_velocity=d_new_relative_velocity.Projected_On_Unit_Direction(col.normal);
                TV new_tangential_velocity=new_relative_velocity-new_normal_velocity;
                TV d_new_tangential_velocity=d_new_relative_velocity-d_new_normal_velocity;
                T new_tangential_velocity_magnitude=new_tangential_velocity.Magnitude();
                T d_new_tangential_velocity_magnitude=TV::Dot_Product(new_tangential_velocity,d_new_tangential_velocity)/new_tangential_velocity_magnitude;
                if(new_tangential_velocity_magnitude>relative_tangential_velocity_magnitude){
                    delta=new_normal_velocity+(relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*new_tangential_velocity-col.relative_velocity;
                    d_delta=d_new_normal_velocity+(-relative_tangential_velocity_magnitude/sqr(new_tangential_velocity_magnitude)*d_new_tangential_velocity_magnitude)*new_tangential_velocity
                        +(relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*d_new_tangential_velocity;}}
            impulse=impulse_factor_inverse*delta;
            d_impulse=impulse_factor_inverse*d_delta;}}
        
//    T ds=(T)1e-6;
//    TV j,jf,jb;
//    rigid_deformable_collisions.Apply_Rigid_Deformable_Collision_Impulse(rigid_body,col.p,particles.X(col.p),col.normal,col.relative_velocity,scale-1,rigid_body.coefficient_of_friction,
//        clamp_friction_magnitude,j,allow_pull,false);
//    rigid_deformable_collisions.Apply_Rigid_Deformable_Collision_Impulse(rigid_body,col.p,particles.X(col.p),col.normal,col.relative_velocity,scale+ds-1,rigid_body.coefficient_of_friction,
//        clamp_friction_magnitude,jf,allow_pull,false);
//    rigid_deformable_collisions.Apply_Rigid_Deformable_Collision_Impulse(rigid_body,col.p,particles.X(col.p),col.normal,col.relative_velocity,scale-ds-1,rigid_body.coefficient_of_friction,
//        clamp_friction_magnitude,jb,allow_pull,false);


    return d_impulse;
}
template class COMBINED_COLLISIONS_COUPLED<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_COUPLED<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_COUPLED<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_COUPLED<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_COUPLED<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_COUPLED<VECTOR<double,3> >;
#endif
