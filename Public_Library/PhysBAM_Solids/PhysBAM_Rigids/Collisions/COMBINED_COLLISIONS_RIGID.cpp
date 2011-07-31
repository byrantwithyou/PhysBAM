//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Utilities/DEBUG_CAST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/COMBINED_COLLISIONS_RIGID_IMPULSE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_PARTICLE_INTERSECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGIDS_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_RIGID<TV>::
COMBINED_COLLISIONS_RIGID(RIGID_BODY_COLLISIONS<TV>& rigid_body_collisions_input,bool do_collision,bool use_saved_pairs_input)
    :rigid_body_collisions(rigid_body_collisions_input),ignore_separating(true),rolling_friction(true),clamp_energy(true),clamp_friction_magnitude(true),desired_separation_distance(0),
    use_coefficient_of_restitution(do_collision),use_pairs_processed_by_collisions(do_collision),use_saved_pairs(use_saved_pairs_input),use_collisions(do_collision)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> COMBINED_COLLISIONS_RIGID<TV>::
~COMBINED_COLLISIONS_RIGID()
{
}
//#####################################################################
// Function Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Discover(const T dt,const T time)
{
    collisions.Remove_All();
    if(use_collisions) Collision_Discover(dt,time);
    else Contact_Discover();
}
//#####################################################################
// Function Contact_Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Contact_Discover()
{
    if(use_saved_pairs) Discover_Helper(rigid_body_collisions.saved_pairs);
    else
        for(int i=1;i<=rigid_body_collisions.precomputed_contact_pairs_for_level.m;i++)
            Discover_Helper(rigid_body_collisions.precomputed_contact_pairs_for_level(i));
}
//#####################################################################
// Function Collision_Discover
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Collision_Discover(const T dt,const T time)
{
    rigid_body_collisions.saved_pairs.Remove_All();
    ARRAY<VECTOR<int,2> > pairs;
    rigid_body_collisions.Get_Bounding_Box_Collision_Pairs(dt,time,pairs,true,true,rigid_body_collisions.parameters.collision_bounding_box_thickness);
    Discover_Helper(pairs);
}
//#####################################################################
// Function Discover_Helper
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Discover_Helper(const ARRAY<VECTOR<int,2> >& pairs)
{
    ARRAY<RIGID_BODY_PARTICLE_INTERSECTION<TV> > particle_intersections;
    for(int i=1;i<=pairs.m;i++){
        particle_intersections.Remove_All();
        rigid_body_collisions.intersections.Append_All_Intersections(pairs(i)(1),pairs(i)(2),particle_intersections,desired_separation_distance);
        for(int j=1;j<=particle_intersections.m;j++){
            const RIGID_BODY_PARTICLE_INTERSECTION<TV>& intersection=particle_intersections(j);
            RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(intersection.particle_body);
            RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(intersection.levelset_body);

            COLLISION col;
            col.levelset_body=intersection.levelset_body;
            col.particle_body=intersection.particle_body;
            col.particle_index=intersection.particle_index;
            col.location=body1.World_Space_Point(intersection.particle_location);
            col.normal=body2.Implicit_Geometry_Normal(col.location);
            col.relative_velocity=RIGID_GEOMETRY<TV>::Relative_Velocity_At_Geometry1_Particle(body1,body2,col.location,intersection.particle_index);

            if(!use_saved_pairs && ignore_separating && TV::Dot_Product(col.relative_velocity,col.normal)>=0) continue;

            if(use_pairs_processed_by_collisions)
                rigid_body_collisions.pairs_processed_by_collisions.Set(VECTOR<int,2>(col.particle_body,col.levelset_body).Sorted());
            else if(!use_saved_pairs){
                rigid_body_collisions.rigid_body_particle_intersections.Set(TRIPLE<int,int,TV>(col.particle_body,col.levelset_body,col.location));
                rigid_body_collisions.saved_pairs.Append(VECTOR<int,2>(col.particle_body,col.levelset_body));}
            collisions.Append(col);}}
}
//#####################################################################
// Function Count
//#####################################################################
template<class TV> int COMBINED_COLLISIONS_RIGID<TV>::
Count() const
{
    return collisions.m;
}
//#####################################################################
// Function Precompute
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Precompute(int e)
{
}
//#####################################################################
// Function Affected_Bodies
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Affected_Bodies(int e,HASHTABLE<COMBINED_BODY_ID>& hash) const
{
    const COLLISION& col=collisions(e);
    hash.Set(Rigid_Body_To_Combined_Body_Id(col.particle_body));
    hash.Set(Rigid_Body_To_Combined_Body_Id(col.levelset_body));
}
//#####################################################################
// Function Accumulate_Impulse
//#####################################################################
template<class TV> bool COMBINED_COLLISIONS_RIGID<TV>::
Accumulate_Impulse(int e,IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return true;
    const COLLISION& col=collisions(e);
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.particle_body);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.levelset_body);

    // TODO: handle coefficient of restitution correctly.
    T cor=0;
    if(use_coefficient_of_restitution) cor=RIGID_BODY<TV>::Coefficient_Of_Restitution(body1,body2);

    TWIST<TV> j=rigid_body_collisions.collision_callbacks.Compute_Collision_Impulse(body1,body2,col.location,col.normal,col.relative_velocity,
        (cor+1)*scale-1,RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),clamp_friction_magnitude,rolling_friction,clamp_energy);

    dynamic_cast<COMBINED_COLLISIONS_RIGID_IMPULSE<TV>&>(impulse).Apply_Rigid_Impulse(col.particle_body,col.levelset_body,col.location,j);
    return true;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_RIGID<TV>::
Diagnose_Impulse(int e,const IMPULSE& impulse,T dt,T time) const
{
    const COLLISION& col=collisions(e);
    const COMBINED_COLLISIONS_RIGID_IMPULSE<TV>& rigid_impulse=dynamic_cast<const COMBINED_COLLISIONS_RIGID_IMPULSE<TV>&>(impulse);
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.particle_body);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.levelset_body);
    TV rel_vel=rigid_impulse.Velocity_At_Point(col.particle_body,col.location)-rigid_impulse.Velocity_At_Point(col.levelset_body,col.location);
    T orig_rel_vel_N=TV::Dot_Product(col.relative_velocity,col.normal);
    T rel_vel_N=TV::Dot_Product(rel_vel,col.normal);
    T cor=0;
    if(use_coefficient_of_restitution) cor=RIGID_BODY<TV>::Coefficient_Of_Restitution(body1,body2);
    return (1-rel_vel_N/orig_rel_vel_N)/(1+cor);
}
//#####################################################################
// Function Moved
//#####################################################################
template<class TV> void COMBINED_COLLISIONS_RIGID<TV>::
Moved(int e)
{
}
//#####################################################################
// Function Compute_Impulse_Derivative
//#####################################################################
template<class TV> TWIST<TV> COMBINED_COLLISIONS_RIGID<TV>::
Compute_Impulse_Derivative(int e,T scale,const T coefficient_of_restitution,const T coefficient_of_friction,const bool clamp_friction_magnitude) const
{
    const COLLISION& col=collisions(e);
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.particle_body);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.levelset_body);

    T adjusted_cor=(coefficient_of_restitution+1)*scale-1;
    T d_adjusted_cor=coefficient_of_restitution+1;

    T ds=(T)1e-5;

    TWIST<TV> impulse;
    TWIST<TV> d_impulse;
    if(body1.Has_Infinite_Inertia() && body2.Has_Infinite_Inertia()) return TWIST<TV>();
    T relative_normal_velocity=min((T)0,TV::Dot_Product(col.relative_velocity,col.normal));
    typename MATRIX_POLICY<TV>::SYMMETRIC_MATRIX impulse_factor=RIGID_BODY<TV>::Impulse_Factor(body1,body2,col.location),impulse_factor_inverse=impulse_factor.Inverse();
    // see if friction stops sliding
    TV sticking_impulse=impulse_factor_inverse*(-adjusted_cor*relative_normal_velocity*col.normal-col.relative_velocity);
    TV d_sticking_impulse=impulse_factor_inverse*(-d_adjusted_cor*relative_normal_velocity*col.normal);
    T normal_component=TV::Dot_Product(sticking_impulse,col.normal);
    T d_normal_component=TV::Dot_Product(d_sticking_impulse,col.normal);
    if((sticking_impulse-normal_component*col.normal).Magnitude()<=coefficient_of_friction*normal_component){
        impulse.linear=sticking_impulse;
        d_impulse.linear=d_sticking_impulse;
        if(rolling_friction){
            impulse+=RIGID_BODY<TV>::Apply_Rolling_Friction(body1,body2,col.location,col.normal,normal_component);
            TWIST<TV> f=RIGID_BODY<TV>::Apply_Rolling_Friction(body1,body2,col.location,col.normal,normal_component+d_normal_component*ds);
            TWIST<TV> b=RIGID_BODY<TV>::Apply_Rolling_Friction(body1,body2,col.location,col.normal,normal_component-d_normal_component*ds);
            d_impulse+=(f-b)/(2*ds);}}
    // friction does not stop sliding
    else{
        TV relative_tangential_velocity=col.relative_velocity.Projected_Orthogonal_To_Unit_Direction(col.normal);
        TV tangential_direction=relative_tangential_velocity;
        T relative_tangential_velocity_magnitude=tangential_direction.Normalize();
        TV impulse_factor_times_direction=impulse_factor*(col.normal-coefficient_of_friction*tangential_direction);
        PHYSBAM_ASSERT(TV::Dot_Product(impulse_factor_times_direction,col.normal));
        TV delta=-(1+adjusted_cor)*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,col.normal)*impulse_factor_times_direction;
        TV d_delta=-d_adjusted_cor*relative_normal_velocity/TV::Dot_Product(impulse_factor_times_direction,col.normal)*impulse_factor_times_direction;
        if(clamp_friction_magnitude && adjusted_cor){ // should only clamp friction magnitude in the elastic case!
            TV new_relative_velocity=col.relative_velocity+delta;
            TV d_new_relative_velocity=d_delta;
            TV new_normal_velocity=new_relative_velocity.Projected_On_Unit_Direction(col.normal);
            TV d_new_normal_velocity=d_new_relative_velocity.Projected_On_Unit_Direction(col.normal);
            TV new_tangential_velocity=new_relative_velocity-new_normal_velocity;
            TV d_new_tangential_velocity=d_new_relative_velocity-d_new_normal_velocity;
            T new_tangential_velocity_magnitude=new_tangential_velocity.Magnitude();
            T d_new_tangential_velocity_magnitude=TV::Dot_Product(new_tangential_velocity,d_new_tangential_velocity)/new_tangential_velocity_magnitude;
            if(new_tangential_velocity_magnitude > relative_tangential_velocity_magnitude){
                delta=new_normal_velocity+(relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*new_tangential_velocity-col.relative_velocity;
                d_delta=d_new_normal_velocity-relative_tangential_velocity_magnitude*d_new_tangential_velocity_magnitude/sqr(new_tangential_velocity_magnitude)*new_tangential_velocity+
                    (relative_tangential_velocity_magnitude/new_tangential_velocity_magnitude)*d_new_tangential_velocity;}}
        impulse.linear=impulse_factor_inverse*delta;
        d_impulse.linear=impulse_factor_inverse*d_delta;}
    PHYSBAM_ASSERT(!clamp_energy);
/*    if(clamp_energy){
        TWIST<TV> impulse_f=impulse+ds*d_impulse;
        TWIST<TV> impulse_b=impulse+ds*d_impulse;
        Compute_Clamped_Impulse(body1,body2,col.location,impulse,saved_rotation_1,saved_rotation_2);
        Compute_Clamped_Impulse(body1,body2,col.location,impulse_f,saved_rotation_1,saved_rotation_2);
        Compute_Clamped_Impulse(body1,body2,col.location,impulse_b,saved_rotation_1,saved_rotation_2);
        d_impulse=(impulse_f-impulse_b)/(2*ds);}*/

//    TWIST<TV> j=rigid_body_collisions.collision_callbacks.Compute_Collision_Impulse(body1,body2,col.location,col.normal,col.relative_velocity,
//        adjusted_cor,RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),clamp_friction_magnitude,rolling_friction,clamp_energy);
//    TWIST<TV> jf=rigid_body_collisions.collision_callbacks.Compute_Collision_Impulse(body1,body2,col.location,col.normal,col.relative_velocity,
//        adjusted_cor+ds*d_adjusted_cor,RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),clamp_friction_magnitude,rolling_friction,clamp_energy);
//    TWIST<TV> jb=rigid_body_collisions.collision_callbacks.Compute_Collision_Impulse(body1,body2,col.location,col.normal,col.relative_velocity,
//        adjusted_cor-ds*d_adjusted_cor,RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),clamp_friction_magnitude,rolling_friction,clamp_energy);

//    LOG::cout<<"r "<<impulse<<"   "<<d_impulse<<"   "<<j<<"   "<<((jf-jb)/(2*ds))<<std::endl;

    return d_impulse;
}
//#####################################################################
// Function Diagnose_Impulse
//#####################################################################
template<class TV> typename TV::SCALAR COMBINED_COLLISIONS_RIGID<TV>::
Kinetic_Energy_Gradient(int e,const IMPULSE& impulse,T scale,T dt,T time) const
{
    if(scale<0) return 0;
    const COLLISION& col=collisions(e);
    RIGID_BODY<TV>& body1=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.particle_body);
    RIGID_BODY<TV>& body2=rigid_body_collisions.rigid_body_collection.Rigid_Body(col.levelset_body);
    T cor=0;
    if(use_coefficient_of_restitution)
        cor=RIGID_BODY<TV>::Coefficient_Of_Restitution(body1,body2);

    TWIST<TV> j=Compute_Impulse_Derivative(e,scale,cor,RIGID_BODY<TV>::Coefficient_Of_Friction(body1,body2),clamp_friction_magnitude);
    const COMBINED_COLLISIONS_RIGID_IMPULSE<TV>& ri=dynamic_cast<const COMBINED_COLLISIONS_RIGID_IMPULSE<TV>&>(impulse);
    T kp=ri.Kinetic_Energy_Change(col.particle_body,col.location,j);
    T kl=ri.Kinetic_Energy_Change(col.levelset_body,col.location,-j);
    return kp+kl;
}
template class COMBINED_COLLISIONS_RIGID<VECTOR<float,1> >;
template class COMBINED_COLLISIONS_RIGID<VECTOR<float,2> >;
template class COMBINED_COLLISIONS_RIGID<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class COMBINED_COLLISIONS_RIGID<VECTOR<double,1> >;
template class COMBINED_COLLISIONS_RIGID<VECTOR<double,2> >;
template class COMBINED_COLLISIONS_RIGID<VECTOR<double,3> >;
#endif
