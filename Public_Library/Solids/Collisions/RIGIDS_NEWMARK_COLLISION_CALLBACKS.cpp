//#####################################################################
// Copyright 2009.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/FRAME.h>
#include <Core/Vectors/TWIST.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <Solids/Collisions/RIGIDS_NEWMARK_COLLISION_CALLBACKS.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
RIGIDS_NEWMARK_COLLISION_CALLBACKS(NEWMARK_EVOLUTION<TV>& evolution_input)
    :evolution(evolution_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
~RIGIDS_NEWMARK_COLLISION_CALLBACKS()
{
}
//#####################################################################
// Function Get_Body_Penetration
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Reevolve_Body_With_Saved_State(const int p,const T dt,const T time)
{
    Restore_Position(p);
    // re-evolve this body using new velocity
    Euler_Step_Position_With_New_Velocity(p,dt,time);
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Restore_Positions()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    for(int i=0;i<rigid_body_collection.dynamic_rigid_body_particles.m;i++)
        Restore_Position(rigid_body_collection.dynamic_rigid_body_particles(i));
    for(int i=0;i<rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++)
        Restore_Position(rigid_body_collection.static_and_kinematic_rigid_bodies(i));
}
//#####################################################################
// Function Restore_Position
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Restore_Position(const int p)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    rigid_body_collection.rigid_body_particles.frame(p)=evolution.rigid_frame_save(p);
    rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();
    rigid_body_collection.Rigid_Body(p).Update_Bounding_Box();
}
//#####################################################################
// Function Save_Position
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Save_Position(const int p)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    evolution.rigid_frame_save(p)=rigid_body_collection.rigid_body_particles.frame(p);
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Restore_Velocity(const int p)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    rigid_body_collection.rigid_body_particles.twist(p).linear=evolution.rigid_velocity_save(p).linear;
    rigid_body_collection.Rigid_Body(p).Angular_Momentum()=evolution.rigid_angular_momentum_save(p);
    rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Save_Velocity(const int p)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    evolution.rigid_velocity_save(p)=rigid_body_collection.rigid_body_particles.twist(p);
    evolution.rigid_angular_momentum_save(p)=rigid_body_collection.Rigid_Body(p).Angular_Momentum();
}
//#####################################################################
// Function Euler_Step_Position
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Euler_Step_Position(const int id,const T dt,const T time)
{
    evolution.Euler_Step_Position(dt,time,id);
    evolution.solid_body_collection.rigid_body_collection.Rigid_Body(id).Update_Bounding_Box();
}
//#####################################################################
// Function Euler_Step_Position_With_New_Velocity
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Euler_Step_Position_With_New_Velocity(const int id,const T dt,const T time)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    Restore_Position(id);
    TV velocity=rigid_body_collection.rigid_body_particles.twist(id).linear;T_SPIN angular_momentum=rigid_body_collection.Rigid_Body(id).Angular_Momentum(); // save velocity
    evolution.Update_Velocity_Using_Stored_Differences(dt,time,id); // temporarily update velocity
    Euler_Step_Position(id,dt,time);
    rigid_body_collection.rigid_body_particles.twist(id).linear=velocity;
    rigid_body_collection.Rigid_Body(id).Angular_Momentum()=angular_momentum; // restore velocity
    rigid_body_collection.Rigid_Body(id).Update_Angular_Velocity(); // re-sync this
}
//#####################################################################
// Function Swap_States
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Swap_State(const int id)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_body_cluster_bindings=rigid_body_collection.rigid_body_cluster_bindings;
    RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(id);
    int parent_id=rigid_body_cluster_bindings.Get_Parent(body).particle_index;
    RIGID_BODY<TV>& parent_body=rigid_body_collection.Rigid_Body(parent_id);
    exchange(rigid_body_collection.rigid_body_particles.frame(id),evolution.rigid_frame_save(body.particle_index));
    if(parent_id!=id) exchange(rigid_body_collection.rigid_body_particles.frame(parent_id),evolution.rigid_frame_save(parent_id));
    parent_body.Update_Angular_Velocity();
}
//#####################################################################
// Function Saved_Particle_To_Levelset_Body_Transform
//#####################################################################
template<class TV> FRAME<TV> RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Saved_Particle_To_Levelset_Body_Transform(const int levelset_body,const int particle_body)
{
    return evolution.rigid_frame_save(levelset_body).Inverse_Times(evolution.rigid_frame_save(particle_body));
}
//#####################################################################
// Function Exchange_Frame
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Exchange_Frame(const int id)
{
    exchange(evolution.solid_body_collection.rigid_body_collection.rigid_body_particles.frame(id),evolution.rigid_frame_save(id));
}
//#####################################################################
// Function Compute_Collision_Impulse
//#####################################################################
template<class TV> TWIST<TV> RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Compute_Collision_Impulse(RIGID_BODY<TV>& body0,RIGID_BODY<TV>& body1,const TV& location,const TV& normal,const TV& relative_velocity,const T coefficient_of_restitution,
    const T coefficient_of_friction,const bool clamp_friction_magnitude,const bool rolling_friction,const bool clamp_energy)
{
    return RIGID_BODY<TV>::Compute_Collision_Impulse(body0,body1,evolution.rigid_frame_save(body0.particle_index).r,evolution.rigid_frame_save(body1.particle_index).r,location,
        normal,relative_velocity,coefficient_of_restitution,coefficient_of_friction,clamp_friction_magnitude,rolling_friction,clamp_energy);
}
//#####################################################################
// Function Subtract_Stored_Differences
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Subtract_Stored_Difference(TV& velocity,T_SPIN& momentum,const int particle_index)
{
    velocity-=evolution.rigid_velocity_difference(particle_index);
    momentum-=evolution.rigid_angular_momentum_difference(particle_index);
}
//#####################################################################
// Function Begin_Fracture
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Begin_Fracture(const int body_id)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    Restore_Position(body_id);Restore_Velocity(body_id);
    old_stored_difference.angular=evolution.rigid_angular_momentum_difference(body_id);
    old_stored_difference.linear=evolution.rigid_velocity_difference(body_id);
    old_position=rigid_body_collection.Rigid_Body(body_id).Frame().t;
}
//#####################################################################
// Function End_Fracture
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
End_Fracture(const int body_id,ARRAY<int>& added_bodies)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=evolution.solid_body_collection.rigid_body_collection;
    int new_size=rigid_body_collection.rigid_body_particles.Size();
    evolution.rigid_frame_save.Resize(new_size);
    evolution.rigid_velocity_save.Resize(new_size);
    evolution.rigid_angular_momentum_save.Resize(new_size);
    evolution.rigid_velocity_difference.Resize(new_size);
    evolution.rigid_angular_momentum_difference.Resize(new_size);

    for(int j=0;j<added_bodies.m;j++){
        Save_Position(added_bodies(j));
        Save_Velocity(added_bodies(j));
        evolution.rigid_angular_momentum_difference(body_id)=old_stored_difference.angular;
        evolution.rigid_velocity_difference(body_id)=
            old_stored_difference.linear+TV::Cross_Product(old_stored_difference.angular,rigid_body_collection.rigid_body_particles.frame(body_id).t-old_position);}
}
//#####################################################################
// Function Begin_Asymmetric_Collisions
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
Begin_Asymmetric_Collisions(const int body_1,const int body_2)
{
    asymmetric_collision_stored_difference.angular=evolution.rigid_angular_momentum_difference(body_2);
    asymmetric_collision_stored_difference.linear=evolution.rigid_velocity_difference(body_2);
    evolution.rigid_angular_momentum_difference(body_2)=T_SPIN();
    evolution.rigid_velocity_difference(body_2)=TV();
}
//#####################################################################
// Function End_Asymmetric_Collisions
//#####################################################################
template<class TV> void RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>::
End_Asymmetric_Collisions(const int body_1,const int body_2,VECTOR<ARRAY<int>,2>& added_bodies)
{
    evolution.rigid_angular_momentum_difference(body_2)=asymmetric_collision_stored_difference.angular;
    evolution.rigid_velocity_difference(body_2)=asymmetric_collision_stored_difference.linear;
    if(added_bodies(1).m){
        for(int body=0;body<added_bodies(1).m;body++){
            int new_body_id=added_bodies(1)(body);
            evolution.rigid_angular_momentum_difference(new_body_id)=asymmetric_collision_stored_difference.angular;
            evolution.rigid_velocity_difference(new_body_id)=asymmetric_collision_stored_difference.linear;}}
    else{
        evolution.rigid_angular_momentum_difference(body_2)=asymmetric_collision_stored_difference.angular;
        evolution.rigid_velocity_difference(body_2)=asymmetric_collision_stored_difference.linear;}

    for(int body=0;body<added_bodies(0).m;body++) Restore_Position(added_bodies(0)(body));
}
//#####################################################################
namespace PhysBAM{
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<float,1> >;
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<float,2> >;
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<float,3> >;
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<double,1> >;
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<double,2> >;
template class RIGIDS_NEWMARK_COLLISION_CALLBACKS<VECTOR<double,3> >;
}
