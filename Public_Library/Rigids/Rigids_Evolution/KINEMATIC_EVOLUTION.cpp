//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Forces_And_Torques/RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigids_Evolution/KINEMATIC_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KINEMATIC_EVOLUTION<TV>::
KINEMATIC_EVOLUTION(RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,RIGIDS_EXAMPLE_FORCES_AND_VELOCITIES<TV>& rigid_body_example_velocities)
    :rigid_body_collection(rigid_body_collection_input),rigid_body_example_velocities(rigid_body_example_velocities)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KINEMATIC_EVOLUTION<TV>::
~KINEMATIC_EVOLUTION()
{
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=rigid_body_collection.kinematic_rigid_bodies(i);Set_External_Velocities(twist(p),velocity_time,p);}
    rigid_body_example_velocities.Set_External_Velocities(twist,velocity_time,current_position_time);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(id);
    if(rigid_body.is_static){twist=TWIST<TV>();return;}
    Set_Kinematic_Velocities(twist,(T)1e-3,time,id);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T frame_dt,const T time,const int id)
{
    RIGID_BODY<TV>* rigid_body=&rigid_body_collection.Rigid_Body(id);
    if(rigid_body_example_velocities.Set_Kinematic_Velocities(twist,time,id)) return;
    RIGID_BODY_STATE<TV> previous_state,current_state;previous_state.time=time-frame_dt;current_state.time=time;
    int new_id=id;
    rigid_body_example_velocities.Set_Kinematic_Positions(previous_state.frame,previous_state.time,new_id);
    rigid_body_example_velocities.Set_Kinematic_Positions(current_state.frame,current_state.time,new_id);
    rigid_body->Compute_Velocity_Between_States(previous_state,current_state,current_state);
    twist=current_state.twist;
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(FRAME<TV>& frame,const T time,const int id)
{
    rigid_body_example_velocities.Set_Kinematic_Positions(frame,time,id);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time)
{
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){
        int p=rigid_body_collection.kinematic_rigid_bodies(i);Set_External_Positions(frame(p),time,p);}
    rigid_body_example_velocities.Set_External_Positions(frame,time);
}
//#####################################################################
// Function Reset_Kinematic_Rigid_Bodies
//#####################################################################
template<class TV> void KINEMATIC_EVOLUTION<TV>::
Reset_Kinematic_Rigid_Bodies(const T time)
{
    // Move kinematic bodies to their position at given time
    for(int i=0;i<rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(rigid_body_collection.kinematic_rigid_bodies(i));
        Set_External_Positions(rigid_body.Frame(),time,rigid_body.particle_index);
        Set_External_Velocities(rigid_body.Twist(),time,rigid_body.particle_index);}
}
//#####################################################################
namespace PhysBAM{
template class KINEMATIC_EVOLUTION<VECTOR<float,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<float,3> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,1> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,2> >;
template class KINEMATIC_EVOLUTION<VECTOR<double,3> >;
}
