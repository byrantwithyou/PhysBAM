//#####################################################################
// Copyright 2004, Eran_Guendelman.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include "KINEMATIC_EXAMPLE.h"
using namespace PhysBAM;
template KINEMATIC_EXAMPLE<double>;
template KINEMATIC_EXAMPLE<float>;
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Initialize_Rigid_Bodies()
{
    std::cout << "*** Using parameter " << parameter << std::endl;

    VECTOR_3D<double> ground_color=RGB_COLORS<T>().Ground_Tan(),color1=RGB_COLORS<T>().Light_Steel_Blue();

    RIGID_BODY<TV> *rigid_body = 0;

    last_frame=600;

    if(parameter==3){
        T baseboxsize=2;
        T dropboxsize=1.2;
        T dropheight=70;
        T smallboxsize=0.5;
        T smallboxmass=1;
        T offset = 0.05;

        T stack_epsilon = 0.3;
        T stack_mu = 0.5;

#if 1
        time1=(T)10/24;
        time2=(T)200/24;
        time3=(T)350/24;
        time4=(T)400/24;

        const char *boxfile = "square";

        T boxsize1=1.0;
        T boxsize2=0.5;

        rigid_body=Initialize_Rigid_Body(boxfile,boxsize1);
        rigid_body->position=VECTOR_2D<T>(0, 2*baseboxsize+boxsize1);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1a");
        Append_Diffuse_Hints(color1);

#if 1
        rigid_body=Initialize_Rigid_Body(boxfile,boxsize2);
        rigid_body->position=VECTOR_2D<T>(0, 2*baseboxsize+2*boxsize1+boxsize2);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1b");
        Append_Diffuse_Hints(color1);
#endif

        rigid_body=Initialize_Rigid_Body(boxfile,baseboxsize);
        rigid_body->Set_Coefficient_Of_Restitution(0.1);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Name("base box");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index) = true;
        kinematic_body_index=rigid_body_parameters.list.rigid_bodies.m;
        kinematic_body_vel=5;
        Append_Diffuse_Hints(color1);

#else
        const char *boxfile = "square_refined";

        time1=(T)100/24;
        time2=(T)200/24;
        time3=(T)350/24;
        time4=(T)400/24;

        rigid_body=Initialize_Rigid_Body(boxfile,baseboxsize);
        rigid_body->Set_Coefficient_Of_Restitution(0.1);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Name("base box");
        rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index) = true;
        kinematic_body_index=rigid_body_parameters.list.rigid_bodies.m;
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_2D<T>(-0.65, baseboxsize+4);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1a");
        Append_Diffuse_Hints(color1);

#if 0
        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0.65, baseboxsize+5, 0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 1b");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(offset, baseboxsize+7, 0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 2");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0, baseboxsize+9, offset);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 3");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(-offset, baseboxsize+14, -offset);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 4");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0, baseboxsize+17, 0);
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 5");
        Append_Diffuse_Hints(color1);

        rigid_body=Initialize_Rigid_Body(boxfile,smallboxsize);
        rigid_body->position=VECTOR_3D<T>(0, baseboxsize+19, 0);
        rigid_body->orientation=QUATERNION<T>(pi/4,VECTOR_3D<T>(0,1,0));
        rigid_body->Set_Coefficient_Of_Restitution(stack_epsilon);
        rigid_body->Set_Coefficient_Of_Friction(stack_mu);
        rigid_body->Set_Mass(smallboxmass);
        rigid_body->Set_Name("stack box 6");
        Append_Diffuse_Hints(color1);
#endif
#endif

        //for(int i=1;i<=rigid_body_parameters.list.rigid_bodies.m;i++) if(!rigid_body_parameters.list.rigid_bodies(i)->is_static) rigid_body_parameters.list.rigid_bodies(i)->is_kinematic=true;
    }

#if 0
    // PLANE -- no gravity
    rigid_body=Initialize_Rigid_Body("ground");
    rigid_body->Set_Coefficient_Of_Restitution(1.0);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;
    Append_Diffuse_Hints(ground_color,false);
#endif

    Update_Animated_Parameters(0);
    for(int i=1;i<=rigid_body_parameters.list.rigid_bodies.m;i++){
        rigid_body=rigid_body_parameters.list.rigid_bodies(i);
        if(rigid_body->rigid_body_collection.rigid_body_particle.kinematic(rigid_body->particle_index)){
            Set_External_Position_And_Orientation(rigid_body->position,rigid_body->orientation,0,i);
            Set_External_Velocities(rigid_body->velocity,rigid_body->angular_velocity,0,i);}}

    rigid_body_parameters.list.Set_External_Forces_And_Velocities(*this);
}
//#####################################################################
// Function Compute_State
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Compute_State(RIGID_BODY_STATE_2D<T>& state,const T time)
{
    T baseboxsize=2;
    state.time=time;
    if(time<time1){
        state.position=VECTOR_2D<T>(0,baseboxsize);
        state.orientation=0;
    }
    else if(time<time2){
        state.position=VECTOR_2D<T>(0,baseboxsize+kinematic_body_vel*(time-time1));
        state.orientation=0;
    }
    else if(time<time3){
        state.position=VECTOR_2D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1));
        state.orientation=0;
    }
    else if(time<time4){
        state.position=VECTOR_2D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1));
        state.orientation=0.4*(time-time3);
    }
    else{
        state.position=VECTOR_2D<T>(0,baseboxsize+kinematic_body_vel*(time2-time1));
        state.orientation=-0.4*(time4-time3);
    }
}
//#####################################################################
// Function Compute_Velocities
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Compute_Velocities(const RIGID_BODY<TV>& rigid_body,RIGID_BODY_STATE_2D<T>& state1,RIGID_BODY_STATE_2D<T>& state2)
{
    T one_over_dt=1/(state2.time-state1.time);
    state2.velocity=one_over_dt*(state2.position-state1.position);
    state2.angular_velocity=one_over_dt*(state2.orientation-state1.orientation);
    rigid_body.Update_Angular_Momentum(state2);
}
//#####################################################################
// Function Update_Animated_Parameters
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Update_Animated_Parameters(const T time)
{
    if(parameter==3){
        T prev_time=time-(T)1/frame_rate,current_time=time,next_time=time+(T)1/frame_rate;
        RIGID_BODY_STATE_2D<T> prev_state;
        Compute_State(prev_state,prev_time);
        Compute_State(current_state,current_time);
        Compute_State(next_state,next_time);
        rigid_body_list(kinematic_body_index).Compute_Velocity_Between_States(prev_state,current_state,current_state);
        rigid_body_list(kinematic_body_index).Compute_Velocity_Between_States(current_state,next_state,next_state);
    }
}    
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Set_External_Velocities(VECTOR_2D<T>& V,T& omega,const T time,const int id_number)
{
    if(!rigid_body_parameters.list.rigid_bodies(id_number)->is_kinematic) return;
    else if(parameter==3){
        if(id_number==kinematic_body_index){
            RIGID_BODY_STATE_2D<T> interpolated_state;
            rigid_body_parameters.list.rigid_bodies(id_number)->Interpolate_Between_States(current_state,next_state,time,interpolated_state);
            V=interpolated_state.velocity;
            omega=interpolated_state.angular_velocity;
        }
    }
}
//#####################################################################
// Function Set_External_Position_And_Orientation
//#####################################################################
template<class T> void KINEMATIC_EXAMPLE<T>::
Set_External_Position_And_Orientation(VECTOR_2D<T>& X,T& orientation,const T time,const int id_number)
{
    if(!rigid_body_parameters.list.rigid_bodies(id_number)->is_kinematic) return;
    if(parameter==3){
        if(id_number==kinematic_body_index){
            RIGID_BODY_STATE_2D<T> interpolated_state;
            rigid_body_parameters.list.rigid_bodies(id_number)->Interpolate_Between_States(current_state,next_state,time,interpolated_state);
            X=interpolated_state.position;
            orientation=interpolated_state.orientation;
        }
    }
}
