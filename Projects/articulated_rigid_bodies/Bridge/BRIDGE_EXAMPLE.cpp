//#####################################################################
// Copyright 2004-2007, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <map>
#include "BRIDGE_EXAMPLE.h"
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T> BRIDGE_EXAMPLE<T>::
BRIDGE_EXAMPLE(const STREAM_TYPE stream_type)
    :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,output_directory,data_directory,solid_body_collection),box1(0),box2(0),start_rolling(100),num_rolling_frames(1000),start_drop(2000),num_rungs(10),
    selection(0)
{
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.rigid_body_collision_parameters.rigid_collisions_particle_partition_size=7;

    restart=false;
    frame_rate=24*4;
    last_frame=4000;
    write_last_frame=true;
}
//#####################################################################
// Destructor
//#####################################################################
template<class T> BRIDGE_EXAMPLE<T>::
~BRIDGE_EXAMPLE()
{
}
template<class T> void BRIDGE_EXAMPLE<T>::
Register_Options()
{
    BASE::Register_Options();
    parse_args->Add("-selection",&selection,"value","selection");
}
template<class T> void BRIDGE_EXAMPLE<T>::
Parse_Options()
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Bridge/%s%s",output_directory.c_str(),(selection==0?"_blocks":"_lathe_chains"));
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Initialize_Bodies()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
//    arb.Set_Write_Extra_Restart_Data(true);
    arb.Set_Use_Shock_Propagation(true);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);

    Make_Bridge();
    tests.Add_Ground((T).5,-2);
    if(selection==0) Make_Blocks();else Make_Lathe_Chains();

    tests.Add_Gravity();

    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Update_Solids_Parameters(const T time)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    int frame=(int)(first_frame+(time-initial_time)*frame_rate);
    if(frame<start_drop){
        T max_linear_velocity=(T).1,max_angular_velocity=(T).1;
        solids_parameters.rigid_body_evolution_parameters.clamp_rigid_body_velocities=true;
        solids_parameters.rigid_body_evolution_parameters.max_rigid_body_linear_velocity=max_linear_velocity;
        solids_parameters.rigid_body_evolution_parameters.max_rigid_body_angular_velocity=max_angular_velocity;
        arb.Set_Iterative_Tolerance((T)1e-3);
        arb.Set_Max_Iterations(100);
        for(int i=box2->particle_index+2;i<=rigid_body_collection.rigid_body_particles.Size();i++) rigid_body_collection.Rigid_Body(i).is_static=true;}
    else{
        solids_parameters.rigid_body_evolution_parameters.clamp_rigid_body_velocities=false;
        arb.Set_Iterative_Tolerance((T)1e-5);
        arb.Set_Max_Iterations(1000);
        for(int i=box2->particle_index+2;i<=rigid_body_collection.rigid_body_particles.Size();i++) rigid_body_collection.Rigid_Body(i).is_static=false;}
}
//#####################################################################
// Function Make_Lathe_Chains
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Make_Lathe_Chains()
{
    for(int k=0;k<=2;k++) for(int i=0;i<20;i++)
        tests.Make_Lathe_Chain(FRAME<TV>(TV((T)(-4+(T).8755*(i+(T).5*(k%2))*(T).5),(T)(5-k),(i+k)%2?-(T).25:(T).25),ROTATION<TV>((T)pi/2,TV(0,1,0))),(T).125);
}
//#####################################################################
// Function Make_Blocks
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Make_Blocks()
{
    for(int i=-4;i<=4;i++) for(int j=-1;j<=1;j++) for(int k=-1;k<=1;k++) Point_Constraint_With_2_Blocks(TV(i*(T).75,4+j*(T).6,k*(T).5));
}
//#####################################################################
// Function Point_Constraint_With_2_Blocks
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Point_Constraint_With_2_Blocks(TV shift)
{   
    int parent_id=Add_Rigid_Body("subdivided_box",(T).1,(T).5,(T).5,FRAME<TV>(shift),"parent").particle_index;
    int child_id=Add_Rigid_Body("subdivided_box",(T).1,(T).5,(T).5,FRAME<TV>(TV((T).2,(T).2,(T).2)+shift),"child").particle_index;
    Add_Joint(parent_id,child_id,new POINT_JOINT<TV>,FRAME<TV>(TV((T)-.1,(T)-.1,(T)-.1)));
}
//#####################################################################
// Function Make_Bridge
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Make_Bridge()
{
    // planks and posts
    T x_shift=-4,rotation=(T)pi/2;
    for(int row=0;row<num_rungs-1;row++){
        Add_Rigid_Body("plank",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*row,(T)-.53125,0)),"bridge");
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*row,0,(T).5),ROTATION<TV>(rotation,TV(0,1,0))),STRING_UTILITIES::string_sprintf("rail%d",(row*5+2)));
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*row,0,(T)-.5),ROTATION<TV>(rotation,TV(0,1,0))),STRING_UTILITIES::string_sprintf("rail%d",(row*5+3)));
        Add_Rigid_Body("plank",(T).06,(T).5,0,FRAME<TV>(TV(x_shift+(T).43775+(T).8755*row,(T)(-.53125+.015625),.5625),ROTATION<TV>((T)pi/2,TV(0,1,0))),"side_bridge",10);
        Add_Rigid_Body("plank",(T).06,(T).5,0,FRAME<TV>(TV(x_shift+(T).43775+(T).8755*row,(T)(-.53125+.015625),(T)-.5625),ROTATION<TV>((T)pi/2,TV(0,1,0))),"side_bridge",10);}

    Add_Rigid_Body("plank",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*(num_rungs-1),(T)-.53125,0)),"bridge");
    Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*(num_rungs-1),0,(T).5),ROTATION<TV>(rotation,TV(0,1,0))),STRING_UTILITIES::string_sprintf("rail%d",(9*5+2)));
    Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).125,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*(num_rungs-1),0,(T)-.5),ROTATION<TV>(rotation,TV(0,1,0))),STRING_UTILITIES::string_sprintf("rail%d",(9*5+3)));

    // ropes
    T rope_length=(T).75/(1+sqrt((T)3)),y1=(T).5-(T).5*rope_length*sin((T)pi/6)-(T).05,y2=(T).5-rope_length*sin((T)pi/6)-(T).05;
    for(int row=0;row<num_rungs-1;row++){
        T x1=x_shift+(T).0625+(T).8755*row+(T).5*rope_length*cos((T)pi/6),x2=x1+(T).5*rope_length*(1+cos((T)pi/6)),x3=x2+(T).5*rope_length*(1+cos((T)pi/6));
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x1,y1,(T).5),ROTATION<TV>((T)pi/3,TV(0,0,1))),"rope",50);
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x2,y2,(T).5),ROTATION<TV>((T)pi/2,TV(0,0,1))),"rope",50);
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x3,y1,(T).5),ROTATION<TV>(-(T)pi/3,TV(0,0,1))),"rope",50);
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x1,y1,(T)-.5),ROTATION<TV>((T)pi/3,TV(0,0,1))),"rope",50);
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x2,y2,(T)-.5),ROTATION<TV>((T)pi/2,TV(0,0,1))),"rope",50);
        Add_Rigid_Body("Rings_Test/cylinder_revolve",(T).03,(T).5,0,FRAME<TV>(TV(x3,y1,(T)-.5),ROTATION<TV>(-(T)pi/3,TV(0,0,1))),"rope",50);}

    // bases
    box1=&Add_Rigid_Body("plank",(T)2.25,(T).5,0,FRAME<TV>(TV(x_shift-(T)2.3755,(T)-11.8125,0),ROTATION<TV>((T)pi/2,TV(1,0,0))),"base1");box1->is_static=true;
    box2=&Add_Rigid_Body("plank",(T)2.25,(T).5,0,FRAME<TV>(TV(x_shift+(T).8755*(num_rungs-1)+(T)2.3755,(T)-11.8125,0),ROTATION<TV>((T)pi/2,TV(1,0,0))),"base2");box2->is_static=true;

    // joints
    LOG::cout<<"adding articulation"<<std::endl;
    for(int i=0;i<num_rungs;i++){int rung=i*5+1;
        // rung to post
        for(int side=0;side<2;side++) Add_Joint(int(rung),int(rung+side),new RIGID_JOINT<TV>,FRAME<TV>(TV(0,(T)-.5,0)));
        // rung to sides
        if(rung>1) for(int side=-2;side<=-1;side++){
            ANGLE_JOINT<TV>& joint=static_cast<ANGLE_JOINT<TV>&>(Add_Joint(int(rung),int(rung+side),new ANGLE_JOINT<TV>,FRAME<TV>(TV(0,(T).015625,(T).3))));
            joint.Set_Angle_Constraints(true,-(T)pi/12,(T)pi/12);}
        if(rung<5*(num_rungs-1)+1) for(int side=3;side<=4;side++){
            ANGLE_JOINT<TV>& joint=static_cast<ANGLE_JOINT<TV>&>(Add_Joint(int(rung),int(rung+side),new ANGLE_JOINT<TV>,FRAME<TV>(TV(0,(T).015625,(T)-.3))));
            joint.Set_Angle_Constraints(true,-(T)pi/12,(T)pi/12);}}

    int num_bodies=5*(num_rungs-1)+3;
    for(int i=0;i<num_rungs-1;i++){
        Add_Joint(int(i*5+2),int(num_bodies+i*6+1),new POINT_JOINT<TV>,FRAME<TV>(TV(0,rope_length*(T).5,0))); // ropes
        Add_Joint(int(num_bodies+i*6+1),int(num_bodies+i*6+2),new POINT_JOINT<TV>,FRAME<TV>(TV(0,rope_length*(T).5,0)));
        Add_Joint(int(num_bodies+i*6+2),int(num_bodies+i*6+3),new POINT_JOINT<TV>,FRAME<TV>(TV(0,-rope_length*(T).5,0)));
        Add_Joint(int(num_bodies+i*6+3),int(i*5+7),new POINT_JOINT<TV>,FRAME<TV>(TV(0,(T).45,(T)-.0625)));
        Add_Joint(int(i*5+3),int(num_bodies+i*6+4),new POINT_JOINT<TV>,FRAME<TV>(TV(0,rope_length*(T).5,0)));
        Add_Joint(int(num_bodies+i*6+4),int(num_bodies+i*6+5),new POINT_JOINT<TV>,FRAME<TV>(TV(0,rope_length*(T).5,0)));
        Add_Joint(int(num_bodies+i*6+5),int(num_bodies+i*6+6),new POINT_JOINT<TV>,FRAME<TV>(TV(0,-rope_length*(T).5,0)));
        Add_Joint(int(num_bodies+i*6+6),int(i*5+8),new POINT_JOINT<TV>,FRAME<TV>(TV(0,(T).45,(T)-.0625)));}

    // attach edge rungs to bases
    JOINT<TV>* joint=&Add_Joint(int(1),box1->particle_index,new ANGLE_JOINT<TV>,FRAME<TV>(TV((T)2.25,0,(T)-11.25),ROTATION<TV>((T)pi/2,TV(0,0,1))));
    static_cast<ANGLE_JOINT<TV>*>(joint)->Set_Angle_Constraints(true,-(T)pi/12,(T)pi/12);
    joint=&Add_Joint(int(5*(num_rungs-1)+1),box2->particle_index,new ANGLE_JOINT<TV>,FRAME<TV>(TV((T)-2.25,0,(T)-11.25),ROTATION<TV>((T)pi/2,TV(0,0,1))));
    static_cast<ANGLE_JOINT<TV>*>(joint)->Set_Angle_Constraints(true,-(T)pi/12,(T)pi/12);
}
//#####################################################################
// Function Add_Rigid_Body
//#####################################################################
template<class T> RIGID_BODY<VECTOR<T,3> >& BRIDGE_EXAMPLE<T>::
Add_Rigid_Body(const std::string& filename,const T scale,const T cof,const T cor,const FRAME<TV>& frame,const std::string& name,const T mass_scale)
{
    RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body(filename,scale,cof);
    rigid_body.Frame()=frame;
    rigid_body.Set_Coefficient_Of_Restitution(cor);
    rigid_body.name=name;
    if(mass_scale) rigid_body.Set_Mass(rigid_body.Mass()*mass_scale);
    return rigid_body;
}
//#####################################################################
// Function Add_Joint
//#####################################################################
template<class T> JOINT<VECTOR<T,3> >& BRIDGE_EXAMPLE<T>::
Add_Joint(const int parent_id,const int child_id,JOINT<TV>* joint,const FRAME<TV>& joint_to_child_frame)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    arb.joint_mesh.Add_Articulation(parent_id,child_id,joint);
    joint->Set_Joint_To_Child_Frame(joint_to_child_frame);
    joint->Set_Joint_To_Parent_Frame(rigid_body_collection.Rigid_Body(parent_id).Frame().Inverse()*rigid_body_collection.Rigid_Body(child_id).Frame()*joint->F_cj());

    return *joint;
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class T> void BRIDGE_EXAMPLE<T>::
Preprocess_Frame(const int frame)
{
    if(frame>start_rolling && frame<=start_rolling+num_rolling_frames){
        box1->Frame().t+=TV((T).1/num_rolling_frames,0,0);
        box2->Frame().t-=TV((T).1/num_rolling_frames,0,0);}
}
//#####################################################################
namespace PhysBAM{
template class BRIDGE_EXAMPLE<float>;
template class BRIDGE_EXAMPLE<double>;
}
