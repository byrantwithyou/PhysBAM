//#####################################################################
// Copyright 2004-2008, Craig Schroeder, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class TANK_EXAMPLE
//##################################################################### 
#ifndef __TANK_EXAMPLE__
#define __TANK_EXAMPLE__

#include <Tools/Read_Write/FILE_UTILITIES.h>
#include <Tools/Read_Write/READ_WRITE_FUNCTIONS.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE.h>
#include <fstream>
#include <iostream>
#include <map>

namespace PhysBAM{

template<class T_input>
class TANK_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef SOLIDS_EXAMPLE<TV> BASE;
public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;using BASE::solids_parameters;using BASE::data_directory;using BASE::stream_type;
    using BASE::test_number;using BASE::solid_body_collection;using BASE::Set_External_Velocities; // silence -Woverloaded-virtual

    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    ARRAY<RIGID_BODY<TV>*> treads[2];
    RIGID_BODY<TV> *gears[6],*body,*gun,*lid;
    ARRAY<JOINT<TV>*> dynamic_joints;
    T gear_radius,speed;
    T start_rolling,accelerate_time;
    bool turn_in_place,start_on_fridge;
    T initial_height;
    int next_slot[68],next_gear[68];
    FRAME<TV> parent_to_joint[10];

    TANK_EXAMPLE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
        :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),turn_in_place(false),start_on_fridge(false)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        last_frame=2000;
        frame_rate=24;
        LOG::cout<<"Frame rate: "<<frame_rate<<std::endl;
        start_rolling=(T).8;
        accelerate_time=(T).8;
        speed=12;
        gear_radius=(T).7;
        parse_args.Add("-tank_on_fridge",&start_on_fridge,"For the tank example, start on fridge blocks");
        parse_args.Add("-tank_turn",&turn_in_place,"For the tank example, turn in place instead of move forward");
        parse_args.Parse();
        tests.data_directory=data_directory;
        output_directory="Tank/output";
    }

    virtual ~TANK_EXAMPLE()
    {
        delete arb;
    }

    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}

void After_Initialization() PHYSBAM_OVERRIDE {BASE::After_Initialization();}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Use_Shock_Propagation(true);
    arb->Set_Contact_Level_Iterations(5);
    arb->Set_Shock_Propagation_Level_Iterations(5);
    arb->Set_Poststabilization_Iterations(5);

    initial_height=start_on_fridge?(T)44:0;
    Make_Tread((T)2.625,0);
    Make_Tread((T)-2.625,1);
    Make_Center();
    Make_Gun();
    Initialize_Tread_Tracking();
    Create_Static_Joints();
    Reset_Joints();

    if(start_on_fridge){
        RIGID_BODY<TV>* rigid_body=0;
        rigid_body=&tests.Add_Rigid_Body("subdivided_box",11,(T)1);
        rigid_body->Frame().t=TV(0,(T)10.21,0);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->name="fridge";
        rigid_body->is_static=true;
        rigid_body=&tests.Add_Rigid_Body("subdivided_box",11,(T)1);
        rigid_body->Frame().t=TV(0,(T)32.21,0);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->name="fridge";
        rigid_body->is_static=true;}

    tests.Add_Ground(1,-(T).79,0);
    
    tests.Add_Gravity();
    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_EXAMPLE<TV>::Initialize_Bodies();
}
//#####################################################################
// Function Make_Gun
//#####################################################################
void Make_Gun()
{
    gun=&tests.Add_Rigid_Body("ARB/gun",1,(T)1);
    gun->Frame().t=TV(0,initial_height+(T)1.3,0);
    gun->Set_Coefficient_Of_Restitution(0);
    gun->name="gun";

    lid=&tests.Add_Rigid_Body("ARB/tanklid",(T)1.05,(T)1);
    lid->Frame().t=TV((T)-.5,initial_height+(T)2.3,0);//+TV(-7,(T)33.79,11);
    lid->Frame().r=ROTATION<TV>((T)pi/4,TV(0,0,1));
    lid->Set_Coefficient_Of_Restitution(0);
    lid->name="lid";
}
//#####################################################################
// Function Make_Center
//#####################################################################
void Make_Center()
{
    // axel - num_joints+2
    body=&tests.Add_Rigid_Body("ARB/tankbody",(T)1.025,(T)1);
    body->Frame().t=TV(0,initial_height+(T).25,0);
    body->Set_Coefficient_Of_Restitution(0);
    body->name="tankbody";
}
//#####################################################################
// Function Make_Tread
//#####################################################################
void Make_Tread(const T z_shift,int tread_side)
{
    RIGID_BODY<TV>* rigid_body=0;
    // tread joints   
    TV move(0,0,0);
    int num_joints=11;
    T split=(T).39;
    T start=-(split*(num_joints+1))*(T).5;
    // tread rigid bodies - num_joints+1
    int tread_num=1;
    T radius=(T).685;
    for(int k=0;k<num_joints+1;k++){
        rigid_body=&tests.Add_Rigid_Body("ARB/tread3_subdivided",(T).5,(T)1);
        rigid_body->Frame().t=TV(start+k*split+(T).5*split,initial_height+radius,z_shift)+move;
        rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi/2,TV(1,0,0));
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Mass(rigid_body->Mass()*10); // NOTE: was rigid_body->mass*=10
        rigid_body->name=LOG::sprintf("tread%d",tread_num);tread_num++;
        treads[tread_side].Append(rigid_body);}
    // end 1
    T radius2=(T).65;
    T angle_start=2*(T)pi/5;T angle=(T)pi/5;
    for(int k=0;k<5;k++){
        rigid_body=&tests.Add_Rigid_Body("ARB/tread3_subdivided",(T).5,(T)1);
        rigid_body->Frame().t=TV(-start+radius2*cos(angle_start-k*angle),initial_height+radius2*sin(angle_start-k*angle),z_shift)+move;
        rigid_body->Frame().r=ROTATION<TV>(-(k+(T).5)*angle,TV(0,0,1))*ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi/2,TV(1,0,0));
        rigid_body->Set_Coefficient_Of_Restitution(0); 
        rigid_body->Set_Mass(rigid_body->Mass()*10); // NOTE: was rigid_body->mass*=10
        rigid_body->name=LOG::sprintf("tread%d",tread_num);tread_num++;
        treads[tread_side].Append(rigid_body);}
    // bottom treads
    for(int k=num_joints;k>=0;k--){
        rigid_body=&tests.Add_Rigid_Body("ARB/tread3_subdivided",(T).5,(T)1);
        rigid_body->Frame().t=TV(start+k*split+(T).5*split,initial_height-radius,z_shift)+move;
        rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>((T)pi/2,TV(1,0,0));
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Mass(rigid_body->Mass()*10); // NOTE: was rigid_body->mass*=10
        rigid_body->name=LOG::sprintf("tread%d",tread_num);tread_num++;
        treads[tread_side].Append(rigid_body);}
    // end 2
    for(int k=4;k>=0;k--){
        rigid_body=&tests.Add_Rigid_Body("ARB/tread3_subdivided",(T).5,(T)1);
        rigid_body->Frame().t=TV(start-radius2*cos(angle_start-k*angle),initial_height+radius2*sin(angle_start-k*angle),z_shift)+move;
        rigid_body->Frame().r=ROTATION<TV>((k+(T).5)*angle,TV(0,0,1))*ROTATION<TV>((T)pi/2,TV(0,1,0))*ROTATION<TV>(-(T)pi/2,TV(1,0,0));
        rigid_body->Set_Coefficient_Of_Restitution(0); 
        rigid_body->Set_Mass(rigid_body->Mass()*10); // NOTE: was rigid_body->mass*=10
        rigid_body->name=LOG::sprintf("tread%d",tread_num);tread_num++;
        treads[tread_side].Append(rigid_body);}
    treads[tread_side].Append(treads[tread_side](1));

    // cogs - num_joints+3
    rigid_body=&tests.Add_Rigid_Body("ARB/gear3",(T).55,(T)1);
    rigid_body->Frame().t=TV(start,initial_height,z_shift)+move;
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Inertia_Tensor()*=10;
    rigid_body->name="wheel1";
    gears[3*tread_side]=rigid_body;

    rigid_body=&tests.Add_Rigid_Body("ARB/gear3",(T).55,(T)1);
    rigid_body->Frame().t=TV(0,initial_height,z_shift)+move;
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="wheel2";
    gears[3*tread_side+1]=rigid_body;

    rigid_body=&tests.Add_Rigid_Body("ARB/gear3",(T).55,(T)1);
    rigid_body->Frame().t=TV(-start,initial_height,z_shift)+move;
    rigid_body->Frame().r=ROTATION<TV>((T)pi/2,TV(0,1,0));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Inertia_Tensor()*=10;
    rigid_body->name="wheel3";
    gears[3*tread_side+2]=rigid_body;
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    Reset_Joints();
    T running_time=max(time-start_rolling,(T)0);

    T current_speed=speed,accumulated_speed=(running_time-(T)0.5*accelerate_time)*speed;
    if(running_time<accelerate_time){
        current_speed*=running_time/accelerate_time;
        accumulated_speed=(T)0.5*running_time*running_time*speed/accelerate_time;}
    TV angular_velocity=TV(0,0,-current_speed);
    ROTATION<TV> orientation=ROTATION<TV>::From_Rotation_Vector(TV(0,0,-accumulated_speed))*ROTATION<TV>((T)pi/2,TV(0,1,0));
    for(int i=0;i<6;i++){//if(i==1 || i==4) continue;
        if(turn_in_place && i==3){
            orientation=ROTATION<TV>::From_Rotation_Vector(TV(0,0,accumulated_speed))*ROTATION<TV>((T)pi/2,TV(0,1,0));
            angular_velocity=TV(0,0,current_speed);}
        gears[i]->Twist().angular=angular_velocity;
        gears[i]->Update_Angular_Momentum();}
}
//#####################################################################
// Function Create_Static_Joints
//#####################################################################
void Create_Static_Joints()
{
    // first tread
    int num_joints=11;
    int total_joints=num_joints*2+12;
    T split=(T).39;
    T start=-(split*(num_joints+1))*(T).5;

    for(int i=0;i<total_joints;i++){JOINT<TV> *left_joint=new POINT_JOINT<TV>(),*right_joint=new POINT_JOINT<TV>();
        arb->joint_mesh.Add_Articulation(treads[0](i)->particle_index,treads[0](i+1)->particle_index,left_joint);
        arb->joint_mesh.Add_Articulation(treads[0](i)->particle_index,treads[0](i+1)->particle_index,right_joint);
        left_joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T)-.425,split/2,(T)-.05)));
        left_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)-.425,-split/2,(T)-.05)));
        right_joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T).425,split/2,(T)-.05)));
        right_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).425,-split/2,(T)-.05)));}

    for(int i=0;i<total_joints;i++){JOINT<TV> *left_joint=new POINT_JOINT<TV>(),*right_joint=new POINT_JOINT<TV>();
        arb->joint_mesh.Add_Articulation(treads[1](i)->particle_index,treads[1](i+1)->particle_index,left_joint);
        arb->joint_mesh.Add_Articulation(treads[1](i)->particle_index,treads[1](i+1)->particle_index,right_joint);
        left_joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T)-.425,split/2,(T)-.05)));
        left_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)-.425,-split/2,(T)-.05)));
        right_joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T).425,split/2,(T)-.05)));
        right_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).425,-split/2,(T)-.05)));}

    // make center
    for(int i=0;i<6;i++){JOINT<TV>* joint=new POINT_JOINT<TV>();int side=i<3?1:-1;
        arb->joint_mesh.Add_Articulation(body->particle_index,gears[i]->particle_index,joint);
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T).62*side,0,0)));
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(start*(1-i%3),(T)-.25,(T)2*side),ROTATION<TV>((T)pi/2,TV(0,1,0))));}

    JOINT<TV>* body_to_gun=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(body->particle_index,gun->particle_index,body_to_gun);
    body_to_gun->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,(T)-.524,0)));
    body_to_gun->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).1,(T).55,0)));

    JOINT<TV>* gun_to_lid=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(gun->particle_index,lid->particle_index,gun_to_lid);
    gun_to_lid->Set_Joint_To_Child_Frame(FRAME<TV>(TV((T)-.575,(T)-.05,0)));
    gun_to_lid->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)-.8,(T).675,0))); 
}
//#####################################################################
// Function Reset_Joints
//#####################################################################
void Reset_Joints()
{
    for(int i=0;i<dynamic_joints.m;i++) arb->joint_mesh.Remove_Articulation(dynamic_joints(i)->id_number);
    dynamic_joints.Remove_All();

    RIGID_BODY<TV>* joint_gears[]={gears[2],gears[0],gears[5],gears[3]};
    for(int s=0;s<=1;s++) for(int i=0;i<34;i++){int index=34*s+i;
        int gear_index=next_gear[index];
        RIGID_BODY<TV> *gear=joint_gears[gear_index],*tread=treads[s](i+1);
        TV gear_axis=gear->Frame().r.Rotated_X_Axis();
        TV up=body->Frame().r.Rotated_Y_Axis();up-=TV::Dot_Product(up,gear_axis)*gear_axis;
        TV forward=TV::Cross_Product(gear_axis,up);
        TV offset=tread->Frame().t-gear->Frame().t;
        T angle_x=TV::Dot_Product(forward,offset);
        T angle_y=TV::Dot_Product(up,offset);
        T angle=atan2(angle_y,angle_x);if(!(gear_index&1)) angle+=(T)pi;angle=fmod(angle+(T)pi*2,(T)pi*2);
        if(angle>(T)1.5*(T)pi || offset.Magnitude_Squared()>sqr(2*gear_radius)) continue;
        if(angle<(T).7*(T)pi){int leading_tread=34*s+(i+1)%34;next_slot[index]=(next_slot[leading_tread]+9)%10;next_gear[index]=next_gear[leading_tread];continue;}
        JOINT<TV>* joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(gear->particle_index,tread->particle_index,joint);
        joint->Set_Joint_To_Parent_Frame(parent_to_joint[next_slot[index]]);
        joint->Set_Joint_To_Child_Frame(FRAME<TV>());
        dynamic_joints.Append(joint);}
}
//#####################################################################
// Function Reset_Joints
//#####################################################################
void Initialize_Tread_Tracking()
{
    T radius=(T).65;
    next_slot[16]=2;next_slot[33]=7;next_slot[50]=2;next_slot[67]=7;
    for(int i=0;i<4;i++) for(int j=15;j>=0;j--){int index=17*i+j;next_slot[index]=(next_slot[index+1]+9)%10;}
    for(int i=0;i<4;i++) for(int j=0;j<17;j++) next_gear[17*i+j]=i;
    for(int i=0;i<10;i++){ROTATION<TV> rotation((T)pi*2*i/10,TV(1,0,0));parent_to_joint[i]=FRAME<TV>(rotation.Rotated_Z_Axis()*radius,rotation);}
}
//#####################################################################
};
}
#endif
