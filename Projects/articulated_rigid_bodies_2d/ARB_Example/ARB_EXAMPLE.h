//#####################################################################
// Copyright 2004-2007, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class ARB_EXAMPLE
//##################################################################### 
#ifndef __ARB_EXAMPLE__
#define __ARB_EXAMPLE__

#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
namespace PhysBAM{

template<class T_input,class RW>
class ARB_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
public:
    typedef VECTOR<T,2> TV;typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solid_body_collection;using BASE::solids_evolution;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::test_number;

    ARTICULATED_RIGID_BODY<TV>* arb;
    int id1,id2,id3;
    PARAMETER_LIST parameter_list;
    SOLIDS_STANDARD_TESTS<TV> tests;
    bool add_ground;
    RIGID_BODY<TV> *ground;

    ARB_EXAMPLE(const STREAM_TYPE stream_type,const int test_number_input=1)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,output_directory,data_directory,solid_body_collection),add_ground(true),ground(0)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;

        last_frame=960;
        frame_rate=24;
        std::cout<<"Frame rate: "<<frame_rate<<std::endl;
    }

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("ARB_Example/output_%d",test_number);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;

    arb->Set_Iterative_Tolerance((T)1e-4);
    arb->Set_Contact_Level_Iterations(10);
    arb->Set_Shock_Propagation_Level_Iterations(10);
    arb->Set_Poststabilization_Iterations(10);
    arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);

    ARB_PARAMETERS::Read_Common_Parameters("ARB_Example/example.param",*this,parameter_list);

    switch(test_number){
      case 1:Five_Blocks();break;
      case 2:Dangling();break;
      case 3:Two_Tumbling();break;
      case 4:Point_Constraint_With_2_Blocks();break;
      case 5:Chain();break;
      case 6:Spring();break;
      default:PHYSBAM_FATAL_ERROR();}

    if(add_ground) ground=&tests.Add_Ground(1,0,(T).5);

    tests.Add_Gravity();
    solid_body_collection.Update_Simulated_Particles();

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();
}
//#####################################################################
// Function Two_Tumbling
//#####################################################################
void Two_Tumbling()
{
    JOINT<TV>* joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
//    JOINT<TV>* joint=new JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
//    joint->Set_Target_Rotation_Angle(2.5);
//    joint->Set_Proportional_Constant(.50);
//    joint->Set_Derivative_Constant(.50);
    
    RIGID_BODY<TV> *rigid_body1=0,*rigid_body2=0;

    T scale_factor=1;

    rigid_body1=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
    rigid_body1->Frame().t=TV(0,4);
    //rigid_body->Frame().r=ROTATION<TV>::From_Rotation_Vector(-(T)pi/4);
    //rigid_body->Angular_Momentum()=-10;
    rigid_body1->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body1->coefficient_of_friction=0;
    rigid_body1->name="square1";
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-1,1)));
    
    rigid_body2=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
    rigid_body2->Frame().t=TV(-2,6);
    rigid_body2->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body2->coefficient_of_friction=0;
    rigid_body2->name="square2";
    rigid_body2->is_static=true;
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1)));

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*rigid_body1,*rigid_body2);
    joint->Set_Joint_Function(jfunc);
    jfunc->Set_Target_Angle(ROTATION<TV>::From_Angle(-(T).3));
    joint->joint_function->Set_k_p(10);

    arb->joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
}
//#####################################################################
// Function Chain
//#####################################################################
void Chain()
{
    add_ground=false;
    FRAME<TV> xform(TV(0,3),ROTATION<TV>::From_Angle((T)0));
    int n=getenv("CHAIN")?atoi(getenv("CHAIN")):7;
    ARRAY<RIGID_BODY<TV>*> rigid_body(n);
    bool all_rigid=parameter_list.Get_Parameter("chain_all_rigid",false);
    for(int i=0;i<n;i++){
        rigid_body(i)=&tests.Add_Rigid_Body("ground",(T).01,1,false);// read only the segmented curves, that way don't get collisions between bodies
        rigid_body(i)->Frame().t=TV((T)(2*(i-1)+1),0);
        rigid_body(i)->Frame()=xform*rigid_body(i)->Frame();
        rigid_body(i)->Set_Mass(1);
        rigid_body(i)->Set_Coefficient_Of_Restitution((T).5);
        rigid_body(i)->coefficient_of_friction=1;
        rigid_body(i)->name=STRING_UTILITIES::string_sprintf("link%d",i);
        rigid_body(i)->is_static=(i==1);
        if(i>1){
            JOINT<TV>* joint=0;
            if(all_rigid || i<=1) joint=new RIGID_JOINT<TV>();
            else{
                joint=new POINT_JOINT<TV>();
#if 0
                joint->Set_Joint_Function(new JOINT_FUNCTION<TV>(joint,rigid_body(i-1),rigid_body(i)));
                joint->joint_function->Set_k_p(100);
                joint->joint_function->Set_Target_Angle(0);
#endif
            }
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0)));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,0)));
            arb->joint_mesh.Add_Articulation(rigid_body(i-1)->particle_index,rigid_body(i)->particle_index,joint);
        }
    }

#if 0
    rigid_body=&tests.Add_Rigid_Body("circle",2,(T).5);// read only the segmented curves, that way don't get collisions between bodies
    rigid_body->Set_Coefficient_Of_Restitution((T).5);
    rigid_body->coefficient_of_friction=1;
    rigid_body->Frame().t=TV(5,5);
    rigid_body->Frame()=xform*rigid_body->Frame();
    rigid_body->Set_Mass(10);
    rigid_body->name="circle");
#endif
}
//#####################################################################
// Function Spring
//#####################################################################
void Spring()
{
    T width=(T).2,vertical_separation=(T).05;
    T slanted_length=sqrt(sqr(width)+sqr(vertical_separation));
//    FRAME<TV> xform(TV(0,1),ROTATION<TV>::From_Rotation_Vector(-(T)pi/2-(T).2));
    FRAME<TV> xform(TV(0,(T).05));
//    FRAME<TV> xform;
    RIGID_BODY<TV>* rigid_body[11]={0};
    for(int i=0;i<10;i++){
        bool horizontal=((i-1)%2==0);int vertical_part=(i-1)/2;
        rigid_body[i]=&tests.Add_Rigid_Body("ground",(T).005*(horizontal?width:slanted_length),1);
        rigid_body[i]->Frame().t=TV(0,(vertical_part+(horizontal?0:(T).5))*vertical_separation);
        if(!horizontal) rigid_body[i]->Frame().r=ROTATION<TV>::From_Complex(COMPLEX<T>(width,vertical_separation));
        rigid_body[i]->Frame()=xform*rigid_body[i]->Frame();
        rigid_body[i]->Set_Coefficient_Of_Restitution(0);
        rigid_body[i]->coefficient_of_friction=1;
        rigid_body[i]->Set_Mass((horizontal?width:slanted_length)*5);
        rigid_body[i]->name=STRING_UTILITIES::string_sprintf("link%d",i);
        if(i>1){
            POINT_JOINT<TV>* joint=new POINT_JOINT<TV>();
            if(horizontal){joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(slanted_length/2,0)));joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(width/2,0)));}
            else{joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-width/2,0)));joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-slanted_length/2,0)));}
            JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*rigid_body[i-1],*rigid_body[i]);
            jfunc->Set_Target_Angle(jfunc->Angle());jfunc->Set_k_p(100);
            joint->Set_Joint_Function(jfunc);
            arb->joint_mesh.Add_Articulation(rigid_body[i-1]->particle_index,rigid_body[i]->particle_index,joint);}}
}
//#####################################################################
// Function Point_Constraint_With_2_Blocks
//#####################################################################
void Point_Constraint_With_2_Blocks()
{
    RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square_refined",1,1);
    rigid_body1->Frame().t=TV(0,16);
    rigid_body1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body1->coefficient_of_friction=1;
    rigid_body1->name="square1";
    
    RIGID_BODY<TV>* rigid_body2=&tests.Add_Rigid_Body("square_refined",1,1);
    rigid_body2->Frame().t=TV(2,18);
    rigid_body2->Twist().angular=VECTOR<T,1>((T)4);
    rigid_body2->Update_Angular_Momentum();
    rigid_body2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body2->coefficient_of_friction=1;
    rigid_body2->name="square2";

    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1)));

    arb->joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint);
}
//#####################################################################
// Function Five_Blocks
//#####################################################################
void Five_Blocks()
{
    T scale_factor=1;

    RIGID_BODY<TV>* rigid_body1=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
    rigid_body1->Frame().t=TV(-1,3);
    rigid_body1->Set_Coefficient_Of_Restitution((T).5);
    rigid_body1->coefficient_of_friction=0;
    rigid_body1->name="square1";
    
    RIGID_BODY<TV>* rigid_body2=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
    rigid_body2->Frame().t=TV(0,2);
    rigid_body2->Twist().linear=TV();
    rigid_body2->Set_Coefficient_Of_Restitution((T).5);
    rigid_body2->coefficient_of_friction=0;
    rigid_body2->name="squarejoint";
    
    RIGID_BODY<TV>* rigid_body3=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
    rigid_body3->Frame().t=TV(1,3);
    rigid_body3->Twist().linear=TV();
    rigid_body3->Set_Coefficient_Of_Restitution((T).5);
    rigid_body3->coefficient_of_friction=0;
    rigid_body3->name="square2";

    JOINT<TV> *joint1=new RIGID_JOINT<TV>(),*joint2=new RIGID_JOINT<TV>();
    joint1->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1)));
    joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1)));
    arb->joint_mesh.Add_Articulation(rigid_body1->particle_index,rigid_body2->particle_index,joint1);
    arb->joint_mesh.Add_Articulation(rigid_body2->particle_index,rigid_body3->particle_index,joint2);
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE{
#if 0
    collision_manager.hash.Insert(PAIR<int,int>(id1,id2));
    collision_manager.hash.Insert(PAIR<int,int>(id2,id1));
#endif
    // collisions
    switch(test_number){
      case 1:{
            if(!solids_evolution->rigid_body_collisions->collision_manager){
                solids_evolution->rigid_body_collisions->collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;}
            RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager=dynamic_cast<RIGID_BODY_COLLISION_MANAGER_HASH*>(solids_evolution->rigid_body_collisions->collision_manager);
            assert(collision_manager);
            for(int i=0;i<solid_body_collection.rigid_body_collection.rigid_body_particles.Size();i++) 
                if(i==ground->particle_index){
                    collision_manager->hash.Insert(PAIR<int,int>(i,ground->particle_index));
                    collision_manager->hash.Insert(PAIR<int,int>(ground->particle_index,i));}}
      default:break;}
}
//#####################################################################
// Function Dangling
//#####################################################################
void Dangling()
{
    bool half_chain=parameter_list.Get_Parameter("dangling_half_chain",false);
    add_ground=parameter_list.Get_Parameter("add_ground",false);

    T epsilon=(T).5;

    JOINT<TV>* joint[11];
    RIGID_BODY<TV> *rigid_body[12]={0};
    T scale_factor=(T).65;
    T shift=sqrt((T).5);

    for(int i=0;i<12-half_chain*6;i++){
        rigid_body[i]=&tests.Add_Rigid_Body("square_refined",scale_factor,(T).5);
        rigid_body[1]->Set_Coefficient_Of_Restitution(epsilon);
        rigid_body[i]->Set_Mass(10);
        if(i>=1) joint[i-1]=new POINT_JOINT<TV>();
        if(i>=2){
            joint[i-2]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,0)));
            if(i!=6 || !half_chain) joint[1-1]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,0)));}}

    rigid_body[0]->Frame().t=TV(+2+9*shift-(T).25,9-8*shift);
    rigid_body[0]->name="static";
    rigid_body[0]->is_static=true;
    rigid_body[1]->Frame().t=TV(-2-9*shift+(T).25,9-8*shift);
    rigid_body[1]->name="static";
    rigid_body[1]->is_static=true;

    for(int i=2;i<=5;i++){
        rigid_body[i]->Frame().r=ROTATION<TV>::From_Angle((T)pi/4);
        rigid_body[i]->Frame().t=TV(-2-(11-2*i)*shift,8-(11-2*i)*shift);}

    rigid_body[6]->Frame().t=TV(-1,8);
    rigid_body[7]->Frame().t=TV(1,8);

    for(int i=8;i<=11;i++){
        rigid_body[i]->Frame().r=ROTATION<TV>::From_Angle(-(T)pi/4);
        rigid_body[i]->Frame().t=TV(2+(2*i-15)*shift,8-(2*i-15)*shift);}

    joint[0]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,-1)));
    if(!half_chain) joint[10]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-1,-1)));

    for(int i=1;i<12-half_chain*6;i++){
        arb->joint_mesh.Add_Articulation(rigid_body[i]->particle_index,rigid_body[i+1]->particle_index,joint[i-1]);}

    if(!half_chain) arb->joint_mesh.Add_Articulation(rigid_body[11]->particle_index,rigid_body[0]->particle_index,joint[10]);
}
//#####################################################################
};
}
#endif
