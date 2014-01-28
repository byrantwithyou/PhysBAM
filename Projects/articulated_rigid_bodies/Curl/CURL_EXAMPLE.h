//#####################################################################
// Copyright 2004-2008, Craig Schroeder, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CURL_EXAMPLE
//#####################################################################
#ifndef __CURL_EXAMPLE__
#define __CURL_EXAMPLE__

#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Math_Tools/wrap.h>
#include <Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <Rigids/Joints/ANGLE_JOINT.h>
#include <Rigids/Joints/JOINT_FUNCTION.h>
#include <Rigids/Joints/POINT_JOINT.h>
#include <Rigids/Joints/RIGID_JOINT.h>
#include <Rigids/Particles/RIGID_BODY_PARTICLES.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_STATE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include "../ARB_PARAMETERS.h"
namespace PhysBAM{

template<class T_input>
class CURL_EXAMPLE:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::stream_type;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solid_body_collection;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::parse_args;
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Position_Nodes; // silence -Woverloaded-virtual

    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    RIGID_BODY<TV>* shelf00,*shelf01,*shelf10,*shelf11;
    int current_frame,start_move,end_move;
    T increment;
    int selection;
    bool add_ground;
    PARAMETER_LIST parameter_list;
    bool traditional_pd;
    ARRAY<TV,int> precomputed_pd_torques;
    bool half_acceleration;
    RIGID_BODY<TV>* jitter_body;
    TV max_jitter_body_move;
    std::string parameter_file;

    CURL_EXAMPLE(const STREAM_TYPE stream_type,std::string parameter_file_input="")
        :BASE(stream_type),tests(stream_type,data_directory,solid_body_collection),parameter_file(parameter_file_input)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;

        last_frame=2000;
        frame_rate=60;
        LOG::cout<<"Frame rate: "<<frame_rate<<std::endl;

        current_frame=0;
        increment=(T).05;
        start_move=5;end_move=40;
        shelf00=shelf01=shelf10=shelf11=0;
        write_last_frame=true;
    }

    virtual ~CURL_EXAMPLE()
    {
        delete arb;
    }

    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

    void Register_Options() PHYSBAM_OVERRIDE
    {
        BASE::Register_Options();
    }
    void Parse_Options() PHYSBAM_OVERRIDE
    {
        BASE::Parse_Options();
        tests.data_directory=data_directory;
        output_directory="Curl/output";
    }
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    add_ground=parameter_list.Get_Parameter("add_ground",false);
    RIGID_BODY<TV>* rigid_body=0;
    int num_joints=0,num_rigid_bodies=0;

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Iterative_Tolerance((T)1e-5);
    arb->Set_Contact_Level_Iterations(100);
    arb->Set_Shock_Propagation_Level_Iterations(100);
    arb->Set_Poststabilization_Iterations(100);
    arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);
    arb->Use_PD_Actuators();

    if(parameter_file.empty()) parameter_file="Curl/example.param";
    ARB_PARAMETERS::Read_Common_Parameters(parameter_file,*this,parameter_list);
    selection=parameter_list.Get_Parameter("selection",5);
    output_directory+=STRING_UTILITIES::string_sprintf("_%d",selection);
    traditional_pd=parameter_list.Get_Parameter("traditional_pd",false);
    if(traditional_pd) arb->Use_No_Actuators();
    half_acceleration=parameter_list.Get_Parameter("half_acceleration",true);
    if(half_acceleration) LOG::cout<<"USING HALF ACCEL"<<std::endl;

    switch(selection){
        case 1:
          PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,0,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p",(T)25));
          break;
        case 2: // point constraint
          rigid_body=&tests.Add_Rigid_Body("sphere",(T).5,(T).5);
          rigid_body->Frame().t=parameter_list.Get_Parameter("ball_position",TV((T)1.5,4,0));
          rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
          rigid_body->name="ball";
          rigid_body->Set_Mass(parameter_list.Get_Parameter("ball_mass",(T)1));
          rigid_body->Twist().linear=parameter_list.Get_Parameter("ball_velocity",TV());
          rigid_body->Angular_Momentum()=TV(0,0,parameter_list.Get_Parameter("ball_momentum",(T)0));
          num_rigid_bodies+=1;

          PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,0,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p_1",(T)5));
          if(parameter_list.Get_Parameter("plank_2",true))
              PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(12,-5,0),ROTATION<TV>((T)pi,TV(0,1,0)),parameter_list.Get_Parameter("k_p_2",(T)15));
          if(parameter_list.Get_Parameter("plank_3",true))
              PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,-9,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p_3",(T)25));
          if(parameter_list.Get_Parameter("plank_4",true))
              PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(12,-12,0),ROTATION<TV>((T)pi,TV(0,1,0)),parameter_list.Get_Parameter("k_p_4",(T)50));
          if(parameter_list.Get_Parameter("plank_5",true))
              PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,-14,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p_5",(T)100));
          break;
      case 3:
          rigid_body=&tests.Add_Rigid_Body("subdivided_box",(T).75,(T)1);
          rigid_body->Frame().t=TV(parameter_list.Get_Parameter("block_x",(T)8),(T).875,0);
          rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
          rigid_body->name="block";
          rigid_body->Set_Mass(parameter_list.Get_Parameter("block_mass",(T)1));
          num_rigid_bodies+=1;

          PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,0,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p",(T)100));

          if(parameter_list.Get_Parameter("apply_block_joint",false)){
              int plank_num(6);
              JOINT<TV>* joint=new POINT_JOINT<TV>();
              arb->joint_mesh.Add_Articulation(plank_num,rigid_body->particle_index,joint);
              joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-(T)0.25,(T)0.125,(T)0.1875)));
              joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)0.75,-(T)0.75,(T)0.1875)));
              JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,arb->rigid_body_collection.Rigid_Body(plank_num),*rigid_body);
              joint->Set_Joint_Function(jfunc);
              joint->joint_function->Set_k_p(parameter_list.Get_Parameter("k_p",(T)100));
              joint->joint_function->Set_Target_Angle(ROTATION<TV>());

              JOINT<TV>* joint2=new POINT_JOINT<TV>();
              arb->joint_mesh.Add_Articulation(plank_num,rigid_body->particle_index,joint2);
              joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-(T)0.25,(T)0.125,-(T)0.1875)));
              joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)0.75,-(T)0.75,-(T)0.1875)));
              JOINT_FUNCTION<TV>* jfunc2=new JOINT_FUNCTION<TV>(*joint,arb->rigid_body_collection.Rigid_Body(plank_num),*rigid_body);
              joint2->Set_Joint_Function(jfunc2);
              joint2->joint_function->Set_k_p(parameter_list.Get_Parameter("k_p",(T)100));
              joint2->joint_function->Set_Target_Angle(ROTATION<TV>());
          }
          break;
      case 4:
          PD_Plank_Curl_Test(num_joints,num_rigid_bodies,TV(0,0,0),ROTATION<TV>(),parameter_list.Get_Parameter("k_p",(T)25));

          rigid_body=&tests.Add_Rigid_Body("miniplank25wide2",parameter_list.Get_Parameter("wall_scale",(T)2),(T).5);
          rigid_body->Frame().t=parameter_list.Get_Parameter("wall_position",TV());
          rigid_body->Frame().r=ROTATION<TV>((T)pi/2,parameter_list.Get_Parameter("wall_rotation_vector",TV(0,0,1)));
          rigid_body->Set_Coefficient_Of_Restitution(1);
          rigid_body->name="wall";
          rigid_body->is_static=true;

        break;
      case 5:
        Raise_To_Horizontal(num_joints,num_rigid_bodies,parameter_list.Get_Parameter("num_links",1),parameter_list.Get_Parameter("k_p",(T)25));
        break;
      case 6:
            PD_Plank_Test(num_joints,num_rigid_bodies,TV(0,0,0),ROTATION<TV>());
            break;
      case 7:
        Jitter_Test(num_joints,num_rigid_bodies,parameter_list.Get_Parameter("num_links",3));
        break;
      default:
          LOG::cout<<"make bodies for testing"<<std::endl;
          rigid_body=&tests.Add_Rigid_Body("short_plank_subdivided",1,(T).5);
          rigid_body->Frame().t=TV(0,0,0);
          rigid_body->Frame().r=ROTATION<TV>(-(T)pi/2,TV(0,1,0));
          rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
          rigid_body->name="parent";
          rigid_body->is_static=true;
          shelf00=rigid_body;

          rigid_body=&tests.Add_Rigid_Body("subdivided_box",1,(T).5);
          rigid_body->Frame().t=TV(1,10,0);
          rigid_body->Twist().linear=TV(0,0,0);
          rigid_body->Angular_Momentum()=TV(2,3,5);
          rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
          rigid_body->name="child";
          break;
    }

    if(add_ground){
        RIGID_BODY<TV>& ground=tests.Add_Ground((T).5,-10);
        ground.Set_Coefficient_Of_Restitution((T)1.0);
        ground.name="ground";}

    tests.Add_Gravity();

    LOG::cout<<"done initializing example"<<std::endl;
}
//#####################################################################
// PD Plank Curl test
//#####################################################################
void PD_Plank_Curl_Test(int& num_joints,int& num_rigid_bodies,TV shift,ROTATION<TV> orient,const T k_p)
{
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>* parent_body=0,*child_body=0;
    T cheight=(T)0;

    T desired_x=parameter_list.Get_Parameter("desired_x",(T)0);
    // Create first body
    parent_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T)1);
    parent_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
    parent_body->Frame().r=orient;
    //parent_body->Angular_Momentum()=TV(0,0,2);
    parent_body->Set_Coefficient_Of_Restitution((T)0.5);
    parent_body->name="parent";
    parent_body->Set_Mass(parameter_list.Get_Parameter("parent_mass",(T)5));
    parent_body->is_static=parameter_list.Get_Parameter("pd_curl_test_static",(bool)true);
    num_rigid_bodies+=1;

    // Add children and joint
    int njoints=parameter_list.Get_Parameter("pd_curl_test_njoints",6);
    if(parameter_list.Get_Parameter("create_polygon",false)) desired_x=(2*(T)pi)/((T)(njoints+1));
    for(int i=0;i<njoints;i++){
        cheight+=(T)1.25;
        child_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T).5);
        child_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Frame().r=orient;
        child_body->Twist().linear=TV(0,0,0);
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);

        if(parameter_list.Get_Parameter("use_bend_joint",false)) joint=new ANGLE_JOINT<TV>();
        else joint=new POINT_JOINT<TV>();
        arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint);

        // Create joint function
        LOG::cout<<"Desired rotation: "<<desired_x<<std::endl;
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x,TV(1,0,0));

        JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*parent_body,*child_body);
        joint->Set_Joint_Function(jfunc);
        joint->joint_function->Set_k_p(k_p);
        joint->joint_function->Set_Target_Angle(desired_rotation);

        if(parameter_list.Get_Parameter("use_double_point_joints",false)){
            JOINT<TV>* joint2=new POINT_JOINT<TV>();
            arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint2);
            joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            num_joints+=1;
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        else{
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        num_joints+=1;
        num_rigid_bodies+=1;

        // Swap now
        parent_body=child_body;
        child_body=0;
        joint=0;
    }
    for(int i=0;i<arb->joint_mesh.Num_Joints();i++){
        LOG::cout<<"Desired rotation for joint "<<i<<" is "<<arb->joint_mesh.Joints(i)->joint_function->target_angle<<std::endl;}

    LOG::cout<<"initializing point joint example"<<std::endl;
}
//#####################################################################
// Raise To Horizontal
//#####################################################################
void Raise_To_Horizontal(int& num_joints,int& num_rigid_bodies,const int num_links,const T k_p)
{
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>* parent_body=0,*child_body=0;
    T cheight=(T)0;

    T desired_x=parameter_list.Get_Parameter("desired_x",(T)0);
    // Create first body
    parent_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T)1);
    parent_body->Frame().t=TV(-(T).625,-(T).625,0);
    parent_body->Set_Coefficient_Of_Restitution((T)0.5);
    parent_body->coefficient_of_friction=1;
    parent_body->name="parent";
    parent_body->Set_Mass(parameter_list.Get_Parameter("parent_mass",(T)5));
    parent_body->is_static=parameter_list.Get_Parameter("pd_curl_test_static",(bool)true);
    num_rigid_bodies+=1;

    // Add children and joint
    for(int i=0;i<num_links;i++){
        cheight-=(T)1.25;
        child_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T).5);
        child_body->Frame().t=TV(0,cheight,0);
        child_body->Frame().r=ROTATION<TV>(-(T)pi/2,TV(0,0,1));
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);

        joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint);

        // Create joint function
        LOG::cout<<"Desired rotation: "<<desired_x<<std::endl;
        ROTATION<TV> desired_rotation=parameter_list.Get_Parameter("desired_rotation",ROTATION<TV>(0,TV(1,0,0)));
        JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*parent_body,*child_body);
        joint->Set_Joint_Function(jfunc);
        joint->joint_function->Set_k_p(k_p);
        joint->joint_function->Set_Target_Angle(desired_rotation);

        if(parameter_list.Get_Parameter("use_double_point_joints",false)){
            JOINT<TV>* joint2=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint2);
            joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            num_joints+=1;
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        else{
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        num_joints+=1;
        num_rigid_bodies+=1;

        // Swap now
        parent_body=child_body;
        child_body=0;
        joint=0;
    }

    LOG::cout<<"initializing point joint example"<<std::endl;
}//#####################################################################
// Raise To Horizontal
//#####################################################################
void Jitter_Test(int& num_joints,int& num_rigid_bodies,const int num_links)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    POINT_JOINT<TV>* joint=0;
    RIGID_BODY<TV>* parent_body=0,*child_body=0;
    T cheight=(T)0;

    T desired_x=parameter_list.Get_Parameter("desired_x",(T)0);
    // Create first body
    parent_body=&tests.Add_Rigid_Body("subdivided_box",1,(T)1,false);
    parent_body->Frame().t=TV(0,0,0);
    parent_body->Frame().r=ROTATION<TV>::From_Components((T).5,-(T).5,-(T).5,-(T).5);
    parent_body->Set_Coefficient_Of_Restitution((T)0.5);
    parent_body->Twist().linear=TV((T).1,0,0);
    //parent_body->Angular_Momentum()=TV(0,0,2);
    parent_body->name="parent";
    parent_body->Set_Mass(parameter_list.Get_Parameter("parent_mass",(T).001));
   //parent_body->is_static=parameter_list.Get_Parameter("pd_curl_test_static",(bool)true);
    rigid_body_collection.rigid_body_particles.kinematic(parent_body->particle_index)=true;
    num_rigid_bodies+=1;
    jitter_body=parent_body;

    // Add children and joint
    for(int i=0;i<num_links;i++){
        cheight-=(T)2.5;
        child_body=&tests.Add_Rigid_Body("subdivided_box",1,(T).5,false);
        child_body->Frame().t=TV(0,cheight,0);
        child_body->Frame().r=ROTATION<TV>::From_Components((T).5,-(T).5,-(T).5,-(T).5);
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);
        child_body->Set_Mass((T).001);

        joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint);
        joint->angular_damping=parameter_list.Get_Parameter("joint_damping_value",(T)0)*i;

        // Create joint function
        LOG::cout<<"Desired rotation: "<<desired_x<<std::endl;
        ROTATION<TV> desired_rotation=parameter_list.Get_Parameter("desired_rotation",ROTATION<TV>(0,TV(1,0,0)));
        JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*parent_body,*child_body);
        joint->Set_Joint_Function(jfunc);
        joint->joint_function->Set_k_p(parameter_list.Get_Parameter("k_p",(T)5));
        joint->joint_function->Set_Target_Angle(desired_rotation);

        if(parameter_list.Get_Parameter("use_double_point_joints",false)){
            JOINT<TV>* joint2=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint2);
            joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,(T)1.25,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-(T)1.25,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            num_joints+=1;
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,(T)1.25,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-(T)1.25,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        else{
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,-(T)1.25),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,0,(T)1.25),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        num_joints+=1;
        num_rigid_bodies+=1;

        // Swap now
        parent_body=child_body;
        child_body=0;
        joint=0;
    }

    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(i);
        if(rigid_body_collection.rigid_body_particles.kinematic(i)){FRAME<TV> frame;Set_Kinematic_Positions(frame,0,i);rigid_body.Frame()=frame;}}

    LOG::cout<<"initializing point joint example"<<std::endl;
}
//#####################################################################
// PD Plank
//#####################################################################
void PD_Plank_Test(int& num_joints,int& num_rigid_bodies,const TV& shift,const ROTATION<TV>& orient)
{
    JOINT<TV>* joint=0;
    RIGID_BODY<TV>* parent_body=0,*child_body=0;
    T cheight=(T)0;

    // Create first body
    parent_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T)1);
    parent_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
    parent_body->Frame().r=orient;
    parent_body->Set_Coefficient_Of_Restitution((T)0.5);
    parent_body->name="parent";
    parent_body->is_static=true;
    num_rigid_bodies+=1;

    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    // Add children and joint
    int njoints=parameter_list.Get_Parameter("pd_curl_test_njoints",(int)6);
    for(int i=0;i<njoints;i++){
        cheight+=(T)1.25;
        child_body=&tests.Add_Rigid_Body("miniplank25wide2",1,(T)1);
        child_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Frame().r=orient;
        child_body->Twist().linear=TV(0,0,0);
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);

        if(parameter_list.Get_Parameter("use_point_joint",true)) joint=new POINT_JOINT<TV>();
        else joint=new ANGLE_JOINT<TV>();

        arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint);
#if 1
        // Create joint function
        ROTATION<TV> desired_rotation=ROTATION<TV>((T)0,TV(1,0,0));
        JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(*joint,*parent_body,*child_body);
        joint->Set_Joint_Function(jfunc);
        joint->joint_function->Set_k_p(k_p);
        joint->joint_function->Set_Target_Angle(desired_rotation);
#endif

        if(typeid(*joint)==typeid(POINT_JOINT<TV>) && parameter_list.Get_Parameter("use_double_point_joints",false)){
            JOINT<TV>* joint2=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_rigid_bodies),int(num_rigid_bodies+1),joint2);
            joint2->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint2->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,-(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            num_joints+=1;
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,(T).1),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        else{
            joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
            joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T).625,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
        }
        num_joints+=1;
        num_rigid_bodies+=1;

        if(i==0) joint->Set_Joint_Frame(FRAME<TV>(ROTATION<TV>(-(T)pi/2,TV(1,0,0))));

        // Swap now
        parent_body=child_body;
        child_body=0;
        joint=0;
    }

    arb->Update_With_Breadth_First_Directed_Graph(int(1));

//    jitter_body->Save_State(1,0);
//    jitter_body->Frame().t=TV(2,0,0);
//    jitter_body->Save_State(2,2);

    LOG::cout<<"initializing point joint example"<<std::endl;
}
//#####################################################################
//
//#####################################################################
/*void Add_External_Forces(TV& F,TV& torque,const T time) PHYSBAM_OVERRIDE
{
    if(traditional_pd){
        LOG::cout<<"Adding torque "<<precomputed_pd_torques(fragment_id)<<" to "<<fragment_id<<std::endl;
        torque+=precomputed_pd_torques(fragment_id);
    }
}
*/
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    INTERPOLATION_CURVE<T,TV> motion_curve;
    motion_curve.Add_Control_Point(0,TV(0,0,0));
    motion_curve.Add_Control_Point((T).833333,TV(5,0,0));
    frame=FRAME<TV>(motion_curve.Value(time),ROTATION<TV>::From_Components((T).5,-(T).5,-(T).5,(T)-.5));
}
//#####################################################################
// Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(selection==7){//jitter test
        //jitter_body->Frame().t=TV(0,(T)-.1+(T).2*sin(fmod(100*time,(T)pi)),0);
/*        if(jitter_body->Frame().t.x<2){
            jitter_body->Frame().t+=TV((T).005,0,0);
        }
*/  }
    PHYSBAM_FATAL_ERROR("pd is now done in the framework");
    if(traditional_pd){
        precomputed_pd_torques.Resize(solid_body_collection.rigid_body_collection.rigid_body_particles.Size());precomputed_pd_torques.Fill(TV());
        for(int i=0;i<arb->joint_mesh.Num_Joints();i++) if(arb->joint_mesh.Joints(i)->joint_function){
            JOINT<TV>* joint=arb->joint_mesh.Joints(i);JOINT_FUNCTION<TV>* jfunc=joint->joint_function;
            RIGID_BODY<TV>* parent=arb->Parent(joint->id_number),*child=arb->Child(joint->id_number);
            ROTATION<TV> Fp_wj=parent->Frame().r*joint->F_pj().r;
            TV rotation_axis_to_current_target=Fp_wj.Rotate((jfunc->Target_Angle(time)*jfunc->Angle().Inverse()).Rotation_Vector());
            T target_minus_current_angle=rotation_axis_to_current_target.Normalize();
            target_minus_current_angle=wrap(target_minus_current_angle,-(T)pi,(T)pi);
            // target angular velocity in parent space; angular velocity already in world space
            T target_minus_current_angular_velocity=TV::Dot_Product(Fp_wj.Rotate(jfunc->Target_Angular_Velocity(time))-jfunc->Angular_Velocity(),rotation_axis_to_current_target);

            TV target_acceleration;
            LOG::cout<<target_minus_current_angle<<", "<<target_minus_current_angular_velocity<<", "<<jfunc->k_p<<" ,"<<rotation_axis_to_current_target<<std::endl;
            TV pd_accel=(jfunc->k_p*(target_minus_current_angle)+jfunc->k_v*(target_minus_current_angular_velocity))*rotation_axis_to_current_target;
            if(half_acceleration) pd_accel*=(T).5;
            TV acceleration=target_acceleration+pd_accel;
            precomputed_pd_torques(parent->particle_index)-=acceleration;
            precomputed_pd_torques(child->particle_index)+=acceleration;}}
}
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    BASE::Write_Output_Files(frame);
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;

    TV total_linear_momentum,total_angular_momentum;
    for(int i=0;i<rigid_body_particles.Size();i++){
        LOG::cout<<i<<": "<<rigid_body_particles.angular_momentum(i)<<std::endl;
        total_linear_momentum+=rigid_body_particles.mass(i)*rigid_body_particles.twist(i).linear;
        total_angular_momentum+=TV::Cross_Product(rigid_body_particles.frame(i).t,rigid_body_particles.mass(i)*rigid_body_particles.twist(i).linear)+rigid_body_particles.angular_momentum(i);}
    LOG::cout<<"MOMENTA === linear "<<total_linear_momentum<<", angular "<<total_angular_momentum<<std::endl;
}
//#####################################################################
};
}
#endif
