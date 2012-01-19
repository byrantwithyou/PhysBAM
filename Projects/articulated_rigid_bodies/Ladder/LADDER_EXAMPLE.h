//#####################################################################
// Copyright 2005, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class LADDER_EXAMPLE
//##################################################################### 
#ifndef __LADDER_EXAMPLE__
#define __LADDER_EXAMPLE__

#include <PhysBAM_Tools/Interpolation/BSPLINE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
#include "../VISIBLE_HUMAN.h"

const float POINT_TOLERANCE=1e-6;

namespace PhysBAM{

template<class T,class RW>
class LADDER_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
public:
    typedef VECTOR<T,3> TV;
    typedef CONSTRAINED_POINT_IN_RIGID_BODY<T,TV> T_CONSTRAINED_POINT_IN_RIGID_BODY;

    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::data_directory;
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    PARAMETER_LIST parameter_list;

    bool add_ground,add_block;

    VISIBLE_HUMAN<T,RW>* da_man;
    FRAME_3D<T> skeleton_frame;
    bool preprocess_complete,test_frame_tracks;
    T ladder_rotate;
    int root;
    ARRAY<PAIR<std::string,FRAME_3D<T> > > starting_joint_positions;
    ARRAY<RIGID_BODY<TV>*> ladder_bodies;
    ARRAY<FRAME_TRACK_3D<T>*> frame_tracks;
    ARRAY<JOINT<TV>*> frame_track_joints;
    ARRAY<RIGID_BODY<TV>*> dynamic_joint_check_bodies;
    bool right_hand_joint,left_hand_joint,right_foot_joint,left_foot_joint;
    ARRAY<JOINT<TV>*> saved_joint;
    ARRAY<int> saved_parent;
    ARRAY<int> saved_child;

    LADDER_EXAMPLE(const std::string parameter_file="Ladder/example.param")
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >::NONE)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        //solids_parameters.gravity_direction=TV(0,0,1);
        //solids_parameters.gravity=0;
        //restart=true;
        //restart_frame=5;
        last_frame=2000;
        frame_rate=60;
        output_directory="Ladder/output";

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        arb->gravity=solids_parameters.gravity;
        arb->gravity_direction=solids_parameters.gravity_direction;
        
        arb->Set_Iterative_Tolerance((T)1e-4);
        //arb->Set_Extra_Iterations_Per_Contact_Level_Factor(100);
        //arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(100);
        //arb->Set_Poststabilization_Iterations(100);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);
        arb->Use_PD_Actuators();
        arb->Use_Balance_Controller(false);
        write_last_frame=true;
        add_ground=true;
        add_block=false;
        right_hand_joint=left_hand_joint=right_foot_joint=left_foot_joint=false;

        std::cout<<"reading parameter file "<<parameter_file<<std::endl;
        ARB_PARAMETERS::Read_Common_Parameters(parameter_file,*this,parameter_list);
        ladder_rotate=parameter_list.Get_Parameter("ladder_rotate",(T)0);

        preprocess_complete=false;

        test_frame_tracks=parameter_list.Get_Parameter("test_frame_tracks",false);
        if(test_frame_tracks) solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;
    }
    
    ~LADDER_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Apply_Constraints
//#####################################################################
// used in test_frame_tracks case to dump the animated body without simulating
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Apply_Constraints(dt,time);
    if(test_frame_tracks) Update_Joints(time);
}
//#####################################################################
// Update_Joints
//#####################################################################
void Update_Joints(const T time)
{
    for(int i=1;i<=arb->joint_mesh.joints.m;i++)
        if(arb->joint_mesh.joints(i)->joint_function)
            arb->joint_mesh.joints(i)->Set_Joint_Frame(FRAME_3D<T>(arb->joint_mesh.joints(i)->joint_function->Target_Angle(time)));
    arb->Update_With_Breadth_First_Directed_Graph(root);
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Save_Joint_Information()
{
    for(int i=1;i<=arb->joint_mesh.joints.m;i++){
        JOINT<TV>* joint=arb->joint_mesh.joints(i);
        saved_joint.Append(joint);
        saved_parent.Append(arb->Parent_Id(joint->id_number));
        saved_child.Append(arb->Child_Id(joint->id_number));}
}
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    if(parameter_list.Get_Parameter("use_ladder",true)){
        Make_Ladder(10,parameter_list.Get_Parameter("ladder_scale",(T)1));}

    T man_x=parameter_list.Get_Parameter("man_x",(T)0);
    T man_y=parameter_list.Get_Parameter("man_y",(T)0);
    T man_z=parameter_list.Get_Parameter("man_z",(T)0);
//ROTATION<TV>(.05,TV(0,0,1))*
    skeleton_frame=FRAME_3D<T>(TV(man_x,man_y,man_z))*FRAME_3D<T>(ROTATION<TV>(ladder_rotate-.4,TV(1,0,0)))*FRAME_3D<T>(TV(0,.05,0),ROTATION<TV>(pi,TV(0,1,0))*ROTATION<TV>((-0.5)*pi,TV(1,0,0)));
    //skeleton_frame=FRAME_3D<T>(TV(0,2.3,-0.2))*FRAME_3D<T>(ROTATION<TV>(0.3,TV(1,0,0)))*FRAME_3D<T>(TV(0,-0.12,0),ROTATION<TV>(pi,TV(0,1,0))*ROTATION<TV>(-0.5*pi,TV(1,0,0)));
    da_man=new VISIBLE_HUMAN<T,RW>(arb,data_directory,skeleton_frame);
    da_man->use_double_point_joints_for_bend=parameter_list.Get_Parameter("use_only_point_joints",false);
    da_man->Initialize_Bodies();

    root=da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_CRANIUM)->id_number;
    dynamic_joint_check_bodies.Append(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_R_PALM));
    dynamic_joint_check_bodies.Append(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_L_PALM));
    dynamic_joint_check_bodies.Append(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_R_ANKLE));
    dynamic_joint_check_bodies.Append(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_L_ANKLE));

    Climb();

    if(add_ground){
        int id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/ground");
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->frame.t=TV(0,.1,0);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;
    }

    for(int i=1;i<=arb->rigid_body_list.Number_Of_Elements();i++) arb->rigid_body_particles.Rigid_Body(i).Set_Coefficient_Of_Friction(1);

    RIGID_BODY_LIST<T,TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) {if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);}

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
    std::cout<<"done initializing example\n";

    Save_Joint_Information();
}
//#####################################################################
// Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    PHYSBAM_FATAL_ERROR("this should be done in kinematic rigid body callbacks");
    if(parameter_list.Get_Parameter("reset_velocities",true) && restart){
        for(int i=1;i<=arb->rigid_body_list.rigid_bodies.m;i++){
            arb->rigid_body_list.rigid_bodies(i)->velocity=TV();
            arb->rigid_body_list.rigid_bodies(i)->angular_velocity=TV();
        }
    }
    //da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_RADIOULNAR)->Set_Joint_Frame(FRAME_3D<T>(Get_Joint_Rotation("joint_r_radioulnar")));
    //da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_RADIOULNAR)->Set_Joint_Frame(FRAME_3D<T>(Get_Joint_Rotation("joint_r_radioulnar_left")));
}
//#####################################################################
// Write_Output_Files
//#####################################################################
/*void Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);

    const RIGID_BODY_LIST<T,TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output("handtransform_"+STRING_UTILITIES::string_sprintf("%d",frame),false);
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++){
        (*output)<<data_directory<<"/Rigid_Bodies/New_Visible_Human_Bones/";
        (*output)<<rigid_body_particles.Rigid_Body(i).name<<"\n"<<rigid_body_particles.Rigid_Body(i).frame<<std::endl;
    }
    delete output;
}*/
//#####################################################################
// Function Make Ladder
//#####################################################################
void Make_Ladder(const int number_of_rungs,const T scale_factor=(T)1)
{
    int id;
    RIGID_BODY<TV>* rigid_body;
    ladder_bodies.Resize(number_of_rungs);
    T ladder_x=parameter_list.Get_Parameter("ladder_x",(T)0);
    T ladder_y=parameter_list.Get_Parameter("ladder_y",(T)0);
    T ladder_z=parameter_list.Get_Parameter("ladder_z",(T)0);
    for(int i=0;i<number_of_rungs;i++){
        if(parameter_list.Get_Parameter("wide_ladder",true))
            id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/wideladder_stretched",scale_factor);
        else
            id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/wideladder",scale_factor);
        rigid_body=arb->rigid_body_list.rigid_bodies(id);
        ROTATION<TV> r=ROTATION<TV>(ladder_rotate,TV(1,0,0));
        TV t=r.Rotate(scale_factor*TV(0,.25+i*.5,0)+TV(ladder_x,ladder_y,ladder_z));
        rigid_body->frame=FRAME_3D<T>(t,r)*rigid_body->frame;
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ladder");
        rigid_body->is_static=true;
        ladder_bodies(i+1)=rigid_body;
    }
}
void Add_Ladder_Joint(const int location,const int ladder_id){
    std::cout<<"ADD LADDER JOINT\n";
    RIGID_BODY<TV>* ladder_body=ladder_id?ladder_bodies(ladder_id):0;
    RIGID_BODY<TV>* bone_body;
    switch(location){
      case 1: //right hand
        bone_body=dynamic_joint_check_bodies(1);
        break;
      case 2: // left hand
        bone_body=dynamic_joint_check_bodies(2);
        break;
      case 3: // right foot
        bone_body=dynamic_joint_check_bodies(3);
        break;
      case 4: // left foot
        bone_body=dynamic_joint_check_bodies(4);
        break;
      default:
        bone_body=0;
        break; // do nothing
    }
    if(bone_body && ladder_body){
        std::cout<<"Adding joint between "<<bone_body->name<<" and "<<ladder_body->name<<std::endl;
        JOINT<TV>* joint;
        if(parameter_list.Get_Parameter("use_only_point_joints",false)) joint=new POINT_JOINT<TV>();
        else joint=new ANGLE_JOINT<TV>();
        TV closest_ladder_point=bone_body->frame.t;
        float distance=ladder_body->Implicit_Geometry_Extended_Value(closest_ladder_point);
        while(fabs(distance)>POINT_TOLERANCE){
            closest_ladder_point+=-ladder_body->Implicit_Geometry_Extended_Normal(closest_ladder_point)*distance;
            distance=ladder_body->Implicit_Geometry_Extended_Value(closest_ladder_point);}
        TV ladder_space_point=ladder_body->Object_Space_Point(closest_ladder_point);
        TV joint_location(ladder_space_point.x,0,0);
        std::cout<<"joint_location: "<<joint_location<<" and lcosest point: "<<closest_ladder_point<<std::endl;
        joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(joint_location));
        joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(TV((location==1||location==2)?-.03:0,0,(location==3||location==4)?-.05:0)));
        arb->joint_mesh.Add_Articulation(ladder_body->id_number,bone_body->id_number,arb->joint_mesh.Add_Joint(joint));
    }
}
//#####################################################################
// Load Splines
//#####################################################################
void Load_Joint_Frames()
{
    std::string filename="Ladder/starting_joint_positions";
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false);
    if(!input) return;
    std::cout << "Loading from " << filename << std::endl;

    ARRAY<T> keyframe_times;
    int number_of_keyframes,number_of_joints;
    *input>>number_of_keyframes;
    keyframe_times.Resize(number_of_keyframes);
    for(int i=0;i<keyframe_times.m;i++) (*input)>>keyframe_times(i);;

    *input>>number_of_joints;
    starting_joint_positions.Resize(number_of_joints);
    for(int i=0;i<starting_joint_positions.m;i++){
        (*input)>>starting_joint_positions(i).x;
        FRAME_3D<T> frame;
        for(int j=1;j<=number_of_keyframes-1;j++)(*input)>>starting_joint_positions(i).y;
        (*input)>>frame;
    }
}

ROTATION<TV> Get_Joint_Rotation(const std::string& joint_name)
{
    std::cout<<"Geting joint: "<<joint_name<<std::endl;
    for(int i=0;i<starting_joint_positions.m;i++) if(starting_joint_positions(i).x==joint_name) return starting_joint_positions(i).y.r;
    std::cout<<"    joint not found!!!"<<std::endl;
    return ROTATION<TV>();
}
FRAME_TRACK_3D<T>* Make_Frame_Track(const ROTATION<TV>& start_frame,const ROTATION<TV>& end_frame,JOINT<TV>* joint,const bool right=false,const int cycle_type=1)
{
    ROTATION<TV> start=start_frame,end=end_frame;TV euler;
//    if(rotate_start){start.Euler_Angles(euler.x,euler.y,euler.z);start=ROTATION<TV>::From_Euler_Angles(euler.x,-euler.y,-euler.z);start=ROTATION<TV>::Switch_Hemisphere(end,start);}
//    else{
    end.Euler_Angles(euler.x,euler.y,euler.z);end=ROTATION<TV>::From_Euler_Angles(euler.x,-euler.y,-euler.z);end=ROTATION<TV>::Switch_Hemisphere(start,end);
    int samples=500;T period=parameter_list.Get_Parameter("cycle_period",(T)5);
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(2*samples+1,0,period);
    frame_track->periodic=true;
    frame_track->name=joint->name+"_track";
    //bool use_splines=true;
    ARRAY<ROTATION<TV> > control_points;
    ARRAY<T> control_points_times;
    if(cycle_type==-1){
        std::cout<<"alternating track found\n"<<std::endl;
        control_points.Append(start);control_points.Append(end);control_points.Append(start);control_points.Append(end); control_points.Append(start);
        control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);control_points_times.Append(4);}
    else{
        if(cycle_type==1){ // ex: hip, shoulder although these are getting their own now . . 
            control_points.Append(start);control_points.Append(right?start:end);control_points.Append(end);control_points.Append(right?start:end);
            control_points.Append(start);control_points.Append(right?start:end);control_points.Append(end);}
        else if(cycle_type==2){ //knee
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("knee_extra",(T)-.98),TV(1,0,0));
            if(right){
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                control_points_times.Append(3);control_points_times.Append(5);control_points_times.Append(6);}
            else{
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(7);}
        }
        else if(cycle_type==3){ // ex: elbow
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("elbow_extra",(T).8),TV(1,0,0));
            if(right){
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                control_points_times.Append(3);control_points_times.Append(5);control_points_times.Append(6);}
            else{
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(7);}}
        
/*
  control_points.Append(start);control_points.Append(right?extra*end:end);control_points.Append(end);control_points.Append(right?start:extra*start);
  control_points.Append(start);control_points.Append(end);control_points.Append(end);
  control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);
  control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(6);}*/
/*
  if(right){
  control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);control_points.Append(start);
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);
                    control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(6);}
                else{
                    control_points.Append(start);control_points.Append(end);control_points.Append(end);control_points.Append(extra*start);
                    control_points.Append(start);control_points.Append(end);control_points.Append(end);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);
                    control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(6);}
*/
/*
                if(right){
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                    control_points_times.Append(6);control_points_times.Append(7);control_points_times.Append(8);}
                else{
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                    control_points_times.Append(6);control_points_times.Append(9);control_points_times.Append(10);}
*/
        else if(cycle_type==4){ //hip
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("hip_extra",(T)-.15),TV(1,0,0));
            control_points.Append(start);control_points.Append(right?extra*end:end);control_points.Append(right?start:extra*start);control_points.Append(right?extra*end:end); control_points.Append(right?start:extra*start);
            control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);control_points_times.Append(4);}
/*
                if(right){
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                    control_points_times.Append(4);control_points_times.Append(5.9);control_points_times.Append(6);}
                else{
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points_times.Append(0);control_points_times.Append(2);control_points_times.Append(2.9);
                    control_points_times.Append(3);control_points_times.Append(5);control_points_times.Append(5.9);}
                    }*/
        else if(cycle_type==5){ // ex: shoulder
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("shoulder_extra",(T).3),TV(1,0,0));
            if(right){
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                control_points_times.Append(3);control_points_times.Append(4);control_points_times.Append(5);}
            else{
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(7);}}
/*
                if(right){
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points.Append(start);control_points.Append(extra*end);control_points.Append(end);
                    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                    control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(6);}
                else{
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points.Append(start);control_points.Append(end);control_points.Append(extra*start);
                    control_points_times.Append(0);control_points_times.Append(2);control_points_times.Append(3);
                    control_points_times.Append(4);control_points_times.Append(6);control_points_times.Append(7);}
                    }*/
        else if(cycle_type==6){ // ex: ankle
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("ankle_extra",(T)-.98),TV(1,0,0));
            if(right){
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*end);
                control_points.Append(start);control_points.Append(end);control_points.Append(extra*end);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                control_points_times.Append(3);control_points_times.Append(5);control_points_times.Append(6);}
            else{
                control_points.Append(start);control_points.Append(end);control_points.Append(start);
                control_points.Append(extra*start);control_points.Append(end);control_points.Append(start);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(7);}}
        else{ // ex: thumb
/*
            control_points.Append(start);control_points.Append(start);control_points.Append(end);control_points.Append(end);
            control_points.Append(start);control_points.Append(start);control_points.Append(end);
            control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);
            control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(6);
*/
            ROTATION<TV> extra=ROTATION<TV>(parameter_list.Get_Parameter("thumb_extra",(T)0),TV(1,0,0));
            ROTATION<TV> mid_right=ROTATION<TV>(0.752962,0.540454,-0.0579377,0.370947);
            ROTATION<TV> mid_left=ROTATION<TV>(0.777346,0.470812,0.0707424,-0.411174);
            if(right){
                control_points.Append(start);control_points.Append(mid_right);control_points.Append(end);
                control_points.Append(start);control_points.Append(mid_right);control_points.Append(end);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
                control_points_times.Append(3);control_points_times.Append(4);control_points_times.Append(5);}
            else{
                control_points.Append(start);control_points.Append(end);control_points.Append(mid_left);
                control_points.Append(start);control_points.Append(end);control_points.Append(mid_left);
                control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(3);
                control_points_times.Append(4);control_points_times.Append(5);control_points_times.Append(7);}}
//            control_points.Append(right?start:end);control_points.Append(right?start:end);control_points.Append(right?start:end);
//            control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);
//}}
    }
    BSPLINE_ROTATION<TV> bspline_q(control_points_times,control_points,parameter_list.Get_Parameter("spline_order",1));bspline_q.Quaternion_Check();
    //bspline_q.Create_Closed_Points();
    bspline_q.Print_Control_Points_And_Times();
    T s_inc=1/((float)2*samples);
    T range=bspline_q.Range();
    T s=bspline_q.Start_Time();
    std::cout<<"start: "<<start<<" range: "<< range<<std::endl;
    for(int i=0;i<2*samples;i++) {
        frame_track->trajectory(i+1)=FRAME_3D<T>(bspline_q.Evaluate(s+s_inc*i*range));
        //if(!right && cycle_type==5) frame_track->trajectory(i+1)=FRAME_3D<T>(bspline_q.Evaluate(s+fmod(s_inc*i*range+range/2,range)));
        if(cycle_type==5)std::cout<<"i: "<<i<<" val: "<<frame_track->trajectory(i+1)<<std::endl;}
    frame_track->trajectory(2*samples+1)=FRAME_3D<T>(bspline_q.Evaluate(s));if(cycle_type==5)std::cout<<"last: "<<frame_track->trajectory(2*samples+1)<<std::endl;
    return frame_track;
}
//#####################################################################
// Function Push_Up
//#####################################################################
void Climb()
{
    Load_Joint_Frames();    
    ROTATION<TV> elbow_angle_min=Get_Joint_Rotation("joint_r_humeroulnar");
    ROTATION<TV> elbow_angle_max=Get_Joint_Rotation("joint_r_humeroulnar_left");
    ROTATION<TV> shoulder_angle_min=Get_Joint_Rotation("joint_r_glenohumeral");
    ROTATION<TV> shoulder_angle_max=Get_Joint_Rotation("joint_r_glenohumeral_left");
    ROTATION<TV> hip_angle_min=Get_Joint_Rotation("joint_r_hip");
    ROTATION<TV> hip_angle_max=Get_Joint_Rotation("joint_r_hip_left");
    ROTATION<TV> knee_angle_min=Get_Joint_Rotation("joint_r_knee");
    ROTATION<TV> knee_angle_max=Get_Joint_Rotation("joint_r_knee_left");
    ROTATION<TV> ankle_angle_min=Get_Joint_Rotation("joint_r_ankle");
    ROTATION<TV> ankle_angle_max=Get_Joint_Rotation("joint_r_ankle_left");
    ROTATION<TV> wrist_angle_min=Get_Joint_Rotation("joint_r_radiocarpal");
    ROTATION<TV> wrist_angle_max=Get_Joint_Rotation("joint_r_radiocarpal_left");
    ROTATION<TV> thumb_angle_min=Get_Joint_Rotation("wrist_to_1mc");
    ROTATION<TV> thumb_angle_max=Get_Joint_Rotation("wrist_to_1mc_left");

    //if(!restart){
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_GLENOHUMERAL)->Set_Joint_Frame(FRAME_3D<T>(shoulder_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HUMEROULNAR)->Set_Joint_Frame(FRAME_3D<T>(elbow_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_GLENOHUMERAL)->Set_Joint_Frame(FRAME_3D<T>(shoulder_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HUMEROULNAR)->Set_Joint_Frame(FRAME_3D<T>(elbow_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HIP)->Set_Joint_Frame(FRAME_3D<T>(hip_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_KNEE)->Set_Joint_Frame(FRAME_3D<T>(knee_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HIP)->Set_Joint_Frame(FRAME_3D<T>(hip_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_KNEE)->Set_Joint_Frame(FRAME_3D<T>(knee_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME_3D<T>(ankle_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_ANKLE)->Set_Joint_Frame(FRAME_3D<T>(ankle_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_1MC)->Set_Joint_Frame(FRAME_3D<T>(thumb_angle_min));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_1MC)->Set_Joint_Frame(FRAME_3D<T>(thumb_angle_max));
        //da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_RADIOCARPAL)->Set_Joint_Frame(FRAME_3D<T>(wrist_angle_min));
        //da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_RADIOCARPAL)->Set_Joint_Frame(FRAME_3D<T>(wrist_angle_max));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_RADIOULNAR)->Set_Joint_Frame(FRAME_3D<T>(Get_Joint_Rotation("joint_r_radioulnar")));
        da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_RADIOULNAR)->Set_Joint_Frame(FRAME_3D<T>(Get_Joint_Rotation("joint_r_radioulnar_left")));
        arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_CRANIUM)->id_number);
        //}

    T k_p=parameter_list.Get_Parameter("k_p",(T)100);
    for(int i=1;i<=da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->Set_k_p(k_p);
        joint_function->Set_Target_Angle(joint_function->Angle());
    }

    if(parameter_list.Get_Parameter("climb_use_track",false)){
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HUMEROULNAR));
        frame_tracks.Append(Make_Frame_Track(elbow_angle_min,elbow_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HUMEROULNAR),true,3));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HUMEROULNAR));
        frame_tracks.Append(Make_Frame_Track(elbow_angle_max,elbow_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HUMEROULNAR),false,3));


        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_GLENOHUMERAL));
        frame_tracks.Append(Make_Frame_Track(shoulder_angle_min,shoulder_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_GLENOHUMERAL),true,5));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_GLENOHUMERAL));
        frame_tracks.Append(Make_Frame_Track(shoulder_angle_max,shoulder_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_GLENOHUMERAL),false,5));

        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HIP));
        frame_tracks.Append(Make_Frame_Track(hip_angle_min,hip_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_HIP),true,4));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HIP));
        frame_tracks.Append(Make_Frame_Track(hip_angle_max,hip_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_HIP),false,4));

        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_KNEE));
        frame_tracks.Append(Make_Frame_Track(knee_angle_min,knee_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_KNEE),true,2));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_KNEE));
        frame_tracks.Append(Make_Frame_Track(knee_angle_max,knee_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_KNEE),false,2));

        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_ANKLE));
        frame_tracks.Append(Make_Frame_Track(ankle_angle_min,ankle_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_ANKLE),true,6));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_ANKLE));
        frame_tracks.Append(Make_Frame_Track(ankle_angle_max,ankle_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_ANKLE),false,6));

        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_1MC));
        frame_tracks.Append(Make_Frame_Track(thumb_angle_min,thumb_angle_max,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_R_1MC),true,7));
        frame_track_joints.Append(da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_1MC));
        frame_tracks.Append(Make_Frame_Track(thumb_angle_max,thumb_angle_min,da_man->joint(VISIBLE_HUMAN<T,RW>::JOINT_L_1MC),false,7));

        for(int i=0;i<frame_tracks.m;i++) frame_track_joints(i)->joint_function->track=frame_tracks(i);
    }
    if(parameter_list.Get_Parameter("manual_joints",false)){
        Add_Ladder_Joint(parameter_list.Get_Parameter("joint_1",0),parameter_list.Get_Parameter("ladder_1",0));
        Add_Ladder_Joint(parameter_list.Get_Parameter("joint_2",0),parameter_list.Get_Parameter("ladder_2",0));
        Add_Ladder_Joint(parameter_list.Get_Parameter("joint_3",0),parameter_list.Get_Parameter("ladder_3",0));
        Add_Ladder_Joint(parameter_list.Get_Parameter("joint_4",0),parameter_list.Get_Parameter("ladder_4",0));}
}

//#####################################################################
// Function Reinit_Joints
//#####################################################################
void Add_Joints_To_ARB()
{
    for(int i=0;i<saved_joint.m;i++) arb->joint_mesh.Add_Articulation(saved_parent(i),saved_child(i),arb->joint_mesh.Add_Joint(saved_joint(i)));
}
void Reinit_Joints()
{
    if(parameter_list.Get_Parameter("rinit_joints",false)){
    std::cout<<"Attempting to reinit joints"<<std::endl;
    arb->Remove_All();
    
    bool force_left_hand=parameter_list.Get_Parameter("force_left_hand_check",false);
    bool force_right_hand=parameter_list.Get_Parameter("force_right_hand_check",false);
    bool force_left_foot=parameter_list.Get_Parameter("force_left_foot_check",false);
    bool force_right_foot=parameter_list.Get_Parameter("force_right_foot_check",false);

    T ladder_radius=parameter_list.Get_Parameter("ladder_scale",(T)1)*parameter_list.Get_Parameter("radius_hand",(T).02); // ladder radius normally
    bool found=false;
    if(force_left_hand || force_right_hand){
        if(force_left_hand) for(int i=0;i<ladder_bodies.m;i++) left_hand_joint=(left_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(2),ladder_radius));
        if(force_right_hand) for(int i=0;i<ladder_bodies.m;i++)  right_hand_joint=(right_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(1),ladder_radius));
    }
    else if(right_hand_joint){
        for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(2),ladder_radius));
        if(found){right_hand_joint=false;left_hand_joint=true;}
        else for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(1),ladder_radius));
    }
    else if(left_hand_joint){
        for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(1),ladder_radius));
        if(found){right_hand_joint=true;left_hand_joint=false;}
        else for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(2),ladder_radius));
    }
    else{
        for(int i=0;i<ladder_bodies.m;i++) right_hand_joint=(right_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(1),ladder_radius));
        for(int i=0;i<ladder_bodies.m;i++) left_hand_joint=(left_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(2),ladder_radius));
        std::cout<<"checking for hand joints. right hand: "<< right_hand_joint<<" left hand: "<<left_hand_joint<<std::endl;
    }
    found=false;
    T ladder_radius_foot=parameter_list.Get_Parameter("ladder_scale",(T)1)*parameter_list.Get_Parameter("radius_foot",(T).02); // ladder radius normally
    if(force_left_foot || force_right_foot){
        if(force_left_foot) for(int i=0;i<ladder_bodies.m;i++) left_foot_joint=(left_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(4),ladder_radius_foot));
        if(force_right_foot) for(int i=0;i<ladder_bodies.m;i++) right_foot_joint=(right_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(3),ladder_radius_foot));
    }
    else if(right_foot_joint){
        for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(4),ladder_radius_foot));
        if(found){right_foot_joint=false;left_foot_joint=true;}
        else for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(3),ladder_radius_foot));
    }
    else if(left_foot_joint){
        for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(3),ladder_radius_foot));
        if(found){right_foot_joint=true;left_foot_joint=false;}
        else for(int i=0;i<ladder_bodies.m;i++) found=(found || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(4),ladder_radius_foot));
    }
    else{
        for(int i=0;i<ladder_bodies.m;i++) right_foot_joint=(right_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(3),ladder_radius_foot));
        for(int i=0;i<ladder_bodies.m;i++) left_foot_joint=(left_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(4),ladder_radius_foot));
    }
    if(parameter_list.Get_Parameter("check_for_all_joints",false)) {
        for(int i=0;i<ladder_bodies.m;i++) right_hand_joint=(right_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(1),ladder_radius));
        for(int i=0;i<ladder_bodies.m;i++) left_hand_joint=(left_hand_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(2),ladder_radius));
        for(int i=0;i<ladder_bodies.m;i++) right_foot_joint=(right_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(3),ladder_radius_foot));
        for(int i=0;i<ladder_bodies.m;i++) left_foot_joint=(left_foot_joint || Check_For_New_Joint(ladder_bodies(i),dynamic_joint_check_bodies(4),ladder_radius_foot));}
    }
}
bool Check_For_New_Joint(RIGID_BODY<TV>* ladder_body,RIGID_BODY<TV>* bone_body,const T ladder_radius)
{
    TV closest_ladder_point=bone_body->frame.t;
    float distance=ladder_body->Implicit_Geometry_Extended_Value(closest_ladder_point);
    while(fabs(distance)>POINT_TOLERANCE){
        closest_ladder_point+=-ladder_body->Implicit_Geometry_Extended_Normal(closest_ladder_point)*distance;
        distance=ladder_body->Implicit_Geometry_Extended_Value(closest_ladder_point);}
    if(fabs((closest_ladder_point-bone_body->frame.t).Magnitude())<ladder_radius){
        // add dynamic joint
        JOINT<TV>* joint;
        if(parameter_list.Get_Parameter("use_only_point_joints",false)){
            // joint=new POINT_JOINT<TV>();
            JOINT<TV>* primary_point_joint=new POINT_JOINT<TV>();
            JOINT<TV>* secondary_point_joint=new POINT_JOINT<TV>();
            int primary_joint_index=arb->joint_mesh.Add_Joint(primary_point_joint);
            int secondary_joint_index=arb->joint_mesh.Add_Joint(secondary_point_joint);
            arb->joint_mesh.Add_Articulation(ladder_body->id_number,bone_body->id_number,primary_joint_index);
            arb->joint_mesh.Add_Articulation(ladder_body->id_number,bone_body->id_number,secondary_joint_index);
            primary_point_joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(bone_body->frame.Inverse()*ladder_body->frame));
            secondary_point_joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(bone_body->frame.Inverse()*(ladder_body->frame.t+TV(.05,0,0))));
            secondary_point_joint->Set_Joint_To_Parent_Frame(FRAME_3D<T>(TV(.05,0,0)));
        }
        else{
            joint=new ANGLE_JOINT<TV>();
            joint->Set_Joint_To_Child_Frame(FRAME_3D<T>(bone_body->frame.Inverse()*ladder_body->frame));
            arb->joint_mesh.Add_Articulation(ladder_body->id_number,bone_body->id_number,arb->joint_mesh.Add_Joint(joint));
        }
        std::cout<<"ADDING DYNAMIC JOINT between "<<bone_body->name<<" and ladder rung "<<std::endl;
        return true;}
    return false;
}

//#####################################################################
// Preprocess_frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    Reinit_Joints();
}
//#####################################################################
};
}
#endif
