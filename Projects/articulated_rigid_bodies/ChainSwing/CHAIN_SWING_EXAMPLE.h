//#####################################################################
// Copyright 2004-2007, Tamar Shinar, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class CHAIN_SWING_EXAMPLE
//##################################################################### 
#ifndef __CHAIN_SWING_EXAMPLE__
#define __CHAIN_SWING_EXAMPLE__

#include <PhysBAM_Tools/Interpolation/BSPLINE_QUATERNION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
#include "../VISIBLE_HUMAN.h"
namespace PhysBAM{

template<class T>
class CHAIN_SWING_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
    typedef VECTOR<T,3> TV;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::restart;using BASE::restart_frame;using BASE::output_directory;using BASE::solids_parameters;
    using BASE::write_last_frame;using BASE::data_directory;using BASE::Write_Output_Files;using BASE::fluids_parameters;using BASE::stream_type;
    
    ARTICULATED_RIGID_BODY<TV>* arb;
    SOLIDS_STANDARD_TESTS<TV> tests;
    bool add_ground;
    int length;
    PARAMETER_LIST parameter_list;
    int mace_joints;
    ARRAY<JOINT<TV>*> joints;
    ARRAY<int> tri_transform_list;
    FRAME<TV> skeleton_frame;
    int root;
    bool test_frame_tracks;
    bool static_hand;
    bool only_sim_mace;
    T track_speedup_factor;
    ARRAY<PAIR<std::string,BSPLINE_QUATERNION<T>* > > motion_splines;
    std::string parameter_file;
    RIGID_BODY<TV>* hand;
    RIGID_BODY<TV>* mace;
    T count;
    FRAME<TV> original_hand_frame;
    VISIBLE_HUMAN<T>* da_man;
    
    CHAIN_SWING_EXAMPLE(const STREAM_TYPE stream_type,const std::string& parameter_file="example.param")
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        //solids_parameters.gravity=0;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

        //solids_parameters.cfl=(T).1;
        //restart=true;
        restart_frame=30;
        last_frame=2000;
        frame_rate=60;

        output_directory="ChainSwing/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(solid_body_collection.deformable_object.particles,solids_parameters.rigid_body_parameters.list);
        solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb); 
        arb->Set_Iterative_Tolerance((T)1e-5);
        arb->Set_Contact_Level_Iterations(100);
        arb->Set_Shock_Propagation_Level_Iterations(100);
        arb->Set_Poststabilization_Iterations(100);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);
        arb->Use_PD_Actuators();

        length=6;
        mace_joints=0;
        ARB_PARAMETERS::Read_Common_Parameters(parameter_file,*this,parameter_list);

        Init_Tri_Transform_List();
        skeleton_frame=FRAME<TV>(TV(),QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T)0.5*(T)pi,TV(1,0,0)));

        test_frame_tracks=parameter_list.Get_Parameter("test_frame_tracks",false);
        if(test_frame_tracks) solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;

        track_speedup_factor=20;

        static_hand=parameter_list.Get_Parameter("static_hand",false);
        only_sim_mace=parameter_list.Get_Parameter("only_sim_mace",false);
        add_ground=parameter_list.Get_Parameter("add_ground",false);
        frame_rate=parameter_list.Get_Parameter("frame_rate",60);

        count=0;
    }

    ~CHAIN_SWING_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Function Read_Frame_Track
//#####################################################################
FRAME_TRACK_3D<T>* Read_Frame_Track(const std::string& filename,const bool periodic=true)
{
    int samples;T tmin,tmax;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false);
    *input>>samples>>tmin>>tmax;
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,tmin,tmin+(tmax-tmin)/track_speedup_factor);
    for(int i=0;i<samples;i++) *input>>frame_track->trajectory(i);
    delete input;
    frame_track->periodic=periodic;
    return frame_track;
}
//#####################################################################
// Function Read_Frame_Track
//#####################################################################
FRAME_TRACK_3D<T>* Read_Frame_Track_From_Function(const std::string& jointname,const bool periodic=true)
{
    int samples=1000;
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,0,parameter_list.Get_Parameter("track_speedup_factor",10));
    frame_track->periodic=periodic;
    frame_track->name=jointname+"_track";
    T increment=1/((float)samples-1);
    if(jointname=="joint_r_glenohumeral"){
        for(int i=0;i<samples;i++){
//                frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles((T)pi/2,(T)pi/4*cos(2*(T)pi*i*increment),0));}
            frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(3*(T)pi/8,0,(T)pi/8)*QUATERNION<T>::From_Euler_Angles(0,0,(T)pi/4*cos(2*(T)pi*i*increment)));}
//                frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles((T)pi/2,0,3*(T)pi/8)*QUATERNION<T>::From_Euler_Angles(0,0,(T)pi/4*cos(2*(T)pi*i*increment)));}
    }
    else if(jointname=="joint_r_humeroulnar"){
        for(int i=0;i<samples;i++){
//                frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(sin(2*(T)pi*i*increment),0,0));}
//                frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles((T)pi/2+((T)pi/6)*sin(2*(T)pi*i*increment),0,0));}
            frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles((T)pi/2+((T)pi/4)*sin(2*(T)pi*i*increment),0,0));}
    }
    else if(jointname=="joint_r_radioulnar"){
        for(int i=0;i<samples;i++){
            frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles((T)pi/2,0,0));}
//                frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(cos(2*(T)pi*i*increment),0,0));}
    }
    else if(jointname=="joint_r_wrist"){
        for(int i=0;i<samples;i++){
            frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(-1,0,1));}
    }
    else std::cout<<"oops, can't find that joint\n";
    return frame_track;
}
//#####################################################################
// Function Read_Frame_Track_From_Spline
//#####################################################################
FRAME_TRACK_3D<T>* Read_Frame_Track_From_Spline(const std::string& jointname,const bool periodic=true,const int samples=1000)
{
    T tmin=0;int index=0;
    for(int i=0;i<motion_splines.m;i++) if(motion_splines(i).x==jointname) index=i;
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,tmin,parameter_list.Get_Parameter("track_speedup_factor",10));
    T time_increment=((T)1./(T)(samples-1))*motion_splines(index).y->Range();
    for(int i=0;i<samples-1;i++) frame_track->trajectory(i+1)=FRAME<TV>(motion_splines(index).y->Evaluate(motion_splines(index).y->Start_Time()+i*time_increment));
    frame_track->trajectory(samples)=FRAME<TV>(motion_splines(index).y->Evaluate(motion_splines(index).y->Start_Time()));
    frame_track->periodic=periodic;
    return frame_track;
}
//#####################################################################
// Function Create_Elbow_Frame_Track
//#####################################################################
FRAME_TRACK_3D<T>* Create_Elbow_Frame_Track()
{
    QUATERNION<T> start(parameter_list.Get_Parameter("elbow_start",(T)pi/2),TV(1,0,0)),end(parameter_list.Get_Parameter("elbow_end",(T).2),TV(1,0,0));
    int samples=1000;T period=parameter_list.Get_Parameter("track_speedup_factor",(T)5);
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples+1,0,period);
    frame_track->periodic=true;
    ARRAY<QUATERNION<T> > control_points;
    ARRAY<T> control_points_times;
    control_points.Append(start);control_points.Append(end);control_points.Append(start);control_points.Append(end); control_points.Append(start);
    control_points_times.Append(0);control_points_times.Append(1);control_points_times.Append(2);control_points_times.Append(3);
    control_points_times.Append(4);
    BSPLINE_QUATERNION<T> bspline_q(control_points_times,control_points,parameter_list.Get_Parameter("elbow_spline_order",1));bspline_q.Quaternion_Check();
    bspline_q.Print_Control_Points_And_Times();
    T s_inc=1/((float)samples),s=bspline_q.Start_Time(),range=bspline_q.Range();std::cout<<"start "<<s<<" range "<<range<<std::endl;
    for(int i=0;i<samples;i++) frame_track->trajectory(i+1)=FRAME<TV>(bspline_q.Evaluate(s+s_inc*i*range));
    frame_track->trajectory(samples+1)=FRAME<TV>(bspline_q.Evaluate(s));
    return frame_track;
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
//    Print_Splines();
    Load_Splines();

    //Make_Mace(length,(T).03,TV(-(T).35,(T)1.2,-(T)0.1),QUATERNION<T>(-(T)pi*(T).5,TV(1,0,0)));
    Make_Mace(length,(T).03,TV(-(T).35,(T)1.165,-(T)0.13),QUATERNION<T>::Rotation_Quaternion(TV(0,0,1),TV(-(T)0.04,(T).04,-(T).04))); // params to fit in hand
    //QUATERNION<T>(-(T)pi*(T).5-(T).2,TV(1,0,0)));
    if(parameter_list.Get_Parameter("use_man",true)){
        da_man=new VISIBLE_HUMAN<T>(stream_type,arb,solid_body_collection.deformable_object.rigid_body_particles,data_directory,skeleton_frame,true);
        da_man->use_only_point_joints=parameter_list.Get_Parameter("use_point_joints_only",false);
        if(parameter_list.Get_Parameter("use_upper_body",true)){
            da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Arm_And_Shoulder_Filter(false,true,true));
            root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index;
            solids_parameters.rigid_body_parameters.list(root)->is_static=true;
            Add_Wrist_Joint(da_man);}
        else{
            da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::All_But_Right_Hand_Filter());
            root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index;

            Add_Wrist_Joint(da_man);
            // put feet flat
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.6,-(T)0.123,(T)0.0477))));
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.689,-(T)0.00996,-(T)0.0822813))));
            arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index);
            da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ANKLE)->is_static=true;
            da_man->bones(VISIBLE_HUMAN<T>::BONE_L_ANKLE)->is_static=true;}

        T k_p=parameter_list.Get_Parameter("k_p",(T)100);
        for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i) && !da_man->joint(i)->joint_function){
            JOINT_FUNCTION<TV>* joint_function=da_man->Create_Joint_Function(i);
            joint_function->Set_k_p(k_p);joint_function->Set_Target_Angle(joint_function->Angle());}
        
        // to get muscles to attach to replaced hand
        da_man->Replace_Bones_With_Fused_Bone("mace_hand_transform",solids_parameters.rigid_body_parameters.list(int(1)));}
    else root=int(1);
//    Articulate_Skeleton_Hand();
    //Articulate_Mace()

    if(add_ground) tests.Add_Ground((T).5,-2,1);

    tests.Add_Gravity();
    solid_body_collection.Update_Fragments();

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Initialize_Bodies();

    std::cout<<"done initializing example\n";

    Update_Joints(0);

    // make chain and mace drop down vertically
    if(mace_joints){
        QUATERNION<T> desired_mace_and_chain_world_orientation(QUATERNION<T>((T)pi,TV(1,0,0)));
        assert(joints(1)==arb->joint_mesh.joints(1) && arb->joint_mesh.Joint_Index_From_Id(joints(1)->id_number)==1);
        joints(1)->Set_Joint_Frame(FRAME<TV>((arb->Parent(joints(1)->id_number)->Frame()*joints(1)->F_pj()).r.Inverse()*desired_mace_and_chain_world_orientation));
        Update_Joints(0);}

    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    if(static_hand) for(int i=0;i<rigid_body_list.rigid_bodies.m;i++)
        if(rigid_body_list.rigid_bodies(i)->name!="mace" && rigid_body_list.rigid_bodies(i)->name!="cyllink") rigid_body_list.rigid_bodies(i)->is_static=true;
}
//#####################################################################
// Function Init_Tri_Transform_List
//#####################################################################
void Init_Tri_Transform_List()
{
    tri_transform_list.Append(1);
    for(int i=15;i<=31;i++) tri_transform_list.Append(i);
}
//#####################################################################
// Function Make_Mace
//#####################################################################
void Make_Mace(const int length,const T& scale_factor=(T)1,const TV & shift,const QUATERNION<T>& orient=QUATERNION<T>())
{
    int num_bodies=solids_parameters.rigid_body_parameters.list.rigid_bodies.m;
    RIGID_BODY<TV>* rigid_body=0;
/*
    rigid_body=&tests.Add_Rigid_Body(data_directory+"/Rigid_Bodies/handle_thin",(T).5*scale_factor,1);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("handle_thin");
    rigid_body->frame=FRAME<TV>(shift,orient)*FRAME<TV>(scale_factor*TV(0,0,-(T)5.75),QUATERNION<T>((T)pi*(T).5,TV(1,0,0)));
*/  


    rigid_body=&tests.Add_Rigid_Body(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/holding",1,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("holding");
    if(!parameter_list.Get_Parameter("use_man",true)) rigid_body->is_static=true;
    rigid_body->Set_Frame(skeleton_frame*rigid_body->Frame());
    hand=rigid_body;
    original_hand_frame=rigid_body->Frame();

    if(parameter_list.Get_Parameter("no_mace",false)) return;
    
    for(int i=0;i<length;i++){
        rigid_body=&tests.Add_Rigid_Body(data_directory+"/Rigid_Bodies/cyllink",(T).25*scale_factor,0);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
        rigid_body->Set_Frame(FRAME<TV>(shift,orient)*FRAME<TV>(scale_factor*TV(0,-(T).5+i,0)));
        rigid_body->Set_Name("cyllink");
        rigid_body->Set_Mass(1);}

    rigid_body=&tests.Add_Rigid_Body(data_directory+"/Rigid_Bodies/mace",(T).4*scale_factor,0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("mace");
//    rigid_body->Set_Mass(parameter_list.Get_Parameter("mace_mass",(T)1)*rigid_body->mass); 
    rigid_body->Set_Mass(parameter_list.Get_Parameter("mace_mass",(T)1));
    rigid_body->Set_Frame(FRAME<TV>(shift,orient)*FRAME<TV>(scale_factor*TV(0,1+length,0)));
    mace=rigid_body;

    mace_joints=length+1;
    joints.Resize(length+2);
    for(int i=0;i<length+1;i++) joints(i)=new POINT_JOINT<TV>();
    //handle to chain
    joints(1)->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(-(T)0.003,-(T)0.213,(T)0.006)));
    joints(1)->Set_Joint_To_Child_Frame(FRAME<TV>(scale_factor*TV(0,-(T).5,0)));
    arb->joint_mesh.Add_Articulation(int(num_bodies+1),int(num_bodies+2),joints(1));
    //chains
    for(int i=2;i<=length;i++){
        joints(i)->Set_Joint_To_Parent_Frame(FRAME<TV>(scale_factor*TV(0,(T).5,0)));
        joints(i)->Set_Joint_To_Child_Frame(FRAME<TV>(scale_factor*TV(0,-(T).5,0)));
        arb->joint_mesh.Add_Articulation(int(num_bodies+i),int(num_bodies+1+i),joints(i));}
    //chain to mace
    joints(length+1)->Set_Joint_To_Parent_Frame(FRAME<TV>(scale_factor*TV(0,(T).5,0)));
    joints(length+1)->Set_Joint_To_Child_Frame(FRAME<TV>(scale_factor*TV(0,-1,0)));
    arb->joint_mesh.Add_Articulation(int(num_bodies+length+1),int(num_bodies+length+2),joints(length+1));

/*
    rigid_body=&tests.Add_Rigid_Body(data_directory+"/Rigid_Bodies/sphere",scale_factor,0);
    rigid_body=solids_parameters.rigid_body_parameters.list.rigid_bodies(id);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->Set_Name("body1");
    rigid_body->Set_Frame(FRAME<TV>(shift,orient)*FRAME<TV>(scale_factor*TV(0,0,-12)));
    rigid_body->Set_Mass(rigid_body->mass*(T).1);
    rigid_body->is_static=true;

    joints.Append(new POINT_JOINT<TV>());
    joints(joints.m)->Set_Joint_To_Parent_Frame(FRAME<TV>());
    joints(joints.m)->Set_Joint_To_Child_Frame(FRAME<TV>(scale_factor*TV(0,-(T)6.25,0)));
    arb->joint_mesh.Add_Articulation(int(num_bodies+length+3),int(num_bodies+1),joints(joints.m));
*/
}
//#####################################################################
// Function Add_Wrist_Joint
//#####################################################################
void Add_Wrist_Joint(VISIBLE_HUMAN<T>* man)
{
    JOINT_FUNCTION<TV>* jfunc;
    ARRAY<JOINT<TV>*>& joints_list=arb->joint_mesh.joints;

    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

#if 1
    JOINT<TV>* joint=new POINT_JOINT<TV>();
    TV RC_FE_origin_in_anatomical_radius_frame(-(T).01834,-(T).00034,-(T).27811);
    TV RC_FE_axis_in_anatomical_radius_frame((T)0.966275,(T)0.0,-(T)0.257513);

    TV RC_FE_axis_in_radius_bone_frame=man->anatomical_frame_in_rigid_body_frame(VISIBLE_HUMAN<T>::BONE_R_RADIUS)*RC_FE_axis_in_anatomical_radius_frame;
    QUATERNION<T> two_dof_wrist_joint=man->Rotation_From_Axis(man->joint(VISIBLE_HUMAN<T>::JOINT_R_RADIOULNAR)->F_cj()*TV(1,0,0),RC_FE_axis_in_radius_bone_frame);
    FRAME<TV> parent_frame(man->anatomical_frame_in_rigid_body_frame(VISIBLE_HUMAN<T>::BONE_R_RADIUS)*RC_FE_origin_in_anatomical_radius_frame,two_dof_wrist_joint);
    ((POINT_JOINT<TV>*)joint)->Use_Twist_Constraint(0,0);

/*    FRAME<TV> RC_FE_in_anatomical_radius_frame(RC_FE_origin_in_anatomical_radius_frame,man->Rotation_From_Axis(RC_FE_axis_in_anatomical_radius_frame,TV(0,1,0)));
    FRAME<TV> parent_frame=man->anatomical_frame_in_rigid_body_frame(VISIBLE_HUMAN<T,RW>::BONE_R_RADIUS)*RC_FE_in_anatomical_radius_frame;
*/
    joint->Set_Joint_To_Parent_Frame(parent_frame); //radius
   //NOTE: we don't have anatomical frame for hand -- instead we assume the initial configuration of the hand is the reference configuration
    FRAME<TV> child_frame=solids_parameters.rigid_body_parameters.list.rigid_bodies(1)->Frame().Inverse()*man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->Frame()*parent_frame;
    joint->Set_Joint_To_Child_Frame(child_frame); //wrist/hand
    arb->joint_mesh.Add_Articulation(man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->particle_index,int(1),joint);


    bool use_function=parameter_list.Get_Parameter("use_function",false);
#if 1
    jfunc=arb->Create_Joint_Function(joints_list(joints_list.m)->id_number);
    //jfunc->Set_Target_Angle(QUATERNION<T>(0,TV(1,0,0)));
    jfunc->Set_Target_Angle(jfunc->Angle());
    joints_list(joints_list.m)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(0,TV(1,0,0))));
    joints_list(joints_list.m)->joint_function->Set_k_p(k_p);
    //jfunc->track=Read_Frame_Track("ChainSwing/hand.track");
    if(only_sim_mace || parameter_list.Get_Parameter("static_wrist",false)) jfunc->Set_Target_Angle(jfunc->Angle());
    else if(use_function) jfunc->track=Read_Frame_Track_From_Function("joint_r_wrist");
    else jfunc->track=Read_Frame_Track_From_Spline("joint_r_wrist");
#else

    int samples=1000;T period=5;
    FRAME_TRACK_3D<T>* wrist_track=new FRAME_TRACK_3D<T>(samples,0,period);shoulder_track->periodic=true;
    for(int i=0;i<samples;i++){
        wrist_track->trajectory(i)=FRAME<TV>(QUATERNION<T>((T)pi/4+((T)pi/4)*(T)0.5*(1-cos(2*(T)pi*(i-1)/(samples-1))),TV(1,0,0)));}
#endif

#if 0
    std::istream* sample_set=FILE_UTILITIES::Safe_Open_Input("../joint_ik/motion_track1",false);
    int samples=700;
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,(T)0,(T)1);
    for(int i=0;i<samples;i++) {
        (*sample_set)>>frame_track->trajectory(i);
        frame_track->trajectory(i).r.Invert();
    }
    delete sample_set;
    joint->joint_function->track=frame_track;
#endif
#endif
    jfunc=man->Create_Joint_Function(VISIBLE_HUMAN<T>::JOINT_R_STERNOCLAVICULAR);
    jfunc->Set_Target_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);

    jfunc=man->Create_Joint_Function(VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR);
    jfunc->Set_Target_Angle(jfunc->Angle());
    jfunc->Set_k_p(k_p);

    jfunc=man->Create_Joint_Function(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL);
    if(only_sim_mace) jfunc->Set_Target_Angle(jfunc->Angle());
    else if(use_function) jfunc->track=Read_Frame_Track_From_Function("joint_r_glenohumeral");
    else jfunc->track=Read_Frame_Track_From_Spline("joint_r_glenohumeral");
    jfunc->Set_k_p(k_p);

    jfunc=man->Create_Joint_Function(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR);
    if(only_sim_mace) jfunc->Set_Target_Angle(jfunc->Angle());
    else if(use_function) jfunc->track=Read_Frame_Track_From_Function("joint_r_humeroulnar");
    else //jfunc->track=Create_Elbow_Frame_Track();
        jfunc->track=Read_Frame_Track_From_Spline("joint_r_humeroulnar");
    jfunc->Set_k_p(k_p);

    jfunc=man->Create_Joint_Function(VISIBLE_HUMAN<T>::JOINT_R_RADIOULNAR);
    if(only_sim_mace) jfunc->Set_Target_Angle(jfunc->Angle());
    else if(use_function) jfunc->track=Read_Frame_Track_From_Function("joint_r_radioulnar");
    else jfunc->track=Read_Frame_Track_From_Spline("joint_r_radioulnar");
    jfunc->Set_k_p(k_p);
}
//#####################################################################
// Function Articulate_Skeleton_Hand
//#####################################################################
void Articulate_Skeleton_Hand()
{
    T finger_angle=(T)pi/2+(T).2,thumb_angle=(T).2,thumb_mc_angle=-(T).2;
    JOINT_FUNCTION<TV>* jfunc=0;
    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    ARRAY<JOINT<TV>*>& joint_mesh=arb->joint_mesh.joints;
    for(int i=mace_joints+7;i<=mace_joints+20;i++){JOINT<TV>* joint=joint_mesh(i);
        jfunc=arb->Create_Joint_Function(joint->id_number);
        jfunc->Set_Target_Angle(QUATERNION<T>(finger_angle,TV(1,0,0)));
        joint->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(finger_angle,TV(1,0,0))));
        joint->joint_function->Set_k_p(100);}
    for(int i=mace_joints+15;i<=mace_joints+16;i++){JOINT<TV>* joint=joint_mesh(i);
        jfunc=arb->Create_Joint_Function(joint->id_number);
        jfunc->Set_Target_Angle(QUATERNION<T>(thumb_angle,TV(1,0,0)));
        joint->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(thumb_angle,TV(1,0,0))));
        joint->joint_function->Set_k_p(100);}
    for(int i=mace_joints+21;i<=mace_joints+21;i++){JOINT<TV>* joint=joint_mesh(i);
        jfunc=arb->Create_Joint_Function(joint->id_number);
        jfunc->Set_Target_Angle(QUATERNION<T>(thumb_mc_angle,TV(1,0,0)));
        joint->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(thumb_mc_angle,TV(1,0,0))));
        joint->joint_function->Set_k_p(100);}
}
//#####################################################################
// Function Articulate_Mace
//#####################################################################
void Articulate_Mace()
{
    T angle=-(T)pi/2;
    JOINT_FUNCTION<TV>* jfunc;
    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    ARRAY<JOINT<TV>*>& joint_mesh=arb->joint_mesh.joints;
    int last_joint=arb->joint_mesh.joints.m;

    jfunc=new JOINT_FUNCTION<TV>(joint_mesh(last_joint),rigid_body_list(length+3),rigid_body_list(1));
    joint_mesh(last_joint)->Set_Joint_Function(jfunc);
    jfunc->Set_Target_Angle(QUATERNION<T>(angle,TV(0,1,0)));
    joint_mesh(last_joint)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(angle,TV(0,1,0))));
    joint_mesh(last_joint)->joint_function->Set_k_p(100);

    int samples=parameter_list.Get_Parameter("chain_swing_samples",(int)1000);
    T angle_magnitude=parameter_list.Get_Parameter("chain_swing_magnitude",(T).5);
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,0,parameter_list.Get_Parameter("chain_swing_period",(T)1));
    frame_track->periodic=true;
    for(int i=0;i<samples;i++) 
        frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>(-(T)pi/2+(T).5*sin(2*(T)pi*(i-1)/(samples-1)),TV(0,1,0)));
    //frame_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Euler_Angles(-(T)pi/2+angle_magnitude*sin(2*(T)pi*(i-1)/(samples-1)),0,angle_magnitude*cos(2*(T)pi*(i-1)/(samples-1))));
    joint_mesh(last_joint)->joint_function->track=frame_track;
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    //T mace_magnitude=(mace)?mace->Twist().angular.Magnitude():0;
    //if(mace_magnitude>parameter_list.Get_Parameter("max_mace_angular_velocity",(T))
    solids_evolution->rigid_body_collisions->collision_manager.Use_Collision_Matrix();
/*    for(int i=0;i<mace_joints;i++){
        solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(i,i+1,false);
        solids_evolution->rigid_body_collisions->collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(i+1,i,false);}
*/
    PHYSBAM_FATAL_ERROR("this should be done in kinematic rigid body callbacks");
    da_man->Turn_Off_Collisions(arb,collision_manager);
    if(hand && !parameter_list.Get_Parameter("use_man",true)){
        //rotate it statically to see if it works
        hand->Rotation()=QUATERNION<T>(time*parameter_list.Get_Parameter("rotation_increment",(T).1),TV(0,1,0))*original_hand_frame.r;
        hand->Twist().angular=TV(0,parameter_list.Get_Parameter("rotation_increment",(T).1),0);
        hand->Update_Angular_Momentum();
        
        //           hand->X()=TV(time*parameter_list.Get_Parameter("rotation_increment",(T).1),0,0)+original_hand_frame.t;
        //    hand->velocity=TV(parameter_list.Get_Parameter("rotation_increment",(T).1),0,0);
    }
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
/*void Write_Output_Files(const int frame) const
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Write_Output_Files(frame);

    const RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;

    std::string outputfile="tri_transform"+frame;

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(outputfile,false);
    for(int i=0;i<tri_transform_list.m;i++){
        (*output)<<data_directory<<"/Rigid_Bodies/";
        if(tri_transform_list(i)>mace_joints+1)(*output)<<"New_Visible_Human_Bones/";
        (*output)<<rigid_body_list(tri_transform_list(i))->name<<"\n"<<rigid_body_list(tri_transform_list(i))->Frame()<<std::endl;
    }
}*/
//#####################################################################
// Function Update_Joints
//#####################################################################
void Update_Joints(const T time)
{
    for(int i=0;i<arb->joint_mesh.joints.m;i++)
        if(arb->joint_mesh.joints(i)->joint_function)
            arb->joint_mesh.joints(i)->Set_Joint_Frame(FRAME<TV>(arb->joint_mesh.joints(i)->joint_function->Target_Angle(time)));
    arb->Update_With_Breadth_First_Directed_Graph(root);
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
// used in test_frame_tracks case to dump the animated body without simulating
void Apply_Constraints(const T dt,const T time)
{
    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >::Apply_Constraints(dt,time);
    if(test_frame_tracks) Update_Joints(time);
}
//#####################################################################
// Function Load_Splines
//#####################################################################
void Load_Splines()
{
    std::string filename="ChainSwing/motion_track.keyframes";
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false);
    if(!input) return;
    std::cout << "Loading from " << filename << std::endl;

    ARRAY<T> keyframe_times;
    int number_of_keyframes,number_of_joints;
    *input>>number_of_keyframes;
    keyframe_times.Resize(number_of_keyframes);
    for(int i=0;i<keyframe_times.m;i++) (*input)>>keyframe_times(i);;

    *input>>number_of_joints;
    motion_splines.Resize(number_of_joints);
    for(int i=0;i<motion_splines.m;i++){
        (*input)>>motion_splines(i).x;
        ARRAY<QUATERNION<T> > control_points;control_points.Resize(number_of_keyframes);
        for(int j=0;j<control_points.m;j++){FRAME<TV> frame;(*input)>>frame;control_points(j)=frame.r;}
        motion_splines(i).y=new BSPLINE_QUATERNION<T>(keyframe_times,control_points,parameter_list.Get_Parameter("spline_order",3));
        if(parameter_list.Get_Parameter("perform_quaternion_check",true)) motion_splines(i).y->Quaternion_Check();
        motion_splines(i).y->Create_Closed_Points();}
}
//#####################################################################
// Function Print_Splines
//#####################################################################
void Print_Splines()
{
    for(int i=0;i<motion_splines.m;i++){
        std::cout<<"Joint: "<<motion_splines(i).x<<std::endl;
        motion_splines(i).y->Print_Control_Points_And_Times();}
}
//#####################################################################
// Function Preprocess_frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
}
//#####################################################################
};
}
#endif
