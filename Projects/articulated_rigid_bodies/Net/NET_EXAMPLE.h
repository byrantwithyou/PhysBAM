//#####################################################################
// Copyright 2004-2007, Craig Schroeder, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NET_EXAMPLE
//#####################################################################
#ifndef __NET_EXAMPLE__
#define __NET_EXAMPLE__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <sstream>
#include "../ARB_PARAMETERS.h"
#include "../VISIBLE_HUMAN.h"
namespace PhysBAM{

template<class T>
class NET_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;typedef VECTOR<T,3> TV;
public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::stream_type;
    using BASE::output_directory;using BASE::fluids_parameters;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;

    SOLIDS_STANDARD_TESTS<TV> tests;
    ARTICULATED_RIGID_BODY<TV>* arb;
    RIGID_BODY<TV> *handle[4];
    TV pos[4];
    T increment;
    bool use_cylinders;
    PARAMETER_LIST parameter_list;
    int num_rows,num_cols;
    VISIBLE_HUMAN<T>* da_man;
    FRAME<TV> skeleton_frame;
    T scale_factor;
    bool zero_angular_velocity_of_net_joints;
    bool kinematic_corners;
    bool static_corners;
    bool move_corners;
    INTERPOLATION_CURVE<T,TV> corner_position_curves[4];
    bool net_global_post_stabilization;
    ARRAY<JOINT_ID> fused_joints;
    bool test_frame_tracks;
    int root;
    T clamp_thorax_speed;
    bool initialized_collision_manager;
    ARRAY<RIGID_BODY<TV>*> man_bones;
    bool verbose;

    NET_EXAMPLE(const STREAM_TYPE stream_type,std::string& parameter_file="")
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solids_parameters)
    {
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.cfl=(T).1;
        fluids_parameters.simulate=false;

        last_frame=1000;
        frame_rate=24;
        output_directory="Net/output";
        LOG::cout<<"Frame rate: "<<frame_rate<<std::endl;
        increment=0;

        arb=new ARTICULATED_RIGID_BODY<TV>(solid_body_collection.deformable_object.particles,solids_parameters.rigid_body_parameters.list);
        solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        //arb->Set_Extra_Iterations_Per_Contact_Level_Factor(2000);
        //arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(2000);
        //arb->Set_Poststabilization_Iterations(2000);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);

        if(parameter_file.empty()) parameter_file="Net/example.param";
        ARB_PARAMETERS::Read_Common_Parameters(parameter_file,*this,parameter_list);
        use_cylinders = parameter_list.Get_Parameter("use_cylinders",false);
        num_rows=parameter_list.Get_Parameter("num_rows",8);
        num_cols=parameter_list.Get_Parameter("num_cols",8);
        scale_factor=parameter_list.Get_Parameter("scale_factor",(T)1);
        zero_angular_velocity_of_net_joints=parameter_list.Get_Parameter("zero_angular_velocity_of_net_joints",true);
        kinematic_corners=parameter_list.Get_Parameter("kinematic_corners",true);
        static_corners=parameter_list.Get_Parameter("static_corners",false);
        move_corners=parameter_list.Get_Parameter("move_corners",true);
        net_global_post_stabilization=parameter_list.Get_Parameter("net_global_post_stabilization",false);
        verbose=parameter_list.Get_Parameter("verbose",false);

        test_frame_tracks=parameter_list.Get_Parameter("test_frame_tracks",false);
        if(test_frame_tracks) solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=false;
        clamp_thorax_speed=parameter_list.Get_Parameter("clamp_thorax_speed",(T)0);

        initialized_collision_manager=false;
        da_man=0;
    }

    ~NET_EXAMPLE()
    {delete arb;}

    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}

//#####################################################################
// Add_Fused_Hands_And_Feet
//#####################################################################
void Add_Fused_Hands_And_Feet(const bool add_left=false)
{
    RIGID_BODY<TV>* left_hand=0,*right_hand=0,*left_foot=0,*right_foot=0;
    RIGID_BODY<TV>* rigid_body=0;

    rigid_body=&tests.Add_Rigid_Body("New_Visible_Human_Bones/fused_hand_right",1,(T)1);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Name("fused_hand_right");
    rigid_body->Set_Frame(skeleton_frame*rigid_body->Frame());
    right_hand=rigid_body;
    man_bones.Append(rigid_body);

    rigid_body=&tests.Add_Rigid_Body("New_Visible_Human_Bones/fused_foot_right",1,(T)1);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Name("fused_foot_right");
    rigid_body->Set_Frame(skeleton_frame*rigid_body->Frame());
    right_foot=rigid_body;
    man_bones.Append(rigid_body);

    if(add_left){
        rigid_body=&tests.Add_Rigid_Body("New_Visible_Human_Bones/fused_hand_left",1,(T)1);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Name("fused_hand_left");
        rigid_body->Set_Frame(skeleton_frame*rigid_body->Frame());
        left_hand=rigid_body;
        man_bones.Append(rigid_body);

        rigid_body=&tests.Add_Rigid_Body("New_Visible_Human_Bones/fused_foot_left",1,(T)1);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Name("fused_foot_left");
        rigid_body->Set_Frame(skeleton_frame*rigid_body->Frame());
        left_foot=rigid_body;
        man_bones.Append(rigid_body);}

    JOINT_ID right_hand_joint,right_foot_joint;

    JOINT<TV>* joint=0;

    // HAND
    if(da_man->use_only_point_joints) joint=new POINT_JOINT<TV>();
    else joint=new ANGLE_JOINT<TV>();
    TV RC_FE_origin_in_anatomical_radius_frame(-(T).01834,-(T).00034,-(T).27811);
    TV RC_FE_axis_in_anatomical_radius_frame((T)0.966275,(T)0.0,-(T)0.257513);
    FRAME<TV> RC_FE_in_anatomical_radius_frame(RC_FE_origin_in_anatomical_radius_frame,da_man->Rotation_From_Axis(RC_FE_axis_in_anatomical_radius_frame,TV(0,1,0)));
    FRAME<TV> parent_frame=da_man->anatomical_frame_in_rigid_body_frame(VISIBLE_HUMAN<T>::BONE_R_RADIUS)*RC_FE_in_anatomical_radius_frame;
    joint->Set_Joint_To_Parent_Frame(parent_frame); //radius
    //NOTE: we don't have anatomical frame for hand -- instead we assume the initial configuration of the hand is the reference configuration
    FRAME<TV> child_frame=right_hand->Frame().Inverse()*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->Frame()*parent_frame;
    joint->Set_Joint_To_Child_Frame(child_frame); //wrist/hand
    joint->name="joint_r_radiocarpal";
    right_hand_joint=da_man->Add_Joint_With_Bend_Fix(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->particle_index,right_hand->particle_index,joint);
    fused_joints.Append(right_hand_joint);

    // ANKLE
    joint=new POINT_JOINT<TV>();
    TV ankle_joint_center((T)0.00380192,(T)0.00355539,-(T)0.259706);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(ankle_joint_center));
    joint->Set_Joint_To_Child_Frame(right_foot->Frame().Inverse()*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TIBIA)->Frame()*FRAME<TV>(ankle_joint_center));
    joint->name="joint_r_ankle";
    right_foot_joint=da_man->Add_Joint_With_Bend_Fix(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TIBIA)->particle_index,right_foot->particle_index,joint);
    fused_joints.Append(right_foot_joint);

    if(add_left){
        JOINT_ID id=da_man->Add_Additional_Reflected_Left_Joint(right_hand_joint,da_man->bones(VISIBLE_HUMAN<T>::BONE_L_RADIUS),left_hand,1);
        fused_joints.Append(id);
        id=da_man->Add_Additional_Reflected_Left_Joint(right_foot_joint,da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TIBIA),left_foot,2);
        fused_joints.Append(id);}

    if(add_left){
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/fused_hand",right_hand,left_hand,1);
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/fused_foot",right_foot,left_foot,2);}
    else{
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/fused_hand",right_hand);
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/fused_foot",right_foot);}
}
//#####################################################################
// Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    int id=0;
    T x_shift,y_shift,z_shift;
    x_shift=parameter_list.Get_Parameter("x_shift",0);
    y_shift=parameter_list.Get_Parameter("y_shift",2);
    z_shift=parameter_list.Get_Parameter("z_shift",0);
    //T scale_factor=parameter_list.Get_Parameter("scale_factor",1);
    RIGID_BODY<TV>* rigid_body=0;
    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=rigid_body_list.rigid_bodies;

    for(int row=0;row<num_rows;row++){
        for(int col=0;col<num_cols;col++){
            if(use_cylinders) rigid_body=&tests.Add_Rigid_Body("skinnycyllink_sub",scale_factor,(T).5);
            else rigid_body=&tests.Add_Rigid_Body("superskinnyhexlink",scale_factor,(T).5);
            rigid_body->X()=scale_factor*TV(x_shift+1+2*col,y_shift,z_shift+2*row);
            rigid_body->Rotation()=QUATERNION<T>(-(T)pi/2,TV(0,0,1));
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->Set_Name("net");}
        for(int col=0;col<num_cols+1;col++){
            rigid_body=&tests.Add_Rigid_Body("sphere",scale_factor,(T).5);
            rigid_body->X()=scale_factor*TV(x_shift+2*col,y_shift,z_shift+2*row);
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->Set_Name("net_joint");} // name is important (see update rigid body params)
        for(int col=0;col<num_cols+1;col++){
            if(use_cylinders) rigid_body=&tests.Add_Rigid_Body("skinnycyllink_sub",scale_factor,(T).5);
            else rigid_body=&tests.Add_Rigid_Body("superskinnyhexlink",scale_factor,(T).5);
            rigid_body->X()=scale_factor*TV(x_shift+2*col,y_shift,z_shift+1+2*row);
            rigid_body->Rotation()=QUATERNION<T>((T)pi/2,TV(1,0,0));
            rigid_body->Set_Coefficient_Of_Restitution(0);
            rigid_body->Set_Name("net");}
        if(row==num_rows-1){
            for(int col=0;col<num_cols+1;col++){
                rigid_body=&tests.Add_Rigid_Body("sphere",(T).1*scale_factor,(T).5);
                rigid_body->X()=scale_factor*TV(x_shift+2*col,y_shift,z_shift+2+2*row);
                rigid_body->Set_Coefficient_Of_Restitution(0);
                rigid_body->Set_Name("net_joint");}}} // name is important (see update rigid body params)

    for(int col=0;col<num_cols;col++){
        if(use_cylinders) rigid_body=&tests.Add_Rigid_Body("skinnycyllink_sub",scale_factor,(T).5);
        else rigid_body=&tests.Add_Rigid_Body("superskinnyhexlink",scale_factor,(T).5);
        rigid_body->X()=scale_factor*TV(x_shift+1+2*col,y_shift,z_shift+2*num_rows);
        rigid_body->Rotation()=QUATERNION<T>(-(T)pi/2,TV(0,0,1));
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Name("net");}

//    int number_of_net_rigid_bodies=rigid_body_list.rigid_bodies.m;

    // normalize their mass
    if(parameter_list.Is_Defined("uniform_mass")){
        T mass=parameter_list.Get_Parameter("uniform_mass",(T)1);
        for(int i=0;i<rigid_bodies.m;i++) rigid_bodies(i)->Set_Mass(mass);}
    else{
        for(int i=0;i<rigid_bodies.m;i++) rigid_bodies(i)->Set_Mass(1000*rigid_bodies(i)->Volume());}

    int num_joints=2*(num_cols*(num_rows+1)+num_rows*(num_cols+1));
    JOINT<TV>** joints;joints=new JOINT<TV>*[num_joints];
    int i=0;
    T link_length=rigid_bodies(1)->Object_Space_Bounding_Box().ymax;

    for(int row=0;row<num_rows+1;row++){
        for(int col=0;col<num_cols+1;col++){
            if(row==0){
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1),int(1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1),int(2*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}
               else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(2*num_cols+1),int(num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(2*num_cols+1),int(3*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(col+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(num_cols+1+col),int(2*num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}}
            else if(row==num_rows){
                int total=num_cols*num_rows+2*(num_cols+1)*num_rows;
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1),int(total-num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1),int(total+num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0))));} // --
                else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+2*num_cols+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0))));} // --
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total-num_cols+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total+num_cols+1+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+1+col),int(total+num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0))));}} // --
            else{
                int total=num_cols*row+2*(num_cols+1)*row;
                if(col==0){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total-num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1),int(total+2*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}
                else if(col==num_cols){
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total+num_cols),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+2*num_cols+1),int(total+3*num_cols+2),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}
                else{
                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+col+1),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))*QUATERNION<T>((T)pi/2,TV(0,1,0)))); // --

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total-num_cols+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));

                    joints[i]=new POINT_JOINT<TV>();arb->joint_mesh.Add_Articulation(int(total+num_cols+1+col),int(total+2*num_cols+2+col),joints[i]);
                    joints[i]->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,-link_length,0),QUATERNION<T>((T)pi/2,TV(0,1,0))));
                    joints[i++]->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(),QUATERNION<T>((T)pi/2,TV(0,1,0))));}}}}
    assert(i==num_joints);

    if(parameter_list.Get_Parameter("constrain_twist",false)) for(int i=0;i<num_joints;i++) ((POINT_JOINT<TV>*)joints[i])->Use_Phi_Constraint(0,0);
    //T net_angular_damping=parameter_list.Get_Parameter("net_angular_damping",(T)0.1);
    for(int i=0;i<num_joints;i++){
        joints[i]->global_post_stabilization=net_global_post_stabilization;}
        //joints[i]->angular_damping=net_angular_damping;}

    //handles
    handle[0]=rigid_bodies(num_cols+1);
    handle[1]=rigid_bodies(2*num_cols+1);
    handle[2]=rigid_bodies(num_cols*num_rows+2*(num_cols+1)*num_rows+1);
    handle[3]=rigid_bodies(num_cols*num_rows+2*(num_cols+1)*num_rows+num_cols+1);
    for(int i=0;i<4;i++) pos[i]=handle[i]->X();
    TV center=(T).25*(pos[0]+pos[1]+pos[2]+pos[3]),top(center.x,center.y,center.z);
    T corners_start_time=parameter_list.Get_Parameter("corners_start_time",(T)1);
    T corners_stop_time=parameter_list.Get_Parameter("corners_stop_time",(T)1);
    T corners_stop_fraction=parameter_list.Get_Parameter("corners_stop_fraction",(T).5);
    for(int i=0;i<4;i++){
        corner_position_curves[i].Add_Control_Point(corners_start_time,pos[i]);
        corner_position_curves[i].Add_Control_Point(corners_stop_time,(T)corners_stop_fraction*top+(T)(1-corners_stop_fraction)*pos[i]);}

    for(int i=0;i<4;i++){
        if(move_corners && kinematic_corners) solid_body_collection.rigid_body_collection.rigid_body_particle.kinematic(handle[i]->particle_index)=true;
        else if(static_corners) handle[i]->is_static=true;}

    // boxes
    if(parameter_list.Get_Parameter("use_boxes",false)){
        T box_scale=parameter_list.Get_Parameter("box_scale",(T)1);
        T box_mass=parameter_list.Get_Parameter("box_mass",(T)1);
        TV box_offset=parameter_list.Get_Parameter("box_offset",TV());

        TV center=(T).25*(handle[0]->X()+handle[1]->X()+handle[2]->X()+handle[3]->X())+box_offset;

        rigid_body=&tests.Add_Rigid_Body("subdivided_box",box_scale,(T).5);
        rigid_body->X()=center+(T)2.2*box_scale*TV::Axis_Vector(1);
        rigid_body->Set_Mass(box_mass);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

        rigid_body=&tests.Add_Rigid_Body("subdivided_box",box_scale,(T).5);
        rigid_body->X()=center-(T)2.2*box_scale*TV::Axis_Vector(1);
        rigid_body->Set_Mass(box_mass*5);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

        rigid_body=&tests.Add_Rigid_Body("subdivided_box",box_scale,(T).5);
        rigid_body->X()=center+(T)2.2*box_scale*TV::Axis_Vector(3);
        rigid_body->Set_Mass(box_mass*10);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

        rigid_body=&tests.Add_Rigid_Body("subdivided_box",box_scale,(T).5);
        rigid_body->X()=center-(T)2.2*box_scale*TV::Axis_Vector(3);
        rigid_body->Set_Mass(box_mass*5);
        rigid_body->Set_Coefficient_Of_Restitution((T)0.5);

        JOINT<TV>* joint=new POINT_JOINT<TV>();
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(box_scale*(T)1.5,0,0)));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-box_scale*(T)1.5,0,0)));
        arb->joint_mesh.Add_Articulation(int(id-3),int(id-2),joint);
        int first_joint_id=arb->joint_mesh.joints.m;

        joint=new POINT_JOINT<TV>();
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(box_scale*(T)1.5,0,0)));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-box_scale*(T)1.5,0,0)));
        arb->joint_mesh.Add_Articulation(int(id-2),int(id-1),joint);

        joint=new POINT_JOINT<TV>();
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(box_scale*(T)1.5,0,0)));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-box_scale*(T)1.5,0,0)));
        arb->joint_mesh.Add_Articulation(int(id-1),int(id),joint);

        arb->Update_With_Breadth_First_Directed_Graph(int(id-3));

        T k_p=parameter_list.Get_Parameter("k_p",(T)100);
        for(int i=0;i<3;i++){
            JOINT_FUNCTION<TV>* joint_function=arb->Create_Joint_Function(JOINT_ID(first_joint_id+i-1));
            joint_function->Set_k_p(k_p);
            joint_function->Set_Target_Angle(joint_function->Angle());}}

    // add the little man

    if(parameter_list.Get_Parameter("add_man",true)){
        TV man_shift=parameter_list.Get_Parameter("man_shift",TV(0,4,-(T)0.2));
        TV man_prerotate=parameter_list.Get_Parameter("man_prerotate",TV());
    //    skeleton_frame=FRAME<TV>(TV(0,(T)2.3,-(T)0.2))*FRAME<TV>(QUATERNION<T>((T)0.3,TV(1,0,0)))*FRAME<TV>(TV(0,-(T)0.12,0),QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T)0.5*(T)pi,TV(1,0,0)));
        skeleton_frame=FRAME<TV>(man_shift)*FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(man_prerotate))*FRAME<TV>(TV(0,-(T)0.12,0),QUATERNION<T>(parameter_list.Get_Parameter("man_rotate_y",(T)pi/2),TV(0,1,0))*QUATERNION<T>(parameter_list.Get_Parameter("man_rotate_x",(T)-pi/2),TV(1,0,0)));
        da_man=new VISIBLE_HUMAN<T>(stream_type,arb,rigid_body_list.rigid_body_particles,data_directory,skeleton_frame,verbose);
        da_man->use_fused_tibia_and_patella=true;
        da_man->use_only_point_joints=parameter_list.Get_Parameter("use_only_point_joints",false);
        da_man->use_double_point_joints_for_bend=parameter_list.Get_Parameter("use_double_point_joints_for_bend",false);
        if(parameter_list.Get_Parameter("bone_phi",false)) da_man->read_phi=true;

        if(parameter_list.Get_Parameter("use_fused_hands_and_feet",true)){
            da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Basic_Filter(true,true,true,false,false,true));

            if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index;
            else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index;
            else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index;

            Add_Fused_Hands_And_Feet(true);}
        else da_man->Initialize_Bodies();

        for(int i=0;i<da_man->bones.m;i++) if(da_man->bones(i)->particle_index) man_bones.Append(da_man->bones(i));

        if(parameter_list.Is_Defined("skeleton_mass")){
            T mass=parameter_list.Get_Parameter("skeleton_mass",(T)100);
            Set_Skeleton_Mass(mass);}

        VISIBLE_HUMAN<T>::Read_Joint_Limits(arb,parameter_list);

        T k_p=parameter_list.Get_Parameter("k_p",(T)100);
        for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
            assert(!da_man->joint(i)->secondary_point_of_bend_joint);
            da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
            joint_function->muscle_control=true;
            joint_function->Set_k_p(k_p);
            joint_function->Set_Target_Angle(joint_function->Angle());}
        for(int i=0;i<fused_joints.m;i++){
            assert(!arb->joint_mesh(fused_joints(i))->secondary_point_of_bend_joint);
            arb->Create_Joint_Function(fused_joints(i));JOINT_FUNCTION<TV>* joint_function=arb->joint_mesh(fused_joints(i))->joint_function;
            joint_function->muscle_control=true;
            joint_function->Set_k_p(k_p);
            joint_function->Set_Target_Angle(joint_function->Angle());}

        if(parameter_list.Is_Defined("soft_k_p")){
            T soft_k_p=parameter_list.Get_Parameter("soft_k_p",(T)100);
            std::string soft_joints="joint_r_glenohumeral joint_r_humeroulnar joint_l1 joint_l2 joint_l3 joint_l4 joint_l5 joint_l5_coccyx joint_r_hip joint_r_knee";
            parameter_list.Get_Parameter("soft_joints",soft_joints);
            std::istringstream sstr(soft_joints);std::string joint;
            LOG::cout<<"Setting soft joints: ";
            while(sstr >> joint){
                for(int i=0;i<arb->joint_mesh.joints.m;i++) if(arb->joint_mesh.joints(i)->joint_function)
                    if(arb->joint_mesh.joints(i)->name==joint || arb->joint_mesh.joints(i)->name==joint+"_left"){
                        LOG::cout<<arb->joint_mesh.joints(i)->name<<" ";
                        arb->joint_mesh.joints(i)->joint_function->Set_k_p(soft_k_p);}}
            LOG::cout<<std::endl;}

        if(parameter_list.Get_Parameter("use_motion_track",true)){
            LOG::cout<<"Reading in motion track"<<std::endl;

            ARRAY<PAIR<std::string,BSPLINE_QUATERNION<T>* > > motion_spline;
            std::string motion_track_name=parameter_list.Get_Parameter("motion_track",std::string("Net/motion_body_repeat"));
            VISIBLE_HUMAN<T>::Load_Splines(motion_track_name,motion_spline,3);

            T track_speedup_factor=parameter_list.Get_Parameter("track_speedup_factor",(T)10);
            for(int i=0;i<arb->joint_mesh.joints.m;i++) if(arb->joint_mesh.joints(i)){JOINT<TV>* joint=arb->joint_mesh.joints(i);
                if(!joint->joint_function || joint->name.find("radiocarpal")!=std::string::npos) continue;
                joint->joint_function->track=VISIBLE_HUMAN<T>::Read_Frame_Track_From_Spline(motion_spline,joint->name,true,1000,track_speedup_factor);
                if(!joint->joint_function->track) continue;
                LOG::cout<<"Got frame track for "<<joint->name<<std::endl;
                // clean up invalid rotations
                if(joint->joint_type==JOINT<TV>::TYPE_ANGLE_JOINT || joint->primary_point_of_bend_joint){
                    LOG::cout<<"Cleaning up frame track..."<<std::endl;
                    for(int i=1;i<=joint->joint_function->track->trajectory.m;i++){
                        FRAME<TV> frame=joint->joint_function->track->trajectory(i);
                        TV euler;frame.r.Euler_Angles(euler.x,euler.y,euler.z);frame.r.From_Euler_Angles(euler.x,0,0);
                        joint->joint_function->track->trajectory(i)=frame;}}}}}

    // the ground!
    if(parameter_list.Get_Parameter("add_ground",true))
        tests.Add_Ground((T).5,0,1);

//    T ether_viscosity=parameter_list.Get_Parameter("ether_viscosity",(T)0);
//    T ether_angular_viscosity=parameter_list.Get_Parameter("ether_angular_viscosity",(T)0);

    // net gets ether drag
/*    LOG::cout<<"Adding ether drag to net but not man"<<std::endl;
    for(int i=0;i<number_of_net_rigid_bodies;i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,ether_viscosity,ether_angular_viscosity);
*/
    tests.Add_Gravity();

    T coefficient_of_friction=parameter_list.Get_Parameter("coefficient_of_friction",(T)1);
    for(int i=0;i<rigid_bodies.m;i++)
        rigid_bodies(i)->Set_Coefficient_Of_Friction(coefficient_of_friction);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);

    if(test_frame_tracks || parameter_list.Get_Parameter("update_joints",false)) Update_Joints(0);
}
//#####################################################################
// Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    RIGID_BODY_COLLISION_MANAGER& collision_manager=solids_evolution->rigid_body_collisions->collision_manager;
    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=rigid_body_list.rigid_bodies;
    if(!initialized_collision_manager){initialized_collision_manager=true;
        collision_manager.Use_Collision_Matrix();

        if(parameter_list.Get_Parameter("no_collisions_between_joined_bodies",false)){
            LOG::cout<<"Turning off collisions between all joined bodies!"<<std::endl;
            VISIBLE_HUMAN<T>::Turn_Off_Collisions(arb,collision_manager);}

        if(da_man && !parameter_list.Get_Parameter("skeleton_collisions",true)){
            LOG::cout<<"turning off self collisions for skeleton"<<std::endl;
            for(int i=0;i<man_bones.m;i++) for(int j=i+1;j<=man_bones.m;j++) if(man_bones(i)->particle_index && man_bones(j)->particle_index){
                collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(man_bones(i)->particle_index,man_bones(j)->particle_index,false);
                collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(man_bones(j)->particle_index,man_bones(i)->particle_index,false);}}

        if(parameter_list.Get_Parameter("no_net_collisions",false)){
            for(int i=0;i<rigid_bodies.m;i++){
                if(rigid_bodies(i)->name=="net" || rigid_bodies(i)->name=="net_joint"){
                    for(int j=i+1;j<=rigid_bodies.m;j++){
                        if(rigid_bodies(j)->name=="net" || rigid_bodies(j)->name=="net_joint"){
                            collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(rigid_bodies(i)->particle_index,rigid_bodies(j)->particle_index,false);
                            collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(rigid_bodies(j)->particle_index,rigid_bodies(i)->particle_index,false);}}}}}
        else if(parameter_list.Get_Parameter("no_collisions",false)){
            for(int i=0;i<rigid_bodies.m;i++) for(int j=i+1;j<=rigid_bodies.m;j++){
                collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(rigid_bodies(i)->particle_index,rigid_bodies(j)->particle_index,false);
                collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(rigid_bodies(j)->particle_index,rigid_bodies(i)->particle_index,false);}}
        else{
            for(int row=0;row<num_rows+1;row++){
                for(int col=0;col<num_cols+1;col++){
                    if(row==0){
                        if(col==0){
                            int ids[2]={int(1),int(2*num_cols+2)};
                            Do_Not_Collide_All_Combinations(collision_manager,2,int(num_cols+1),ids);}
                        else if(col==num_cols){
                            int ids[2]={int(num_cols),int(3*num_cols+2)};
                            Do_Not_Collide_All_Combinations(collision_manager,2,int(2*num_cols+1),ids);}
                        else{
                            int ids[3]={int(col),int(col+1),int(2*num_cols+2+col)};
                            Do_Not_Collide_All_Combinations(collision_manager,3,int(num_cols+1+col),ids);}}
                    else if(row==num_rows){
                        int total=num_cols*num_rows+2*(num_cols+1)*num_rows;
                        if(col==0){
                            int ids[2]={int(total-num_cols),int(total+num_cols+2)};
                            Do_Not_Collide_All_Combinations(collision_manager,2,int(total+1),ids);}
                        else if(col==num_cols){
                            int ids[2]={int(total),int(total+2*num_cols+1)};
                            Do_Not_Collide_All_Combinations(collision_manager,2,int(total+num_cols+1),ids);}
                        else{
                            int ids[3]={int(total-num_cols+col),int(total+num_cols+1+col),int(total+num_cols+2+col)};
                            Do_Not_Collide_All_Combinations(collision_manager,3,int(total+1+col),ids);}}
                    else{
                        int total=num_cols*row+2*(num_cols+1)*row;
                        if(col==0){
                            int ids[3]={int(total+1),int(total-num_cols),int(total+2*num_cols+2)};
                            Do_Not_Collide_All_Combinations(collision_manager,3,int(total+num_cols+1),ids);}
                        else if(col==num_cols){
                            int ids[3]={int(total),int(total+num_cols),int(total+3*num_cols+2)};
                            Do_Not_Collide_All_Combinations(collision_manager,3,int(total+2*num_cols+1),ids);}
                        else{
                            int ids[4]={int(total+col),int(total+col+1),int(total-num_cols+col),int(total+2*num_cols+2+col)};
                            Do_Not_Collide_All_Combinations(collision_manager,4,int(total+num_cols+1+col),ids);}}}}}}

    PHYSBAM_FATAL_ERROR("this should be done in kinematic rigid body callbacks");
    if(move_corners && static_corners){
        LOG::cout<<"moving corners in update rigid body parameters...time="<<time<<"\n";
        for(int i=0;i<4;i++) handle[i]->X()=corner_position_curves[i].Value(time);}

    // Zero out some angular velocities for improved stability
    if(zero_angular_velocity_of_net_joints){
        for(int i=0;i<rigid_bodies.m;i++)
            if(rigid_bodies(i)->name=="net_joint"){
                rigid_bodies(i)->Twist().angular=rigid_bodies(i)->Angular_Momentum()=TV();}
            else if(rigid_bodies(i)->name=="net"){
                rigid_bodies(i)->Twist().angular.Project_Orthogonal_To_Unit_Direction(rigid_bodies(i)->Rotation().Rotated_Y_Axis());
                rigid_bodies(i)->Update_Angular_Momentum();}}

    if(clamp_thorax_speed && da_man && da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)){
        TV& vel=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->Twist().linear;
        vel.y=max(vel.y,-clamp_thorax_speed);}
}

void Do_Not_Collide(RIGID_BODY_COLLISION_MANAGER& collision_manager,const int id1,const int id2)
{collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(id1,id2,false);
collision_manager.Set_Rigid_Body_Collides_With_Other_Rigid_Body(id2,id1,false);}

void Do_Not_Collide_All_Combinations(RIGID_BODY_COLLISION_MANAGER& collision_manager,const int num_bodies,const int center,const int* ids)
{for(int i=0;i<num_bodies;i++) Do_Not_Collide(collision_manager,center,ids[i]);
 for(int i=0;i<num_bodies-1;i++) for(int j=i+1;j<num_bodies;j++) Do_Not_Collide(collision_manager,ids[i],ids[j]);}

//#####################################################################
// Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    RIGID_BODY<TV>* rigid_body=solids_parameters.rigid_body_parameters.list(id);
    if(kinematic_corners && move_corners && solid_body_collection.rigid_body_collection.rigid_body_particle.kinematic(id)){
        int i;for(i=0;i<4;i++) if(handle[i]->particle_index==id) break;if(i>=4) return;
        frame.t=corner_position_curves[i].Value(time);}
}
//#####################################################################
// Update_Joints
//#####################################################################
void Update_Joints(const T time)
{
    for(int i=0;i<arb->joint_mesh.joints.m;i++){
        if(arb->joint_mesh.joints(i)->joint_function){
            FRAME<TV> desired_frame=FRAME<TV>(arb->joint_mesh.joints(i)->joint_function->Target_Angle(time));
            arb->joint_mesh.joints(i)->Set_Joint_Frame(desired_frame);}
        else if(arb->joint_mesh.joints(i)->secondary_point_of_bend_joint){
            std::string name=arb->joint_mesh.joints(i)->name;std::string::size_type pos=name.find("_secondary_point");
            assert(pos!=std::string::npos);std::string reference_name=name.substr(0,pos);
            for(int j=0;j<arb->joint_mesh.joints.m;j++) if(arb->joint_mesh.joints(j)->name==reference_name && arb->joint_mesh.joints(j)->joint_function){
                if(verbose) LOG::cout<<"Found reference '"<<reference_name<<"' for '"<<name<<"'"<<std::endl;
                FRAME<TV> desired_frame=FRAME<TV>(arb->joint_mesh.joints(j)->joint_function->Target_Angle(time));
                arb->joint_mesh.joints(i)->Set_Joint_Frame(desired_frame);}}}
    if(root) arb->Update_With_Breadth_First_Directed_Graph(root);
}
//#####################################################################
// Apply_Constraints
//#####################################################################
// used in test_frame_tracks case to dump the animated body without simulating
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Apply_Constraints(dt,time);
    if(test_frame_tracks) Update_Joints(time);
}
//#####################################################################
// Set_Skeleton_Mass
//#####################################################################
void Set_Skeleton_Mass(const T new_mass)
{
    T old_mass=0;for(int i=0;i<man_bones.m;i++) old_mass+=man_bones(i)->Mass();
    T scale=old_mass?new_mass/old_mass:0;
    LOG::cout<<"Setting skeleton mass to "<<new_mass<<" (old was "<<old_mass<<")"<<std::endl;
    for(int i=0;i<man_bones.m;i++) man_bones(i)->Set_Mass(scale*man_bones(i)->Mass());
}
//#####################################################################
// Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout<<"*** IN Read_Output_Files_Solids"<<std::endl;
    BASE::Read_Output_Files_Solids(frame);
    RIGID_BODY_LIST<TV>& rigid_body_list=solids_parameters.rigid_body_parameters.list;
    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=rigid_body_list.rigid_bodies;
    if(parameter_list.Get_Parameter("zero_velocities_after_restart",false)){
        LOG::cout<<"Zeroing velocities after restart"<<std::endl;
        for(int i=0;i<rigid_bodies.m;i++){
            rigid_bodies(i)->Twist().linear=rigid_bodies(i)->Twist().angular=TV();
            rigid_bodies(i)->Update_Angular_Momentum();}}
    if(parameter_list.Is_Defined("uniform_mass")){
        T mass=parameter_list.Get_Parameter("uniform_mass",(T)1);
        LOG::cout<<"Updating net mass after restart to "<<mass<<std::endl;
        for(int i=0;i<rigid_bodies.m;i++)
            if(rigid_bodies(i)->name=="net" || rigid_bodies(i)->name=="net_joint")
                rigid_bodies(i)->Set_Mass(mass);}
    if(parameter_list.Is_Defined("skeleton_mass")){
        T mass=parameter_list.Get_Parameter("skeleton_mass",(T)100);
        Set_Skeleton_Mass(mass);}

    if(parameter_list.Get_Parameter("reset_motion_curve_after_restart",false)){
        //handles
        for(int i=0;i<4;i++) pos[i]=handle[i]->X();
        TV center=(T).25*(pos[0]+pos[1]+pos[2]+pos[3]),top(center.x,center.y,center.z);
        T corners_start_time=parameter_list.Get_Parameter("corners_start_time",(T)1);
        T corners_stop_time=parameter_list.Get_Parameter("corners_stop_time",(T)1);
        T corners_stop_fraction=parameter_list.Get_Parameter("corners_stop_fraction",(T).5);
        for(int i=0;i<4;i++){
            corner_position_curves[i]=INTERPOLATION_CURVE<T,TV>();
            corner_position_curves[i].Add_Control_Point(corners_start_time,pos[i]);
            corner_position_curves[i].Add_Control_Point(corners_stop_time,(T)corners_stop_fraction*top+(T)(1-corners_stop_fraction)*pos[i]);}}
}
};
}
#endif
