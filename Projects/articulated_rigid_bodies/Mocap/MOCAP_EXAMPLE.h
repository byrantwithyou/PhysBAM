//#####################################################################
// Copyright 2004, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MOCAP_EXAMPLE
//##################################################################### 
#ifndef __MOCAP_EXAMPLE__
#define __MOCAP_EXAMPLE__

#include <PhysBAM_Tools/Grids_Uniform/GRID.h>
#include <PhysBAM_Tools/Grids_Uniform_Arrays/ARRAYS_ND.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Vectors/VECTOR_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Dynamics/Motion/FRAME_TRACK_3D.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
#include "../VISIBLE_HUMAN.h"
namespace PhysBAM{

template<class T,class RW>
class MOCAP_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>
{
    typedef VECTOR<T,3> TV;
public:
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::first_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::frame_rate;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::restart_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::output_directory;
    using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::solids_parameters;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::write_last_frame;using SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::data_directory;
    
    ARTICULATED_RIGID_BODY<TV>* arb;

    MOCAP_EXAMPLE()
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>(0,FLUIDS_PARAMETERS_3D<T>::NONE)
    {

        solids_parameters.perform_self_collision=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        //solids_parameters.cfl=.1;
        solids_parameters.gravity=0;
        solids_parameters.perform_collision_body_collisions=false;

        //restart=true;
        restart_frame=30;
        last_frame=2000;
        frame_rate=24;
        output_directory="Mocap/output";
        std::cout << "Frame rate: "<<frame_rate<<std::endl;

        arb=new ARTICULATED_RIGID_BODY<TV>(this->solids_parameters.rigid_body_parameters.list);
        this->solids_parameters.rigid_body_parameters.Set_Articulated_Rigid_Body(arb);
        
        arb->Set_Iterative_Tolerance((T)1e-4);
        arb->Set_Extra_Iterations_Per_Contact_Level_Factor(10);
        arb->Set_Extra_Iterations_Per_Shock_Propagation_Level_Factor(10);
        arb->Set_Poststabilization_Iterations(10);
        arb->Set_Use_Shock_Propagation(false);
        arb->Set_Do_Final_Pass(false);
        write_last_frame=true;

        arb->Use_PD_Actuators();        
    }

    ~MOCAP_EXAMPLE()
    {
        delete arb;
    }

// Find a joint
int Find_Joint(std::string name)
{
    std::cout<<"Looking for "<<name<<": ";
    for (int i = 1; i <= arb->joint_mesh.joints.m; i++) {
        if (arb->joint_mesh.joints(i)->name == name) {
            std::cout<<"Found!\n";
            return i;
        }
    }
    std::cout<<"Not Found!\n";
    return -1;
}
JOINT<TV>* Get_Joint(std::string name)
{
    std::cout<<"Looking for "<<name<<": ";
    for (int i = 1; i <= arb->joint_mesh.joints.m; i++) {
        if (arb->joint_mesh.joints(i)->name == name) {
            std::cout<<"Found!\n";
            return arb->joint_mesh.joints(i);
        }
    }
    std::cout<<"Not Found!\n";
    return NULL;
}
void WorldRotateJoint(std::string name, QUATERNION<T> rotation)
{
    JOINT<TV> *j = arb->joint_mesh.joints(Find_Joint(name));
    FRAME<TV> jw = j->F_pj().Inverse() * j->joint_function->parent->frame;    
    j->Set_Joint_Frame(FRAME<TV>(jw.r.Inverse() * rotation * jw.r));
}
QUATERNION<T> WorldRotationInJointSpace(std::string name, QUATERNION<T> rotation)
{
    JOINT<TV> *j = arb->joint_mesh.joints(Find_Joint(name));
    FRAME<TV> jw = j->F_pj().Inverse() * j->joint_function->parent->frame;    
    return jw.r.Inverse() * rotation * jw.r;    
}
FRAME<TV> GetJointWorldFrame(std::string name)
{
    JOINT<TV> *j = arb->joint_mesh.joints(Find_Joint(name));
    FRAME<TV> jw = j->joint_function->parent->frame * j->F_pj();    
    return jw;
}
void TransformTrack(std::string mocap, std::string newskeleton)
{
    // Get both joints
    JOINT<TV> *newjoint = Get_Joint(newskeleton);
    JOINT<TV> *oldjoint = Get_Joint(mocap);

    // Get new parent and child frames
    QUATERNION<T> Rpw = newjoint->joint_function->parent->frame.r * newjoint->F_pj().r;
    QUATERNION<T> Rcw = newjoint->joint_function->child->frame.r * newjoint->F_cj().r;    

    // Get joint world frames
    QUATERNION<T> newjw = GetJointWorldFrame(newskeleton).r, oldjw = GetJointWorldFrame(mocap).r;   

    // Build transform
    QUATERNION<T> pretransform(Rpw.Inverse() * oldjw);        
    QUATERNION<T> posttransform(oldjw.Inverse() * Rcw);              

    // Get mocap track and create new matching track
    FRAME_TRACK_3D<T> *mocaptrack = oldjoint->joint_function->track;
    FRAME_TRACK_3D<T> *newtrack = new FRAME_TRACK_3D<T>(mocaptrack->time_grid.m, mocaptrack->time_grid.xmin, mocaptrack->time_grid.xmax);
    
    // Transform data
    for (int i = 1; i <= mocaptrack->trajectory.m; i++){
        newtrack->trajectory(i) = FRAME<TV>(pretransform) * mocaptrack->trajectory(i) * FRAME<TV>(posttransform);        
    }
    
    // Set new track
    newjoint->joint_function->track = newtrack;    
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    int id=0;
    RIGID_BODY<TV> *rigid_body=0;

    // Grab list ref
    RIGID_BODY_LIST_3D<T>& rigid_body_list=solids_parameters.rigid_body_parameters.list;

    // Load the test and set static for now
    Load_Motion_Test();
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) rigid_body_particles.Rigid_Body(i).is_static = true;
    int offset = arb->joint_mesh.joints.m;    

    // Load visible human
    VISIBLE_HUMAN<T,RW>* da_man=new VISIBLE_HUMAN<T,RW>(arb,data_directory,FRAME<TV>(TV(),QUATERNION<T>(pi,TV(0,1,0))*QUATERNION<T>(-0.5*pi,TV(1,0,0))));
    da_man->Initialize_Bodies();    

    // Now walk hierarchy by using topological sort and set joint frames            
    arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_CRANIUM)->id_number);
    //for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) rigid_body_particles.Rigid_Body(i).is_static = true;

    /*id=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"ground");
    rigid_body=arb->rigid_bodies_list.rigid_bodies(id);
    rigid_body->frame.t=TV(0,-2,0);
    rigid_body->velocity=TV(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(1);
    rigid_body->Set_Coefficient_Of_Friction(.5);
    rigid_body->Set_Name("ground");
    rigid_body->is_static=true;
    rigid_body->add_to_spatial_partition=false;*/

    // Create joint functions now
    for (int i = 1; i <= da_man->joint.m; i++) da_man->Create_Joint_Function(i);    

    // Now establish correspondences between joints
    // So we want to rotate here such that we align the two components
    
     WorldRotateJoint("joint_r_ankle", 
        QUATERNION<T>((T)0.38, TV(0,0,1)) * QUATERNION<T>((T)0.43, TV(1,0,0)));
    
    WorldRotateJoint("joint_r_ankle_left", 
        QUATERNION<T>((T)-0.38, TV(0,0,1)) * QUATERNION<T>((T)0.43, TV(1,0,0)));

    WorldRotateJoint("joint_r_hip", 
        QUATERNION<T>(0.3, TV(0,0,1)));
    
    WorldRotateJoint("joint_r_hip_left", 
        QUATERNION<T>(-0.3, TV(0,0,1)));
    
    Get_Joint("joint_r_glenohumeral")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(pi/2, TV(0,1,0))));    
    Get_Joint("joint_r_glenohumeral_left")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-pi/2, TV(0,1,0))));    

    Get_Joint("joint_r_radiocarpal")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(pi/2, TV(0,0,1))));    
    Get_Joint("joint_r_radiocarpal_left")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(pi/2, TV(0,0,1))));    
    
    Get_Joint("joint_r_radioulnar")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(pi/2, TV(1,0,0))));    
    Get_Joint("joint_r_radioulnar_left")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(pi/2, TV(1,0,0))));    

    // Note, aligning the spine causes difficulties - here we get good results by keeping it fixed.
    // we'll need more precise 
    /*Get_Joint("joint_c3")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.010, TV(1,0,0))));    
    Get_Joint("joint_c4")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.010, TV(1,0,0))));    
    Get_Joint("joint_c5")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.10, TV(1,0,0))));    
    Get_Joint("joint_c6")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.10, TV(1,0,0))));    
    Get_Joint("joint_c7")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.030, TV(1,0,0))));   
    Get_Joint("joint_c7_thorax")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(0.05, TV(1,0,0))));   

    Get_Joint("joint_l1")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(0.38, TV(1,0,0))));  
    Get_Joint("joint_l2")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.04, TV(1,0,0))));        
    Get_Joint("joint_l3")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.04, TV(1,0,0))));    
    Get_Joint("joint_l4")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-0.04, TV(1,0,0))));    
    Get_Joint("joint_l5")->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(0.02, TV(1,0,0))));        */
     

    // Update configuration so both characters are correctly positioned in the same pose
    // in world space - the next functions use this fact to determine track transforms    
    arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_CRANIUM)->id_number);

    // Set lower extremity joints
    TransformTrack("rfemur", "joint_r_hip");
    TransformTrack("lfemur", "joint_r_hip_left");

    // Set lower extremity joints
    TransformTrack("rtibia", "joint_r_knee");
    TransformTrack("ltibia", "joint_r_knee_left");    

    // Ankle joints
    TransformTrack("rfoot", "joint_r_ankle");
    TransformTrack("lfoot", "joint_r_ankle_left");

    // Right toes (move together for now)
    TransformTrack("rtoes", "joint_r_toe_1");
    TransformTrack("rtoes", "joint_r_toe_2");
    TransformTrack("rtoes", "joint_r_toe_3");
    TransformTrack("rtoes", "joint_r_toe_4");
    TransformTrack("rtoes", "joint_r_toe_5");

    // Left toes (move together for now)
    TransformTrack("ltoes", "joint_r_toe_1_left");
    TransformTrack("ltoes", "joint_r_toe_2_left");
    TransformTrack("ltoes", "joint_r_toe_3_left");
    TransformTrack("ltoes", "joint_r_toe_4_left");
    TransformTrack("ltoes", "joint_r_toe_5_left");

    // Shoulder joints
    TransformTrack("rhumerus", "joint_r_glenohumeral");
    TransformTrack("lhumerus", "joint_r_glenohumeral_left");

    // Arm joints
    TransformTrack("rradius", "joint_r_humeroulnar");
    TransformTrack("lradius", "joint_r_humeroulnar_left");

    // Arm joints
    TransformTrack("rwrist", "joint_r_radiocarpal");
    TransformTrack("lwrist", "joint_r_radiocarpal_left");

    // Spine joints
    //TransformTrack("lowerback", "joint_l5");
    //TransformTrack("upperback", "joint_l2");

    // Head joints
    //TransformTrack("upperneck", "joint_c7_thorax");
    //TransformTrack("head", "joint_axis");
       

    // Now - set both characters into same first pose
    for (int i = 1; i <= arb->joint_mesh.joints.m; i++) {
        JOINT<TV> *j = arb->joint_mesh.joints(i);
        //j->Set_Joint_Frame(FRAME<TV>());
        if (j->joint_function && j->joint_function->track) {
            //j->Set_Joint_Frame(j->joint_function->track->Frame(0));
        }
    }

    // Now - update both to initial config
    arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T,RW>::BONE_CRANIUM)->id_number);

    std::cout<<"First body is: "<<rigid_body_list.rigid_bodies(1)->name<<"\n";
    arb->Update_With_Breadth_First_Directed_Graph(1);    
    
    for(int i=1;i<=rigid_body_list.Number_Of_Elements();i++) if(!rigid_body_particles.Rigid_Body(i).is_static)
        rigid_body_particles.Rigid_Body(i).Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);

    solids_parameters.collision_body_list.Add_Bodies(solids_parameters.rigid_body_parameters.list);
    arb->use_pd_actuators = true;

    SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV>,RW>::Initialize_Bodies();
    std::cout<<"done initializing example\n";
}
//#####################################################################
// Scene creation helpers
//#####################################################################
int CreateRigidBody(std::string name, FRAME<TV> frame, T scale=1)
{
    int pid=solids_parameters.rigid_body_parameters.list.template Add_Rigid_Body<RW>(data_directory+"/Rigid_Bodies/"+"subdivided_box",scale);
    RIGID_BODY<TV> *parent_body=arb->rigid_bodies_list.rigid_bodies(pid);
    parent_body->frame = frame;
    parent_body->Set_Coefficient_Of_Restitution(0.05);
    parent_body->Set_Coefficient_Of_Friction(0.7);
    parent_body->Set_Name(name); 
    return pid;
}
int CreateRigidJoint(int b1, int b2, FRAME<TV> jf)
{
    JOINT<TV>* joint=new RIGID_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
    RIGID_BODY<TV> *parent=arb->rigid_bodies_list.rigid_bodies(b1);
    RIGID_BODY<TV> *child=arb->rigid_bodies_list.rigid_bodies(b2);
    
    FRAME<TV> joint_to_parent = parent->frame.Inverse() * jf;
    FRAME<TV> joint_to_child = child->frame.Inverse() * jf;
    joint->Set_Joint_To_Parent_Frame(joint_to_parent);
    joint->Set_Joint_To_Child_Frame(joint_to_child);

    int jid = arb->joint_mesh.joints.m;
    arb->joint_mesh.Add_Articulation(b1,b2,jid);
    return jid;
}
int CreateBendJointXAxis(int b1, int b2, FRAME<TV> jf, T angle)
{
    JOINT<TV>* joint=new ANGLE_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
    RIGID_BODY<TV> *parent=arb->rigid_bodies_list.rigid_bodies(b1);
    RIGID_BODY<TV> *child=arb->rigid_bodies_list.rigid_bodies(b2);

    FRAME<TV> joint_to_parent = parent->frame.Inverse() * jf;
    FRAME<TV> joint_to_child = child->frame.Inverse() * jf;    
    joint->Set_Joint_To_Parent_Frame(joint_to_parent);
    joint->Set_Joint_To_Child_Frame(joint_to_child);

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,parent,child);
    joint->Set_Joint_Function(jfunc);
    joint->joint_function->Set_k_p(200);
    joint->joint_function->Set_Target_Angle(QUATERNION<T>(angle,TV(1,0,0)));

    int jid = arb->joint_mesh.joints.m;
    arb->joint_mesh.Add_Articulation(b1,b2,jid);
    return jid;
}
int CreatePointJoint(int b1, int b2, FRAME<TV> jf, QUATERNION<T> da)
{
    JOINT<TV>* joint=new POINT_JOINT<TV>();arb->joint_mesh.Add_Joint(joint);
    RIGID_BODY<TV> *parent=arb->rigid_bodies_list.rigid_bodies(b1);
    RIGID_BODY<TV> *child=arb->rigid_bodies_list.rigid_bodies(b2);

    FRAME<TV> joint_to_parent = parent->frame.Inverse() * jf;
    FRAME<TV> joint_to_child = child->frame.Inverse() * jf;
    joint->Set_Joint_To_Parent_Frame(joint_to_parent);
    joint->Set_Joint_To_Child_Frame(joint_to_child);

    JOINT_FUNCTION<TV>* jfunc=new JOINT_FUNCTION<TV>(joint,parent,child);
    joint->Set_Joint_Function(jfunc);
    joint->joint_function->Set_k_p(200);
    joint->joint_function->Set_Target_Angle(da);

    int jid = arb->joint_mesh.joints.m;
    arb->joint_mesh.Add_Articulation(b1,b2,jid);
    return jid;
}
//#####################################################################
// Test for loading a motion from parsed asf/amc files
//####################################################################
T read_T(std::ifstream &s)
{
    T t;
    s >> t;
    return t;
}
TV read_vector_3d(std::ifstream &s)
{
    TV vec; 
    s >> vec;    
    return vec;
}
FRAME_TRACK_3D<T>* read_rotation_track_3d(std::ifstream &s)
{
    int count = 0;s >> count;
    std::cout << "Reading motion array with " << count << " samples.\n";

    if (count==0) return NULL;
    
    FRAME_TRACK_3D<T>* track = new FRAME_TRACK_3D<T>(count, (T)0.0, (T)(count/120.0));
    for (int i = 1; i <= count; i++) {

        // Read keys
        T time = read_T(s);
        TV euler; s >> euler;
        QUATERNION<T> direct = QUATERNION<T>::From_Euler_Angles(euler.x, euler.y, euler.z);

        // Store as frame in the track
        //track->time_grid(i) = time;
        track->trajectory(i) = FRAME<TV>(direct);
    }
    return track;
}
void Load_Motion_Test()
{    
    // Global transform for mocap skeleton   
    T scale = (T)(0.45 * 0.13);
    FRAME<TV> transform(
        TV(2,1,0),
        QUATERNION<T>(pi,TV(0,1,0)));    

    // Load the mocap skeleton
    std::ifstream is("MocapData/13.asf.conv");    

    // Read in bodies
    int nbodies; is >> nbodies;
    for (int i = 0; i < nbodies; i++)
    {
        std::string name; is >> name;
        std::cout<<"Parsing body "<<i<<"\n";

        TV center_of_mass = read_vector_3d(is);
        TV euler_world_space = read_vector_3d(is);

        T length; is >> length;
        TV direction = read_vector_3d(is);

        // Override size for certain data
        T size = (T)0.2 * scale;
        /*if (name == "lfoot" || name == "rfoot") size = (T)0.6;
        else if (name == "ltoes" || name == "rtoes") size = (T)0.3;
        else if (name == "lhipjoint" || name == "rhipjoint") size = (T)0.30;
        else if (name == "ltibia" || name == "rtibia") size = (T)0.52;
        else if (name == "lfemur" || name == "rfemur") size = (T)0.58;
        else if (name == "lowerback" || name == "upperback") size = (T)0.48;
        else if (name == "thorax") size = (T)0.45;
        else if (name == "lowerneck") size = (T)0.42;
        else if (name == "upperneck") size = (T)0.2;
        else if (name == "head") size = (T)0.85;
        else if (name == "root") size = (T)0.30;
        else if (name == "lclavicle" || name == "rclavicle") size = (T)0.35;
        else if (name == "lhumerus" || name == "rhumerus") size = (T)0.39;
        else if (name == "lradius" || name == "rradius") size = (T)0.36;
        else if (name == "lwrist" || name == "rwrist") size = (T)0.30;*/
        
        // Get global rotation transform to bone setup        
        QUATERNION<T> c = QUATERNION<T>::From_Euler_Angles(euler_world_space.x, euler_world_space.y, euler_world_space.z);
        FRAME<TV> frame(center_of_mass * scale, c); frame = transform * frame;        
        
        // Orient body in this way        
        int bid = CreateRigidBody(name, frame, size);
        RIGID_BODY<TV> *body=arb->rigid_bodies_list.rigid_bodies(bid);        
    }

    // Read in articulation points
    int njoints; is >> njoints;    
    for (int i = 0; i < njoints; i++)
    {
        int jid = 0;is >> jid;
        std::string name;is >> name;
        std::cout<<"Loaded joint "<<jid<<", name is "<<name<<"\n";
        TV joint_pos = read_vector_3d(is);
        TV axis = read_vector_3d(is);            
        int parent, child;is >> parent;is >> child;

        // Read in degrees of freedom
        int dof; is >> dof;
        for (int j = 0; j < dof; j++) { std::string degree; is >> degree; }

        // Build bone axis again
        QUATERNION<T> c = QUATERNION<T>::From_Euler_Angles(axis.x, axis.y, axis.z);
        FRAME<TV> frame(joint_pos * scale, c); frame = transform * frame; 

        // Align joints with child(bone) (does this matter?)
        // TODO: We have no case for 2d joints yet..            
        JOINT<TV> *j = NULL; int id = 0;
        if (dof == 0)      { id = CreateRigidJoint(parent, child, frame); } 
        else if (dof == 1) { id = CreateBendJointXAxis(parent, child, frame, (T).0); } 
        else               { id = CreatePointJoint(parent, child, frame, QUATERNION<T>(0, TV(1,0,0))); }

        // Store the name of the joint as corresponding to the bone name
        j = arb->joint_mesh.joints(id);
        j->name = name;
    }

    // Close the skeleton
    is.close();
  

    // Now load motion sequence and use to set angles?
    std::ifstream motions("MocapData/13_40.amc.conv");    

    // Read motion data
    int njointdata = 0; motions >> njointdata;
    for (int i = 0; i < njointdata; i++) 
    {
        // Get id of this joint
        int id; motions >> id;        
        JOINT<TV> *j = arb->joint_mesh.joints(id);        

        // Read in rx,ry,rz tracks
        FRAME_TRACK_3D<T> *track = read_rotation_track_3d(motions);
        if (track == NULL) continue;
                
        // Set inital joint angles                
        j->Set_Joint_Frame( FRAME<TV>() );  // Set to zero so we can match native pose
        if (j->joint_function) { j->joint_function->track = track; }

        // Set name
        track->name = j->name;
    }

    motions.close();

    // Now walk hierarchy by using topological sort and set joint frames    
    arb->Update_With_Breadth_First_Directed_Graph(1);
    
}
//#####################################################################
};
}
#endif
