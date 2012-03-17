//#####################################################################
// Copyright 2007, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISIBLE_HUMAN_MUSCLE_EXAMPLE
//##################################################################### 
#ifndef __VISIBLE_HUMAN_MUSCLE_EXAMPLE__
#define __VISIBLE_HUMAN_MUSCLE_EXAMPLE__

#include <PhysBAM_Tools/Data_Structures/PAIR.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Math_Tools/clamp.h>
#include <PhysBAM_Tools/Matrices/MATRIX_POLICY.h>
#include <PhysBAM_Tools/Matrices/UPPER_TRIANGULAR_MATRIX_3X3.h>
#include <PhysBAM_Tools/Read_Write/FILE_UTILITIES.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/FACE_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_FORCE_CURVE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
#include "../../articulated_rigid_bodies/VISIBLE_HUMAN.h"
#include "../ARB_PARAMETERS.h"
namespace PhysBAM{

template<class T_input>
class VISIBLE_HUMAN_MUSCLE_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE; 
    using BASE::Set_External_Velocities;using BASE::Zero_Out_Enslaved_Velocity_Nodes; // silence -Woverloaded-virtual
public:
    typedef typename MUSCLE_SEGMENT<TV>::MUSCLE_SEGMENT_TYPE T_MUSCLE_SEGMENT_TYPE;
    typedef typename ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_TYPE T_MUSCLE_SEGMENT_CURVE_TYPE;
    typedef TRIPLE<T_MUSCLE_SEGMENT_TYPE,T_MUSCLE_SEGMENT_CURVE_TYPE,ARRAY<T> > T_MUSCLE_SEGMENT_DATA;

    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::output_directory;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;
    using BASE::stream_type;using BASE::restart;using BASE::solid_body_collection;using BASE::parse_args;

    SOLIDS_STANDARD_TESTS<TV> tests;

    ARTICULATED_RIGID_BODY<TV>* arb;
    MUSCLE<TV>* muscle;
    PARAMETER_LIST parameter_list;
    T peak_force;
    int root;
    bool use_muscle_control,draw_flesh;
    VISIBLE_HUMAN<T>* da_man;

    VISIBLE_HUMAN_MUSCLE_EXAMPLE(const STREAM_TYPE stream_type)
        :SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),muscle(0){
        last_frame=3000;
        frame_rate=96;
        output_directory="Visible_Human_Muscle/output";

        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-4;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;

        ARB_PARAMETERS::Read_Common_Parameters("Surface_Muscle/example.param",*this,parameter_list);
        peak_force=parameter_list.Get_Parameter("peak_force",(T)10);
        use_muscle_control=parameter_list.Get_Parameter("use_muscle_control",false);
    }

    virtual ~VISIBLE_HUMAN_MUSCLE_EXAMPLE()
    {}

    // unused callbacks
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
  
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
}
//#####################################################################
// Function Parse_Late_Options
//#####################################################################
void Parse_Late_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Late_Options();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    // Initialize deformable bodies
    Initialize_Rigid_Bodies();

    arb=&solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Do_Final_Pass(false);
    arb->Use_Epsilon_Scale(false);
    arb->Set_Use_Shock_Propagation(false);
//  arb->Use_Muscle_Actuators();
    arb->Use_PD_Actuators();

    if(draw_flesh){
/*
        DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection;
        DEFORMABLE_PARTICLES<TV>& particles=deformable_object.particles;    
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        solid_body_collection.Add_Force(Create_Quasistatic_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,true));
        solid_body_collection.Update_Fragments();}
*/}
    
    std::cout<<"done initializing example\n";

//    if(!restart && tracks_start_time==0) Add_Tracks(0);
//   if(!restart && parameter_list.Get_Parameter("update_joints",true)) Update_Joints(0);
//    if(draw_flesh) Get_Constrained_Particle_Data();
    
    std::cout<<"SOLIDS_PARAMTER "<< solids_parameters.cfl<<std::endl;
}
//#####################################################################
// Function Initialize_Rigd_Bodies
//#####################################################################
void Initialize_Rigid_Bodies()
{
    // Initialize muscle list for future muscle creation
    arb->muscle_list->muscle_force_curve.Initialize(data_directory);

    Skeleton_In_Flesh();

    tests.Add_Ground(.5,-2);

//    da_man->Initialize_Muscle_Segments();
}
//#####################################################################
// Function Skeleton_In_Flesh
//#####################################################################
void Skeleton_In_Flesh()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    FRAME<TV> skeleton_frame=FRAME<TV>(TV(0,(T).05,0),ROTATION<TV>((T)pi/2,TV(1,0,0))*ROTATION<TV>((T)pi,TV(0,1,0)));
    da_man=new VISIBLE_HUMAN<T>(stream_type,arb,rigid_body_collection,data_directory,skeleton_frame);
    da_man->Initialize_Bodies();
    da_man->Initialize_Muscle_Segments();
    if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index;

    T k_p=parameter_list.Get_Parameter("k_p",(T)10000);

    if(parameter_list.Get_Parameter("set_static",true)){
        da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->is_static=true;
    }

    Adjust_Joints_For_Skeleton_In_Flesh();

    arb->Update_With_Breadth_First_Directed_Graph(root);

    // make joint functions to try to keep this pose
    for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->muscle_control=use_muscle_control;
        joint_function->Set_k_p(k_p);
        joint_function->Set_Target_Angle(joint_function->Angle());
    }

    for(int i=0;i<da_man->muscles.m;i++) if(da_man->muscles(i)){
        T length=da_man->muscles(i)->Total_Length();
        T total_rest_length=da_man->muscles(i)->optimal_length + da_man->muscles(i)->tendon_slack_length;
        if(total_rest_length>1.5*length){
            std::cout<<"Muscle " <<da_man->muscles(i)->name<< ": (" <<da_man->muscles(i)->optimal_length<< "+" <<da_man->muscles(i)->tendon_slack_length<< ") " <<total_rest_length<< " vs " <<length<<std::endl;}}

    // Determine which joints will get PD as opposed to muscle actuation
/*
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_STERNOCLAVICULAR)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_STERNOCLAVICULAR)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_STERNOCLAVICULAR)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_STERNOCLAVICULAR)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ACROMIOCLAVICULAR)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ACROMIOCLAVICULAR)->joint_function->muscle_control=false;
    for(int i=VISIBLE_HUMAN<T>::JOINT_AXIS;i<=VISIBLE_HUMAN<T>::JOINT_L5_COCCYX;i++) if(da_man->joint(i)) da_man->joint(i)->joint_function->muscle_control=false;
    for(int i=VISIBLE_HUMAN<T>::JOINT_R_TOE_1;i<=VISIBLE_HUMAN<T>::JOINT_R_TOE_5;i++) if(da_man->joint(i)) da_man->joint(i)->joint_function->muscle_control=false;
    for(int i=VISIBLE_HUMAN<T>::JOINT_L_TOE_1;i<=VISIBLE_HUMAN<T>::JOINT_L_TOE_5;i++) if(da_man->joint(i)) da_man->joint(i)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_PATELLA)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_PATELLA)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_PATELLA)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_PATELLA)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->joint_function->muscle_control=false;
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)) da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)->joint_function->muscle_control=false;
*/
    if(arb->use_muscle_actuators){
        LOG::cout<<"Joints with only muscle control"<<std::endl;
        for(int i=0;i<arb->joint_mesh.joints.m;i++) if(arb->joint_mesh.joints(i)->joint_function && arb->joint_mesh.joints(i)->joint_function->muscle_control)
            LOG::cout<<"\t"<<arb->joint_mesh.joints(i)->name<<std::endl;}
    for(int i=0;i<rigid_body_collection.rigid_body_particle.array_collection->Size();i++) if(rigid_body_collection.Is_Active(i)) rigid_body_collection.Rigid_Body(i).Set_Coefficient_Of_Friction(1);
}
//#####################################################################
// Function Adjust_Joints_For_Skeleton_In_Flesh
//#####################################################################
void Adjust_Joints_For_Skeleton_In_Flesh()
{
    // Un-tweak right-joints.
    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX),
                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_CLAVICLE),
                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_STERNOCLAVICULAR));
    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_CLAVICLE),
                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_SCAPULA),
                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR));
    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_SCAPULA),
                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_HUMERUS),
                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL));
    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_HUMERUS),
                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ULNA),
                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR));
    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ULNA),
                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS),
                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_RADIOULNAR));
//    Untweak_Joint(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS),
//                  da_man->bones(VISIBLE_HUMAN<T>::BONE_R_WRIST),
//                  da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_RADIOCARPAL));    


    // Update all Left Joints to be symmetric
    for(int i=VISIBLE_HUMAN<T>::medial_joints+VISIBLE_HUMAN<T>::Individual_Side_Joints()+1;i<VISIBLE_HUMAN<T>::JOINT_END;i++) da_man->Update_Reflected_Left_Joint(i);
        
    // Tweak individual left joints where flesh-mesh is not semetric (still requires additional refinement).
    Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),(T)(.035*pi));
    Augment_Child_Angle_Y(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),0); 
    Augment_Child_Angle_Z(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),0);
}
//#####################################################################
// Function Untweak_Joint
//#####################################################################
void Untweak_Joint(RIGID_BODY<TV>* parent,RIGID_BODY<TV>* child,JOINT<TV>* joint)
{
    joint->Set_Joint_To_Child_Frame(child->Frame().Inverse()*parent->Frame()*joint->F_pj());
}
//#####################################################################
// Function Augment_Angle - NOTE: NOT COMUTATIVE
//#####################################################################
void Augment_Child_Angle_X(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_twist=ROTATION<TV>(angle,TV(1,0,0));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_twist));
}
void Augment_Child_Angle_Y(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_rot=ROTATION<TV>(angle,TV(0,1,0));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_rot));
}
void Augment_Child_Angle_Z(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_bend=ROTATION<TV>(angle,TV(0,0,1));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_bend));
}
void Augment_Parent_Angle_X(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_twist=ROTATION<TV>(angle,TV(1,0,0));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
void Augment_Parent_Angle_Y(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_twist=ROTATION<TV>(angle,TV(0,1,0));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
void Augment_Parent_Angle_Z(JOINT<TV>* joint,T angle)
{
    ROTATION<TV> pre_twist=ROTATION<TV>(angle,TV(0,0,1));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
//##################$##################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
//    enlarge_nodes.Resize(solid_body_collection.particles.array_collection->Size());
//    enlarge_nodes.Fill(false);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{

}
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{

}

};
}
#endif
