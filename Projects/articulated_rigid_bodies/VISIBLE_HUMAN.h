//#####################################################################
// Copyright 2004-2007, Kevin Der, Eran Guendelman, Ranjitha Kumar, Mike Rodgers, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISIBLE_HUMAN
//#####################################################################
#ifndef __VISIBLE_HUMAN__
#define __VISIBLE_HUMAN__

#include <PhysBAM_Tools/Interpolation/BSPLINE_QUATERNION.h>
#include <PhysBAM_Tools/Matrices/MATRIX_4X4.h>
#include <PhysBAM_Tools/Parsing/PARAMETER_LIST.h>
#include <PhysBAM_Tools/Read_Write/Matrices_And_Vectors/READ_WRITE_FRAME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ANALYTIC_SURFACE_MUSCLE_SEGMENT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/ATTACHMENT_POINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Muscles/MUSCLE_SEGMENT.h>
namespace PhysBAM{

template<class T>
class VISIBLE_HUMAN
{
    typedef VECTOR<T,3> TV;
    typedef typename MUSCLE_SEGMENT<TV>::MUSCLE_SEGMENT_TYPE T_MUSCLE_SEGMENT_TYPE;
    typedef typename ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_TYPE T_MUSCLE_SEGMENT_CURVE_TYPE;
    typedef TRIPLE<T_MUSCLE_SEGMENT_TYPE,T_MUSCLE_SEGMENT_CURVE_TYPE,ARRAY<T> > T_MUSCLE_SEGMENT_DATA;
public:
    const STREAM_TYPE stream_type;
    ARTICULATED_RIGID_BODY<TV>* arb;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;
    ARRAY<RIGID_BODY<TV>*> bones;
    ARRAY<MUSCLE<TV>*> muscles;
    ARRAY<std::string> bone_names;
    ARRAY<FRAME<TV> > anatomical_frame_in_rigid_body_frame;
    ARRAY<JOINT<TV>*> joint;
    ARRAY<JOINT_ID> local_joint_index_to_arb_joint_id;
    ARRAY<PAIR<int,int> > joint_parent_and_child;
    HASHTABLE<std::string,ARRAY<T_MUSCLE_SEGMENT_DATA>*> segment_data_map;

    std::string data_directory,upper_body_file_jnt,upper_body_file_msl,lower_body_file_jnt,lower_body_file_msl,simm_directory;

    enum JOINTS{JOINT_ATLAS=1,JOINT_AXIS,JOINT_C3,JOINT_C4,JOINT_C5,JOINT_C6,JOINT_C7,JOINT_C7_THORAX,JOINT_L1,JOINT_L2,JOINT_L3,JOINT_L4,JOINT_L5,JOINT_L5_COCCYX,
               JOINT_R_STERNOCLAVICULAR,JOINT_R_ACROMIOCLAVICULAR,JOINT_R_GLENOHUMERAL,JOINT_R_HUMEROULNAR,JOINT_R_RADIOULNAR,JOINT_R_RADIOCARPAL,
               JOINT_R_HIP,JOINT_R_KNEE,JOINT_R_PATELLA,JOINT_R_ANKLE,JOINT_R_TOE_1,JOINT_R_TOE_2,JOINT_R_TOE_3,JOINT_R_TOE_4,JOINT_R_TOE_5,
               JOINT_R_5PROXPH,JOINT_R_4PROXPH,JOINT_R_3PROXPH,JOINT_R_2PROXPH,JOINT_R_5MIDPH,JOINT_R_4MIDPH,JOINT_R_3MIDPH,JOINT_R_2MIDPH,JOINT_R_1DISTPH,
               JOINT_R_1PROXPH,JOINT_R_5DISTPH,JOINT_R_4DISTPH,JOINT_R_3DISTPH,JOINT_R_2DISTPH,JOINT_R_1MC,JOINT_R_PALM,
               JOINT_L_STERNOCLAVICULAR,JOINT_L_ACROMIOCLAVICULAR,JOINT_L_GLENOHUMERAL,JOINT_L_HUMEROULNAR,JOINT_L_RADIOULNAR,JOINT_L_RADIOCARPAL,
               JOINT_L_HIP,JOINT_L_KNEE,JOINT_L_PATELLA,JOINT_L_ANKLE,JOINT_L_TOE_1,JOINT_L_TOE_2,JOINT_L_TOE_3,JOINT_L_TOE_4,JOINT_L_TOE_5,
               JOINT_L_5PROXPH,JOINT_L_4PROXPH,JOINT_L_3PROXPH,JOINT_L_2PROXPH,JOINT_L_5MIDPH,JOINT_L_4MIDPH,JOINT_L_3MIDPH,JOINT_L_2MIDPH,JOINT_L_1DISTPH,
               JOINT_L_1PROXPH,JOINT_L_5DISTPH,JOINT_L_4DISTPH,JOINT_L_3DISTPH,JOINT_L_2DISTPH,JOINT_L_1MC,JOINT_L_PALM,JOINT_END};
    enum BONE{BONE_THORAX=1,BONE_HIP,BONE_CRANIUM,BONE_ATLAS,BONE_AXIS,BONE_C3,BONE_C4,BONE_C5,BONE_C6,BONE_C7,BONE_L1,BONE_L2,BONE_L3,BONE_L4,BONE_L5,
              BONE_R_CLAVICLE,BONE_R_SCAPULA,BONE_R_HUMERUS,BONE_R_ULNA,BONE_R_RADIUS,
              BONE_R_1DISTPH,BONE_R_2DISTPH,BONE_R_3DISTPH,BONE_R_4DISTPH,BONE_R_5DISTPH,BONE_R_1MC,BONE_R_2MIDPH,BONE_R_3MIDPH,BONE_R_4MIDPH,BONE_R_5MIDPH,
              BONE_R_1PROXPH,BONE_R_2PROXPH,BONE_R_3PROXPH,BONE_R_4PROXPH,BONE_R_5PROXPH,BONE_R_PALM,BONE_R_WRIST,
              BONE_R_PATELLA,BONE_R_FEMUR,BONE_R_TIBIA,BONE_R_ANKLE,BONE_R_TOE_1,BONE_R_TOE_2,BONE_R_TOE_3,BONE_R_TOE_4,BONE_R_TOE_5,
              BONE_L_CLAVICLE,BONE_L_SCAPULA,BONE_L_HUMERUS,BONE_L_ULNA,BONE_L_RADIUS,
              BONE_L_1DISTPH,BONE_L_2DISTPH,BONE_L_3DISTPH,BONE_L_4DISTPH,BONE_L_5DISTPH,BONE_L_1MC,BONE_L_2MIDPH,BONE_L_3MIDPH,BONE_L_4MIDPH,BONE_L_5MIDPH,
              BONE_L_1PROXPH,BONE_L_2PROXPH,BONE_L_3PROXPH,BONE_L_4PROXPH,BONE_L_5PROXPH,BONE_L_PALM,BONE_L_WRIST,
              BONE_L_PATELLA,BONE_L_FEMUR,BONE_L_TIBIA,BONE_L_ANKLE,BONE_L_TOE_1,BONE_L_TOE_2,BONE_L_TOE_3,BONE_L_TOE_4,BONE_L_TOE_5,BONE_END};
    static const int medial_joints=14;
    static const int medial_bones=15;
    static const int num_bones=77;

    MATRIX<T,4> thorax_left_side_reflection,pelvis_left_side_reflection; // world space reflections to get left side bones from right
                                                                           // thorax reflection is used for upper body
                                                                           // pelvis reflection is used for lower body
    FRAME<TV> transform;
    bool verbose;
    bool use_fused_tibia_and_patella;
    bool read_phi;
    bool use_only_point_joints;
    bool use_double_point_joints_for_bend;
    T secondary_point_joint_offset;

    VISIBLE_HUMAN(const STREAM_TYPE stream_type,ARTICULATED_RIGID_BODY<TV>* arb,RIGID_BODY_COLLECTION<TV>& rigid_body_collection_input,const std::string& directory,
        const FRAME<TV> transform=FRAME<TV>(),const bool verbose=false)
        :stream_type(stream_type),arb(arb),rigid_body_collection(rigid_body_collection_input),data_directory(directory),transform(transform),verbose(verbose)
    {
        // from VH_Bones/ribs_reflect_transform.txt
        pelvis_left_side_reflection=MATRIX<T,4>((T)-0.999,(T)0.03998,(T)0.01999,(T)0,(T)0.03998,(T)0.9992,(T)-0.0003998,(T)0,(T)0.01999,(T)-0.0003998,(T)0.9998,(T)0,(T)0.575378,
            (T)-0.0115075,(T)-0.00575377,(T)1);
        thorax_left_side_reflection=MATRIX<T,4>((T)-0.996921,(T)-0.0501875,(T)-0.0602302,(T)0,(T)-0.0499594,(T)0.998737,(T)-0.00528605,(T)0,(T)-0.0604195,(T)0.00226067,(T)0.998169,
            (T)0,(T)0.68125,(T)0.0110093,(T)0.0206681,(T)1);

        upper_body_file_jnt=data_directory+"/SIMM_Data/UpperExtremityModel/Stanford VA upper limb model.jnt";
        upper_body_file_msl=data_directory+"/SIMM_Data/UpperExtremityModel/Stanford VA upper limb model.msl";
        lower_body_file_jnt=data_directory+"/SIMM_Data/LowerExtremityModel/delpjnt.jnt";
        lower_body_file_msl=data_directory+"/SIMM_Data/LowerExtremityModel/delpmusc.msl";
        simm_directory=data_directory+"/SIMM_Data/";
        use_fused_tibia_and_patella=false;
        read_phi=false;
        use_only_point_joints=false;
        use_double_point_joints_for_bend=false;
        secondary_point_joint_offset=(T).1;
    }

//#####################################################################
// Function Filter_Union
//#####################################################################
static ARRAY<bool> Filter_Union(const ARRAY<bool>& f1,const ARRAY<bool>& f2)
{
    ARRAY<bool> filter(num_bones);for(int i=0;i<num_bones;i++) filter(i)=f1(i) || f2(i);return filter;
}
//#####################################################################
// Function Reflected_Bones_Filter
//#####################################################################
static ARRAY<bool> Reflected_Bones_Filter(const ARRAY<bool>& filter)
{
    ARRAY<bool> reflected_filter(num_bones);for(int i=0;i<filter.m;i++) if(filter(i)) reflected_filter(Reflected_Bone(i))=true;return reflected_filter;
}
//#####################################################################
// Function Add_Reflected_Bones_To_Filter
//#####################################################################
static void Add_Reflected_Bones_To_Filter(ARRAY<bool>& filter)
{
    for(int i=0;i<filter.m;i++) if(filter(i)) filter(Reflected_Bone(i))=true;
}
//#####################################################################
// Function Right_Hand_Filter
//#####################################################################
static ARRAY<bool> Right_Hand_Filter()
{
    ARRAY<bool> filter(num_bones);
    filter(BONE_R_1DISTPH)=filter(BONE_R_2DISTPH)=filter(BONE_R_3DISTPH)=filter(BONE_R_4DISTPH)=filter(BONE_R_5DISTPH)=true;
    filter(BONE_R_1MC)=filter(BONE_R_2MIDPH)=filter(BONE_R_3MIDPH)=filter(BONE_R_4MIDPH)=filter(BONE_R_5MIDPH)=true;
    filter(BONE_R_1PROXPH)=filter(BONE_R_2PROXPH)=filter(BONE_R_3PROXPH)=filter(BONE_R_4PROXPH)=filter(BONE_R_5PROXPH)=true;
    filter(BONE_R_PALM)=filter(BONE_R_WRIST)=true;
    return filter;
}
//#####################################################################
// Function Arm_And_Shoulder_Filter
//#####################################################################
static ARRAY<bool> Arm_And_Shoulder_Filter(const bool with_left_side=true,const bool with_hand=true,const bool no_right_hand=false)
{
    ARRAY<bool> filter(num_bones);
    filter(BONE_THORAX)=filter(BONE_R_CLAVICLE)=filter(BONE_R_SCAPULA)=filter(BONE_R_HUMERUS)=filter(BONE_R_ULNA)=filter(BONE_R_RADIUS)=true;
    if(with_left_side) Add_Reflected_Bones_To_Filter(filter);
    if(with_hand){
        if(!no_right_hand) filter=Filter_Union(filter,Right_Hand_Filter());
        if(with_left_side) filter=Filter_Union(filter,Reflected_Bones_Filter(Right_Hand_Filter()));}
    return filter;
}
//#####################################################################
// Function All_But_Right_Hand_Filter
//#####################################################################
static ARRAY<bool> All_But_Right_Hand_Filter()
{
    ARRAY<bool> filter=Basic_Filter(true,true,true,false,true,true);filter=Filter_Union(filter,Reflected_Bones_Filter(Right_Hand_Filter()));return filter;
}
//#####################################################################
// Function No_Hands_Or_Toes_Filter
//#####################################################################
static ARRAY<bool> No_Hands_Or_Toes_Filter(const bool lower_body=true)
{
    ARRAY<bool> filter=Basic_Filter(true,true,true,false,false,true);filter(BONE_R_ANKLE)=true;filter(BONE_L_ANKLE)=true;
    filter(BONE_R_WRIST)=true;filter(BONE_L_WRIST)=true;filter(BONE_R_PALM)=true;filter(BONE_L_PALM)=true;
    if(!lower_body){filter(BONE_R_ANKLE)=filter(BONE_L_ANKLE)=filter(BONE_R_FEMUR)=filter(BONE_L_FEMUR)=filter(BONE_R_PATELLA)=filter(BONE_L_PATELLA)=false;
    filter(BONE_R_TIBIA)=filter(BONE_L_TIBIA)=filter(BONE_HIP)=false;
    filter(BONE_L1)=filter(BONE_L2)=filter(BONE_L3)=filter(BONE_L4)=filter(BONE_L5)=false;}
    return filter;
}
//#####################################################################
// Function Humerus_And_Ulna_Filter
//#####################################################################
static ARRAY<bool> Humerus_And_Ulna_Filter(const bool with_radius=true,const bool with_left_side=true)
{
    ARRAY<bool> filter(num_bones);
    filter(BONE_R_HUMERUS)=filter(BONE_R_ULNA)=true;filter(BONE_R_RADIUS)=with_radius;
    //filter(BONE_R_ULNA)=filter(BONE_R_RADIUS)=true;
    if(with_left_side) Add_Reflected_Bones_To_Filter(filter);
    return filter;
}
//#####################################################################
// Function Upper_Body_Filter
//#####################################################################
static ARRAY<bool> Upper_Body_Filter()
{
    ARRAY<bool> filter(num_bones);
    filter(BONE_R_HUMERUS)=filter(BONE_R_ULNA)=filter(BONE_R_RADIUS)=filter(BONE_R_CLAVICLE)=filter(BONE_R_SCAPULA)=true;
    for(int i=BONE_R_1DISTPH;i<=BONE_R_WRIST;i++) filter(i)=true;
    Add_Reflected_Bones_To_Filter(filter);
    filter(BONE_THORAX)=true;
    return filter;
}
//#####################################################################
// Function Full_Body_Filter
//#####################################################################
static ARRAY<bool> Full_Body_Filter()
{
    ARRAY<bool> filter(num_bones);
    for(int i=0;i<num_bones;i++) filter(i)=true;
    return filter;
}
//#####################################################################
// Function Spine_Filter
//#####################################################################
static ARRAY<bool> Spine_Filter()
{
    ARRAY<bool> filter(num_bones);
    for(int i=BONE_THORAX;i<=BONE_L5;i++) filter(i)=true;
    return filter;
}
//#####################################################################
// Function Skeleton_In_Flesh_Filter
//#####################################################################
static ARRAY<bool> Skeleton_In_Flesh_Filter()
{
    //return Full_Body_Filter();
    //return Spine_Filter();
    return Upper_Body_Filter();
    //return Filter_Union(Spine_Filter(),Upper_Body_Filter());
}
//#####################################################################
// Function Basic_Filter
//#####################################################################
static ARRAY<bool> Basic_Filter(const bool upper_body,const bool lower_body,const bool head_and_neck,const bool real_hands,const bool feet,const bool with_left_side)
{
    ARRAY<bool> filter(num_bones);
    if(upper_body){filter(BONE_THORAX)=filter(BONE_R_CLAVICLE)=filter(BONE_R_SCAPULA)=filter(BONE_R_HUMERUS)=filter(BONE_R_ULNA)=filter(BONE_R_RADIUS)=true;}
    if(lower_body){filter(BONE_HIP)=filter(BONE_R_PATELLA)=filter(BONE_R_FEMUR)=filter(BONE_R_TIBIA)=true;}
    if(lower_body && upper_body){filter(BONE_L1)=filter(BONE_L2)=filter(BONE_L3)=filter(BONE_L4)=filter(BONE_L5)=true;}
    if(real_hands){filter(BONE_R_1DISTPH)=filter(BONE_R_2DISTPH)=filter(BONE_R_3DISTPH)=filter(BONE_R_4DISTPH)=filter(BONE_R_5DISTPH)=true;
              filter(BONE_R_1MC)=filter(BONE_R_2MIDPH)=filter(BONE_R_3MIDPH)=filter(BONE_R_4MIDPH)=filter(BONE_R_5MIDPH)=true;
              filter(BONE_R_1PROXPH)=filter(BONE_R_2PROXPH)=filter(BONE_R_3PROXPH)=filter(BONE_R_4PROXPH)=filter(BONE_R_5PROXPH)=true;
              filter(BONE_R_PALM)=filter(BONE_R_WRIST)=true;}
    if(feet){filter(BONE_R_ANKLE)=filter(BONE_R_TOE_1)=filter(BONE_R_TOE_2)=filter(BONE_R_TOE_3)=filter(BONE_R_TOE_4)=filter(BONE_R_TOE_5)=true;}
    if(head_and_neck){filter(BONE_CRANIUM)=filter(BONE_ATLAS)=filter(BONE_AXIS)=filter(BONE_C3)=filter(BONE_C4)=filter(BONE_C5)=filter(BONE_C6)=filter(BONE_C7)=true;}
    if(with_left_side) Add_Reflected_Bones_To_Filter(filter);
    return filter;
}
//#####################################################################
// Function ARB_Joint_Id
//#####################################################################
JOINT_ID ARB_Joint_Id(const int local_joint_index) const
{
    return local_joint_index_to_arb_joint_id(local_joint_index);
}
//#####################################################################
// Function Create_Joint_Function
//#####################################################################
JOINT_FUNCTION<TV>* Create_Joint_Function(const int local_joint_index) const
{
    return arb->Create_Joint_Function(ARB_Joint_Id(local_joint_index));
}
//#####################################################################
// Function Rotated_Track
//#####################################################################
static INTERPOLATION_CURVE<T,ROTATION<TV> >* Rotated_Track(const INTERPOLATION_CURVE<T,ROTATION<TV> >& frame_track)
{
    PHYSBAM_FATAL_ERROR("Use interplation curve instead of track.");
#if 0
    INTERPOLATION_CURVE<T,ROTATION<TV> >* rotated_track=new INTERPOLATION_CURVE<T,ROTATION<TV> >(frame_track.trajectory.m,frame_track.time_grid.xmin,frame_track.time_grid.xmax);
    rotated_track->name=frame_track.name+"_rotated";
    rotated_track->periodic=frame_track.periodic;
    for(int i=0;i<frame_track.trajectory.m;i++){
        FRAME<TV> f=frame_track.trajectory(i);TV euler;f.r.Euler_Angles(euler.x,euler.y,euler.z);
        rotated_track->trajectory(i)=FRAME<TV>(TV(f.t.x,-f.t.y,-f.t.z),ROTATION<TV>::From_Euler_Angles(euler.x,-euler.y,-euler.z));}
#endif
    return &frame_track;
}
//#####################################################################
// Function Get_Double_Point_Joints_For_Bend_Joint
//#####################################################################
static void Get_Double_Point_Joints_For_Bend_Joint(const RIGID_BODY<TV>* parent,const RIGID_BODY<TV>* child,const JOINT<TV>* bend_joint,const T secondary_point_joint_offset,
    JOINT<TV>*& primary_point_joint,JOINT<TV>*& secondary_point_joint)
{
    assert(typeid(*bend_joint)==typeid(ANGLE_JOINT<TV>) && !((ANGLE_JOINT<TV>*)bend_joint)->constrain_angle);
    primary_point_joint=new POINT_JOINT<TV>;secondary_point_joint=new POINT_JOINT<TV>;
    primary_point_joint->frame_pj=bend_joint->frame_pj;primary_point_joint->frame_jp=bend_joint->frame_jp;
    primary_point_joint->frame_cj=bend_joint->frame_cj;primary_point_joint->frame_jc=bend_joint->frame_jc;
    primary_point_joint->J=bend_joint->J;primary_point_joint->J_inverse=bend_joint->J_inverse;
    primary_point_joint->joint_function=bend_joint->joint_function;primary_point_joint->name=bend_joint->name;
    primary_point_joint->global_post_stabilization=bend_joint->global_post_stabilization;
    primary_point_joint->primary_point_of_bend_joint=true;

    // location of secondary point is shifted along parent's x axis by some amount
    TV secondary_point_parent_location=bend_joint->F_pj()*TV(secondary_point_joint_offset,0,0);
    TV secondary_point_child_location=child->Frame().Inverse()*parent->Frame()*secondary_point_parent_location;
    secondary_point_joint->Set_Joint_To_Parent_Frame(FRAME<TV>(secondary_point_parent_location,bend_joint->frame_pj.r));
    secondary_point_joint->Set_Joint_To_Child_Frame(FRAME<TV>(secondary_point_child_location,bend_joint->frame_cj.r));
    // rest can be as usual I think
    secondary_point_joint->J=bend_joint->J;secondary_point_joint->J_inverse=bend_joint->J_inverse;
    secondary_point_joint->name=bend_joint->name+"_secondary_point";
    secondary_point_joint->global_post_stabilization=bend_joint->global_post_stabilization;
    secondary_point_joint->secondary_point_of_bend_joint=true;
}
//#####################################################################
// Function Add_Joint_With_Bend_Fix
//#####################################################################
JOINT_ID Add_Joint_With_Bend_Fix(const int parent_id,const int child_id,JOINT<TV>*& joint)
{
    if(use_double_point_joints_for_bend && typeid(*joint)==typeid(ANGLE_JOINT<TV>)){
        JOINT<TV> *primary_point_joint,*secondary_point_joint;
        RIGID_BODY<TV> *parent=&arb->rigid_body_collection.Rigid_Body(parent_id),*child=&arb->rigid_body_collection.Rigid_Body(child_id);
        Get_Double_Point_Joints_For_Bend_Joint(parent,child,joint,secondary_point_joint_offset,primary_point_joint,secondary_point_joint);
        delete joint;joint=primary_point_joint;
        arb->joint_mesh.Add_Articulation(parent_id,child_id,joint);
        if(verbose) LOG::cout<<"Adding secondary point joint for "<<joint->name<<std::endl;
        arb->joint_mesh.Add_Articulation(parent_id,child_id,secondary_point_joint);
        return secondary_point_joint->id_number;}
    arb->joint_mesh.Add_Articulation(parent_id,child_id,joint);
    return joint->id_number;
}
//#####################################################################
// Function Add_Joints_To_ARB
//#####################################################################
void Add_Joints_To_ARB()
{
    // Add joints to articulated rigid body (build local_joint_id_to_arb_joint_id)
    local_joint_index_to_arb_joint_id.Clean_Memory();
    local_joint_index_to_arb_joint_id.Resize(joint.m);
    for(int i=0;i<joint.m;i++) if(joint(i) && bones(joint_parent_and_child(i).x) && bones(joint_parent_and_child(i).y)){
        if(!bones(joint_parent_and_child(i).x)->particle_index || !bones(joint_parent_and_child(i).y)->particle_index){
            if(verbose) LOG::cout<<"Skipping joint "<<joint(i)->name<<std::endl;
            delete joint(i);joint(i)=0;}
        else{
            if(verbose) LOG::cout<<"Adding joint "<<joint(i)->name<<" to joint list"<<std::endl;
            local_joint_index_to_arb_joint_id(i)=Add_Joint_With_Bend_Fix(bones(joint_parent_and_child(i).x)->particle_index,bones(joint_parent_and_child(i).y)->particle_index,joint(i));}}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies(const ARRAY<bool>& bone_filter=ARRAY<bool>())
{
    Make_Bones(bone_filter);
    Make_Joints();
    Make_Muscles();

    Add_Joints_To_ARB();

    // Add muscles to muscle list
    for(int i=0;i<muscles.m;i++){
        bool has_all_rigid_bodies=true;
//         if(!muscles(i)->attachment_point_1 || !((muscles(i)->attachment_point_1))->particle_index) has_all_rigid_bodies=false;
//         if(!muscles(i)->attachment_point_2 || !((muscles(i)->attachment_point_2))->particle_index) has_all_rigid_bodies=false;
//         for(int j=1;j<=muscles(i)->via_points.m;j++) if(!muscles(i)->via_points(j)->particle_index) has_all_rigid_bodies=false;
        if(has_all_rigid_bodies){
            if(verbose) LOG::cout<<"Adding muscle "<<muscles(i)->name<<" to muscle list"<<std::endl;
            arb->muscle_list->Add_Muscle(muscles(i));}
        else{if(verbose) LOG::cout<<"Skipping muscle "<<muscles(i)->name<<std::endl;}}
}
//#####################################################################
// Function Make Bones
//#####################################################################
// makes bones in order of enum
void Make_Bones(const ARRAY<bool>& bone_filter=ARRAY<bool>())
{
    ARRAY<bool> the_filter=bone_filter;
    std::string bone_files[]={"sternum_ribs_back",   "hip",                "cranium",        "atlas",          "axis",
                              "vertebrae_c3",        "vertebrae_c4",       "vertebrae_c5",   "vertebrae_c6",   "vertebrae_c7",
                              "vertebrae_l1",        "vertebrae_l2",       "vertebrae_l3",   "vertebrae_l4",   "vertebrae_l5",
                              "clavicle_right",      "scapula_right",      "humerus_right",  "ulna_right",     "radius_right",
                              "1distph_right",       "2distph_right",      "3distph_right",  "4distph_right",  "5distph_right",
                              "1mc_right",           "2midph_right",       "3midph_right",   "4midph_right",   "5midph_right",
                              "1proxph_right",       "2proxph_right",      "3proxph_right",  "4proxph_right",  "5proxph_right",
                              "palm_right",          "wrist_right",        "patella_right",  "femur_right",    "tibia_right",
                              "foot_no_toes_right",  "toe_1_right",        "toe_2_right",    "toe_3_right",    "toe_4_right",
                              "toe_5_right",         "clavicle_left",      "scapula_left",   "humerus_left",   "ulna_left",
                              "radius_left",         "1distph_left",       "2distph_left",   "3distph_left",   "4distph_left",
                              "5distph_left",        "1mc_left",           "2midph_left",    "3midph_left",    "4midph_left",
                              "5midph_left",         "1proxph_left",       "2proxph_left",   "3proxph_left",   "4proxph_left",
                              "5proxph_left",        "palm_left",          "wrist_left",     "patella_left",   "femur_left",
                              "tibia_left",          "foot_no_toes_left",  "toe_1_left",     "toe_2_left",     "toe_3_left",
                              "toe_4_left",          "toe_5_left"};
    bones.Resize(num_bones);bone_names.Resize(num_bones);
    anatomical_frame_in_rigid_body_frame.Resize(num_bones); // only use these for the thorax and the 5 right side bones

    if(use_fused_tibia_and_patella){
        LOG::cout<<"Using fused tibia and patella"<<std::endl;
        if(!the_filter.m){the_filter.Resize(num_bones);the_filter.Fill(true);}
        the_filter(BONE_R_PATELLA)=the_filter(BONE_L_PATELLA)=false;
        bone_files[BONE_R_TIBIA-1]="tibia_and_patella_right";
        bone_files[BONE_L_TIBIA-1]="tibia_and_patella_left";}

    // Got these by looking at *original* bones that used to be in Rigid_Bodies/Visible_Human_Bones which were in an anatomical frame (rather than a rigid body frame
    // as they are now).  Joey found these anatomical frames for them.  Here they are given with respect to world frame (the variable is called in_rigid_body_frame
    // because we transform them to rigid body frame just below)
    anatomical_frame_in_rigid_body_frame(BONE_THORAX)=FRAME<TV>(TV((T)0.2931,(T)0.12740003,(T)1.5409998),ROTATION<TV>::From_Components((T)-0.0089113209,(T)-0.02153822,(T)0.0058329366,(T)0.99971116));
    anatomical_frame_in_rigid_body_frame(BONE_R_CLAVICLE)=FRAME<TV>(TV((T)0.26047158,(T)0.14116499,(T)1.531361),ROTATION<TV>::From_Components((T)0.17886177,(T)0.19062576,(T)-0.15840296,(T)0.95214438));
    anatomical_frame_in_rigid_body_frame(BONE_R_SCAPULA)=FRAME<TV>(TV((T)0.12737156,(T)0.18436509,(T)1.5960611),ROTATION<TV>::From_Components((T)-0.36945692,(T)0.16853584,(T)-0.11753589,(T)0.90624654));
    anatomical_frame_in_rigid_body_frame(BONE_R_HUMERUS)=FRAME<TV>(TV((T)0.11017161,(T)0.18825513,(T)1.555921),ROTATION<TV>::From_Components((T)-0.40632755,(T)0.073737398,(T)-0.16054413,(T)0.89648551));
    anatomical_frame_in_rigid_body_frame(BONE_R_ULNA)=FRAME<TV>(TV((T)0.047721587,(T)0.24119507,(T)1.2628108),ROTATION<TV>::From_Components((T)-0.3545996,(T)-0.049906723,(T)0.39113101,(T)0.84781182));
    anatomical_frame_in_rigid_body_frame(BONE_R_RADIUS)=FRAME<TV>(TV((T)0.030011594,(T)0.22594509,(T)1.267401),ROTATION<TV>::From_Components((T)0.84402317,(T)-0.21794307,(T)-0.38632673,(T)-0.30145824));

    for(int i=0;i<num_bones;i++){
        RIGID_BODY<TV>* rigid_body=0;
        std::string filename=data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/"+bone_files[i-1];
        if(!the_filter.Valid_Index(i) || the_filter(i)){
            if(verbose) LOG::cout<<"Reading rigid body "<<filename<<std::endl;
            int id=arb->rigid_body_collection.Add_Rigid_Body(stream_type,filename,(T)1,true,true,false);
            rigid_body=&arb->rigid_body_collection.Rigid_Body(id);}
        else{
            if(verbose) LOG::cout<<"Reading rigid body (but not adding to rigid body list) "<<filename<<std::endl;
            FRAME<TV> frame;
            rigid_body=new RIGID_BODY<TV>(rigid_body_collection);FILE_UTILITIES::Read_From_File(stream_type,filename+".rgd",rigid_body->Mass(),frame);
            rigid_body->Set_Frame(frame);}
        rigid_body->Set_Coefficient_Of_Restitution((T).5);
        rigid_body->Set_Coefficient_Of_Friction((T).5);
        bone_names(i)=bone_files[i-1];
        if(bone_files[i-1]=="tibia_and_patella_right") bone_names(i)="tibia_right";
        else if(bone_files[i-1]=="tibia_and_patella_left") bone_names(i)="tibia_left";
        rigid_body->Set_Name(bone_files[i-1]);
        // get anatomical frame w.r.t. rigid body frame instead of world frame
        anatomical_frame_in_rigid_body_frame(i)=rigid_body->Frame().Inverse()*anatomical_frame_in_rigid_body_frame(i);
        rigid_body->Set_Frame(transform*rigid_body->Frame());
        bones(i)=rigid_body;}
}
//#####################################################################
// Function Individual_Side_Bones
//#####################################################################
static int Individual_Side_Bones()
{
    return ((int)BONE_END-1-medial_bones)/2;
}
//#####################################################################
// Function Individual_Side_Joints
//#####################################################################
static int Individual_Side_Joints()
{
    return ((int)JOINT_END-1-medial_joints)/2;
}
//#####################################################################
// Function Reflected_Bone
//#####################################################################
static int Reflected_Bone(const int bone_id)
{
    int individual_side_bones=Individual_Side_Bones();
    if(bone_id<=medial_bones) return bone_id;
    else if(bone_id<=medial_bones+individual_side_bones) return bone_id+individual_side_bones;
    else return bone_id-individual_side_bones;
}
//#####################################################################
// Function Reflected_Joint
//#####################################################################
static int Reflected_Joint(const int joint_id)
{
    int individual_side_joints=Individual_Side_Joints();
    if(joint_id<=medial_joints) return joint_id;

    else if(joint_id<=medial_joints+individual_side_joints) return joint_id+individual_side_joints;
    else return joint_id-individual_side_joints;
}
//#####################################################################
// Function Reflection_Matrix
//#####################################################################
// returns reflection matrix to get to this bone from its right counterpart (works for medial bones too)
MATRIX<T,4> Reflection_Matrix(const int bone_id,int specific_side=0)
{
    if(bone_id==BONE_THORAX || (BONE_CRANIUM<=bone_id && bone_id<=BONE_L5) || (BONE_L_CLAVICLE<=bone_id && bone_id<=BONE_L_WRIST) || specific_side==1) return thorax_left_side_reflection;
    else{assert(bone_id==BONE_HIP || (BONE_L_PATELLA<=bone_id && bone_id<=BONE_L_TOE_5) || specific_side==2);return pelvis_left_side_reflection;}
}
//#####################################################################
// Function Rotation_From_Axis
//#####################################################################
// Figures out some orientation (quaternion) given one axis which will be the x axis, and another axis used to create an orthogonal frame
static ROTATION<TV> Rotation_From_Axis(const TV& axis,const TV& other_axis)
{
    TV x=axis.Normalized(),y=(other_axis.Projected_Orthogonal_To_Unit_Direction(x)).Normalized(),z=TV::Cross_Product(x,y);
    return ROTATION<TV>(MATRIX<T,3>(x,y,z));
}
//#####################################################################
// Function Reflected_Frame
//#####################################################################
// Reflects a frame from the right bones (which appear on his left side) to the left bones (which appear on the right side)
// Essentially this reflection flips the x axis (so we stay a LHS) but the other axes are reflected
static FRAME<TV> Reflected_Frame(const FRAME<TV>& frame,const MATRIX<T,4>& reflection)
{
    MATRIX<T,3> r=frame.r.Rotation_Matrix();r.Column(1)*=-1;
    return FRAME<TV>(reflection.Homogeneous_Times(frame.t),ROTATION<TV>(reflection.Extract_Rotation()*r));
}
//#####################################################################
// Function Reflected_Frame_Relative_To_Body
//#####################################################################
// Given a left side bone_id, returns reflected joint frame relative to that bone (using given frame relative to right bone)
FRAME<TV> Reflected_Frame_Relative_To_Body(const int bone_id,const FRAME<TV>& frame)
{
    return bones(bone_id)->Frame().Inverse()*transform*Reflected_Frame(transform.Inverse()*bones(Reflected_Bone(bone_id))->Frame()*frame,Reflection_Matrix(bone_id));
}
//#####################################################################
// Function Reflected_Frame_Relative_To_Body_From_Given_Bone
//#####################################################################
// Given a bone, returns reflected joint frame relative to that bone (using given frame relative to right bone)
FRAME<TV> Reflected_Frame_Relative_To_Body_From_Given_Bone(const RIGID_BODY<TV>* this_bone,const RIGID_BODY<TV>* opposite_bone,const FRAME<TV>& frame,const int specific_side=0)
{
    return this_bone->Frame().Inverse()*transform*Reflected_Frame(transform.Inverse()*opposite_bone->Frame()*frame,Reflection_Matrix(0,specific_side));
}
//#####################################################################
// Function Add_Reflected_Left_Joint
//#####################################################################
// Given parent and child (which are left bones) and joint_id (the left joint's id) we compute the reflected transforms from the right side counterpart
void Add_Reflected_Left_Joint(const int local_left_joint_id)
{
    int local_right_joint_id=Reflected_Joint(local_left_joint_id);
    JOINT<TV>* right_joint=joint(local_right_joint_id);if(!right_joint) return; // joint doesn't exist
    int left_parent_bone=Reflected_Bone(joint_parent_and_child(local_right_joint_id).x),
        left_child_bone=Reflected_Bone(joint_parent_and_child(local_right_joint_id).y);
    JOINT<TV>* new_joint=0;
    assert(!right_joint->primary_point_of_bend_joint);
    if(typeid(*right_joint)==typeid(POINT_JOINT<TV>)) new_joint=new POINT_JOINT<TV>();
    else if(typeid(*right_joint)==typeid(ANGLE_JOINT<TV>)) new_joint=new ANGLE_JOINT<TV>();
    else PHYSBAM_NOT_IMPLEMENTED();
    if(Add_Joint_Safe(left_parent_bone,left_child_bone,local_left_joint_id,new_joint)){
        joint(local_left_joint_id)->Set_Name(joint(local_right_joint_id)->name+"_left");
        if(verbose) LOG::cout<<"Reflecting joint "<<joint(local_right_joint_id)->name<<" -> "<<joint(local_left_joint_id)->name<<std::endl;
        FRAME<TV> reflected_parent_frame=Reflected_Frame_Relative_To_Body(left_parent_bone,right_joint->frame_pj);
        FRAME<TV> reflected_child_frame=Reflected_Frame_Relative_To_Body(left_child_bone,right_joint->frame_cj);
        joint(local_left_joint_id)->Set_Joint_To_Parent_Frame(reflected_parent_frame);
        joint(local_left_joint_id)->Set_Joint_To_Child_Frame(reflected_child_frame);}
}
//#####################################################################
// Function Update_Reflected_Left_Joint
//#####################################################################
void Update_Reflected_Left_Joint(const int local_left_joint_id)
{
    int local_right_joint_id=Reflected_Joint(local_left_joint_id);
    JOINT<TV>* right_joint=joint(local_right_joint_id);if(!right_joint) return; // joint doesn't exist
    int left_parent_bone=Reflected_Bone(joint_parent_and_child(local_right_joint_id).x),
        left_child_bone=Reflected_Bone(joint_parent_and_child(local_right_joint_id).y);
    assert(!right_joint->primary_point_of_bend_joint);
    if(joint(local_left_joint_id)){
        joint(local_left_joint_id)->Set_Name(joint(local_right_joint_id)->name+"_left");
        if(verbose) LOG::cout<<"Reflecting joint "<<joint(local_right_joint_id)->name<<" -> "<<joint(local_left_joint_id)->name<<std::endl;
        FRAME<TV> reflected_parent_frame=Reflected_Frame_Relative_To_Body(left_parent_bone,right_joint->frame_pj);
        FRAME<TV> reflected_child_frame=Reflected_Frame_Relative_To_Body(left_child_bone,right_joint->frame_cj);
        joint(local_left_joint_id)->Set_Joint_To_Parent_Frame(reflected_parent_frame);
        joint(local_left_joint_id)->Set_Joint_To_Child_Frame(reflected_child_frame);}
}
//#####################################################################
// Function Add_Additional_Reflected_Left_Joint
//#####################################################################
JOINT_ID Add_Additional_Reflected_Left_Joint(const JOINT_ID right_joint_id,RIGID_BODY<TV>* parent,RIGID_BODY<TV>* child,const int specific_side)
{
    JOINT<TV>* right_joint=arb->joint_mesh(right_joint_id);
    RIGID_BODY<TV> *left_parent_bone=parent,*left_child_bone=child;
    RIGID_BODY<TV> *right_parent_bone=arb->Parent(right_joint_id),*right_child_bone=arb->Child(right_joint_id);
    JOINT<TV>* new_joint=0;
    if(right_joint->joint_type==JOINT<TV>::TYPE_ANGLE_JOINT || right_joint->primary_point_of_bend_joint) new_joint=new ANGLE_JOINT<TV>();
    else if(right_joint->joint_type==JOINT<TV>::TYPE_POINT_JOINT) new_joint=new POINT_JOINT<TV>();
    else assert(false);

    new_joint->Set_Name(right_joint->name+"_left");
    if(verbose) LOG::cout<<"Reflecting joint "<<right_joint->name<<" -> "<<new_joint->name<<std::endl;
    FRAME<TV> reflected_parent_frame=Reflected_Frame_Relative_To_Body_From_Given_Bone(left_parent_bone,right_parent_bone,right_joint->frame_pj,specific_side);
    FRAME<TV> reflected_child_frame=Reflected_Frame_Relative_To_Body_From_Given_Bone(left_child_bone,right_child_bone,right_joint->frame_cj,specific_side);
    new_joint->Set_Joint_To_Parent_Frame(reflected_parent_frame);
    new_joint->Set_Joint_To_Child_Frame(reflected_child_frame);
    return Add_Joint_With_Bend_Fix(left_parent_bone->particle_index,left_child_bone->particle_index,new_joint);
}
//#####################################################################
// Function Add_Joint_Safe2
//#####################################################################
bool Add_Joint_Safe2(const int local_parent_index,const int local_child_index,const int local_joint_index,JOINT<TV>* new_joint)
{
    JOINT<TV>* new_joint_tmp=new_joint;return Add_Joint_Safe(local_parent_index,local_child_index,local_joint_index,new_joint_tmp);
}
//#####################################################################
// Function Add_Joint_Safe
//#####################################################################
bool Add_Joint_Safe(const int local_parent_index,const int local_child_index,const int local_joint_index,JOINT<TV>*& new_joint)
{
    if(bones(local_parent_index) && bones(local_child_index)){
        if(use_only_point_joints && typeid(*new_joint)!=typeid(POINT_JOINT<TV>)){delete new_joint;new_joint=new POINT_JOINT<TV>();}
        joint(local_joint_index)=new_joint;
        joint_parent_and_child(local_joint_index)=PAIR<int,int>(local_parent_index,local_child_index);
        return true;}
    else{delete new_joint;new_joint=0;return false;}
}
//#####################################################################
// Function Get_Joint_Orientation_In_Parent_Anatomical_Frame
//#####################################################################
ROTATION<TV> Get_Joint_Orientation_In_Parent_Anatomical_Frame(const int parent_index,const int child_index,const T rot1_deg,const T rot2_deg,const T rot3_deg,
    const ROTATION<TV>& joint_orientation_in_child_frame)
{
    ROTATION<TV> anatomical_to_world_frame_parent=bones(parent_index)->Rotation()*anatomical_frame_in_rigid_body_frame(parent_index).r;
    ROTATION<TV> anatomical_to_world_frame_child=bones(child_index)->Rotation()*anatomical_frame_in_rigid_body_frame(child_index).r;
    ROTATION<TV> rotation=ROTATION<TV>(rot1_deg*(T)(pi/180),TV(1,0,0))*ROTATION<TV>(rot2_deg*(T)(pi/180),TV(0,1,0))*ROTATION<TV>(rot3_deg*(T)(pi/180),TV(0,0,1));
    return anatomical_to_world_frame_parent.Inverse()*anatomical_to_world_frame_child*joint_orientation_in_child_frame*rotation.Inverse();
}
//#####################################################################
// Function Make_Joints
//#####################################################################
void Make_Joints()
{
    std::string joint_names[]={"joint_atlas","joint_axis","joint_c3","joint_c4","joint_c5","joint_c6","joint_c7","joint_c7_thorax","joint_l1","joint_l2",
                               "joint_l3","joint_l4","joint_l5","joint_l5_coccyx",
                               "joint_r_sternoclavicular","joint_r_acromioclavicular","joint_r_glenohumeral","joint_r_humeroulnar","joint_r_radioulnar","joint_r_radiocarpal",
                               "joint_r_hip","joint_r_knee","joint_r_patella","joint_r_ankle","joint_r_toe_1","joint_r_toe_2","joint_r_toe_3","joint_r_toe_4","joint_r_toe_5",
                               "joint_l_sternoclavicular","joint_l_acromioclavicular","joint_l_glenohumeral","joint_l_humeroulnar","joint_l_radioulnar","joint_l_radiocarpal",
                               "joint_l_hip","joint_l_knee","joint_l_patella","joint_l_ankle","joint_l_toe_1","joint_l_toe_2","joint_l_toe_3","joint_l_toe_4","joint_l_toe_5"};

    bool do_left=true;
    FRAME<TV> child_frame,parent_frame;
    int num_joints=medial_joints+(do_left?2:1)*Individual_Side_Joints();
    joint.Clean_Memory();joint_parent_and_child.Clean_Memory();
    joint.Resize(num_joints);joint_parent_and_child.Resize(num_joints);

    if(Add_Joint_Safe2(BONE_CRANIUM,BONE_ATLAS,JOINT_ATLAS,new POINT_JOINT<TV>())){
        TV atlas_joint_center;
        joint(JOINT_ATLAS)->Set_Joint_To_Parent_Frame(bones(BONE_CRANIUM)->Frame().Inverse()*bones(BONE_ATLAS)->Frame()*FRAME<TV>(atlas_joint_center));
        joint(JOINT_ATLAS)->Set_Joint_To_Child_Frame(FRAME<TV>(atlas_joint_center));}

    if(Add_Joint_Safe2(BONE_ATLAS,BONE_AXIS,JOINT_AXIS,new POINT_JOINT<TV>())){
        TV axis_joint_center((T)0.00414979,(T)-0.0112626,(T)0.00800778);
        joint(JOINT_AXIS)->Set_Joint_To_Parent_Frame(bones(BONE_ATLAS)->Frame().Inverse()*bones(BONE_AXIS)->Frame()*FRAME<TV>(axis_joint_center));
        joint(JOINT_AXIS)->Set_Joint_To_Child_Frame(FRAME<TV>(axis_joint_center));}

    if(Add_Joint_Safe2(BONE_AXIS,BONE_C3,JOINT_C3,new POINT_JOINT<TV>())){
        TV c3_joint_center(0,0,(T)0.00673424);
        joint(JOINT_C3)->Set_Joint_To_Parent_Frame(bones(BONE_AXIS)->Frame().Inverse()*bones(BONE_C3)->Frame()*FRAME<TV>(c3_joint_center));
        joint(JOINT_C3)->Set_Joint_To_Child_Frame(FRAME<TV>(c3_joint_center));}

    if(Add_Joint_Safe2(BONE_C3,BONE_C4,JOINT_C4,new POINT_JOINT<TV>())){
        TV c4_joint_center(0,0,(T)0.00585404);
        joint(JOINT_C4)->Set_Joint_To_Parent_Frame(bones(BONE_C3)->Frame().Inverse()*bones(BONE_C4)->Frame()*FRAME<TV>(c4_joint_center));
        joint(JOINT_C4)->Set_Joint_To_Child_Frame(FRAME<TV>(c4_joint_center));}

    if(Add_Joint_Safe2(BONE_C4,BONE_C5,JOINT_C5,new POINT_JOINT<TV>())){
        TV c5_joint_center(0,0,(T)0.00573239);
        joint(JOINT_C5)->Set_Joint_To_Parent_Frame(bones(BONE_C4)->Frame().Inverse()*bones(BONE_C5)->Frame()*FRAME<TV>(c5_joint_center));
        joint(JOINT_C5)->Set_Joint_To_Child_Frame(FRAME<TV>(c5_joint_center));}

    if(Add_Joint_Safe2(BONE_C5,BONE_C6,JOINT_C6,new POINT_JOINT<TV>())){
        TV c6_joint_center(0,(T)-0.00139795,(T)0.00569737);
        joint(JOINT_C6)->Set_Joint_To_Parent_Frame(bones(BONE_C5)->Frame().Inverse()*bones(BONE_C6)->Frame()*FRAME<TV>(c6_joint_center));
        joint(JOINT_C6)->Set_Joint_To_Child_Frame(FRAME<TV>(c6_joint_center));}

    if(Add_Joint_Safe2(BONE_C6,BONE_C7,JOINT_C7,new POINT_JOINT<TV>())){
        TV c7_joint_center((T)0.00157246,(T)-0.004641241,(T)0.00609826);
        joint(JOINT_C7)->Set_Joint_To_Parent_Frame(bones(BONE_C6)->Frame().Inverse()*bones(BONE_C7)->Frame()*FRAME<TV>(c7_joint_center));
        joint(JOINT_C7)->Set_Joint_To_Child_Frame(FRAME<TV>(c7_joint_center));}

    if(Add_Joint_Safe2(BONE_C7,BONE_THORAX,JOINT_C7_THORAX,new POINT_JOINT<TV>())){
        TV c7_thorax_joint_center((T)0.000384506,(T)-0.00167134,(T)-0.0105332);
        joint(JOINT_C7_THORAX)->Set_Joint_To_Parent_Frame(FRAME<TV>(c7_thorax_joint_center));
        joint(JOINT_C7_THORAX)->Set_Joint_To_Child_Frame(bones(BONE_THORAX)->Frame().Inverse()*bones(BONE_C7)->Frame()*FRAME<TV>(c7_thorax_joint_center));}

    if(Add_Joint_Safe2(BONE_THORAX,BONE_L1,JOINT_L1,new POINT_JOINT<TV>())){
        TV l1_joint_center(0,(T)0.0255006,(T)0.0160207);
        joint(JOINT_L1)->Set_Joint_To_Parent_Frame(bones(BONE_THORAX)->Frame().Inverse()*bones(BONE_L1)->Frame()*FRAME<TV>(l1_joint_center));
        joint(JOINT_L1)->Set_Joint_To_Child_Frame(FRAME<TV>(l1_joint_center));}

    if(Add_Joint_Safe2(BONE_L1,BONE_L2,JOINT_L2,new POINT_JOINT<TV>())){
        TV l2_joint_center(0,(T)0.0253012,(T)0.0172893);
        joint(JOINT_L2)->Set_Joint_To_Parent_Frame(bones(BONE_L1)->Frame().Inverse()*bones(BONE_L2)->Frame()*FRAME<TV>(l2_joint_center));
        joint(JOINT_L2)->Set_Joint_To_Child_Frame(FRAME<TV>(l2_joint_center));}

    if(Add_Joint_Safe2(BONE_L2,BONE_L3,JOINT_L3,new POINT_JOINT<TV>())){
        TV l3_joint_center(0,(T)0.0247191,(T)0.0169881);
        joint(JOINT_L3)->Set_Joint_To_Parent_Frame(bones(BONE_L2)->Frame().Inverse()*bones(BONE_L3)->Frame()*FRAME<TV>(l3_joint_center));
        joint(JOINT_L3)->Set_Joint_To_Child_Frame(FRAME<TV>(l3_joint_center));}

    if(Add_Joint_Safe2(BONE_L3,BONE_L4,JOINT_L4,new POINT_JOINT<TV>())){
        TV l4_joint_center(0,(T)0.0251905,(T)0.0161453);
        joint(JOINT_L4)->Set_Joint_To_Parent_Frame(bones(BONE_L3)->Frame().Inverse()*bones(BONE_L4)->Frame()*FRAME<TV>(l4_joint_center));
        joint(JOINT_L4)->Set_Joint_To_Child_Frame(FRAME<TV>(l4_joint_center));}

    if(Add_Joint_Safe2(BONE_L4,BONE_L5,JOINT_L5,new POINT_JOINT<TV>())){
        TV l5_joint_center((T)0.00722887,(T)0.0232633,(T)0.0136113);
        joint(JOINT_L5)->Set_Joint_To_Parent_Frame(bones(BONE_L4)->Frame().Inverse()*bones(BONE_L5)->Frame()*FRAME<TV>(l5_joint_center));
        joint(JOINT_L5)->Set_Joint_To_Child_Frame(FRAME<TV>(l5_joint_center));}

    if(Add_Joint_Safe2(BONE_L5,BONE_HIP,JOINT_L5_COCCYX,new POINT_JOINT<TV>())){
        TV l5_coccyx_joint_center((T)0.0053103,(T)0.0261411,(T)-0.0161265);
        joint(JOINT_L5_COCCYX)->Set_Joint_To_Parent_Frame(FRAME<TV>(l5_coccyx_joint_center));
        joint(JOINT_L5_COCCYX)->Set_Joint_To_Child_Frame(bones(BONE_HIP)->Frame().Inverse()*bones(BONE_L5)->Frame()*FRAME<TV>(l5_coccyx_joint_center));}

    // Sternoclavicular - gimbal joint (ball-and-socket)
    if(Add_Joint_Safe2(BONE_THORAX,BONE_R_CLAVICLE,JOINT_R_STERNOCLAVICULAR,new POINT_JOINT<TV>())){
        TV SC_origin_in_anatomical_thorax_frame((T)0.03254,(T)-0.01438,(T)-0.00801);
        TV SC_origin_in_anatomical_clavicle_frame;
        ROTATION<TV> SC_orientation_in_anatomical_clavicle_frame(MATRIX<T,3>(0,-1,0,0,0,1,-1,0,0));
        ROTATION<TV> SC_orientation_in_anatomical_thorax_frame=Get_Joint_Orientation_In_Parent_Anatomical_Frame(BONE_THORAX,BONE_R_CLAVICLE,0,0,0,SC_orientation_in_anatomical_clavicle_frame);
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_THORAX)*FRAME<TV>(SC_origin_in_anatomical_thorax_frame,SC_orientation_in_anatomical_thorax_frame);
        joint(JOINT_R_STERNOCLAVICULAR)->Set_Joint_To_Parent_Frame(parent_frame); //thorax
        child_frame=anatomical_frame_in_rigid_body_frame(BONE_R_CLAVICLE)*FRAME<TV>(SC_origin_in_anatomical_clavicle_frame,SC_orientation_in_anatomical_clavicle_frame);
        joint(JOINT_R_STERNOCLAVICULAR)->Set_Joint_To_Child_Frame(child_frame);} //clavicle

    // Acromioclavicular - gimbal joint
    if(Add_Joint_Safe2(BONE_R_CLAVICLE,BONE_R_SCAPULA,JOINT_R_ACROMIOCLAVICULAR,new POINT_JOINT<TV>())){
        TV AC_origin_in_anatomical_clavicle_frame(TV((T)0.15417,0,0));
        TV AC_origin_in_anatomical_scapula_frame;
        ROTATION<TV> AC_orientation_in_anatomical_scapula_frame(MATRIX<T,3>(0,-1,0,0,0,1,-1,0,0));
        ROTATION<TV> AC_orientation_in_anatomical_clavicle_frame=Get_Joint_Orientation_In_Parent_Anatomical_Frame(BONE_R_CLAVICLE,BONE_R_SCAPULA,0,0,0,
            AC_orientation_in_anatomical_scapula_frame);
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_R_CLAVICLE)*FRAME<TV>(AC_origin_in_anatomical_clavicle_frame,AC_orientation_in_anatomical_clavicle_frame);
        joint(JOINT_R_ACROMIOCLAVICULAR)->Set_Joint_To_Parent_Frame(parent_frame); //clavicle
        child_frame=anatomical_frame_in_rigid_body_frame(BONE_R_SCAPULA)*FRAME<TV>(AC_origin_in_anatomical_scapula_frame,AC_orientation_in_anatomical_scapula_frame);
        joint(JOINT_R_ACROMIOCLAVICULAR)->Set_Joint_To_Child_Frame(child_frame);} //scapula

    // Glenohumeral - gimbal joint
    if(Add_Joint_Safe2(BONE_R_SCAPULA,BONE_R_HUMERUS,JOINT_R_GLENOHUMERAL,new POINT_JOINT<TV>())){
        TV GH_origin_in_anatomical_scapula_frame(TV(0,0,(T)-0.04384));
        TV GH_origin_in_anatomical_humerus_frame;
        ROTATION<TV> GH_orientation_in_anatomical_humerus_frame(MATRIX<T,3>(1,0,0,0,-1,0,0,0,-1));
        ROTATION<TV> GH_orientation_in_anatomical_scapula_frame=Get_Joint_Orientation_In_Parent_Anatomical_Frame(BONE_R_SCAPULA,BONE_R_HUMERUS,(T)-14.799,(T)17.371,(T)-45.266,
            GH_orientation_in_anatomical_humerus_frame);
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_R_SCAPULA)*FRAME<TV>(GH_origin_in_anatomical_scapula_frame,GH_orientation_in_anatomical_scapula_frame);
        joint(JOINT_R_GLENOHUMERAL)->Set_Joint_To_Parent_Frame(parent_frame); //scapula
        child_frame=anatomical_frame_in_rigid_body_frame(BONE_R_HUMERUS)*FRAME<TV>(GH_origin_in_anatomical_humerus_frame,GH_orientation_in_anatomical_humerus_frame);
        joint(JOINT_R_GLENOHUMERAL)->Set_Joint_To_Child_Frame(child_frame);} //humerus

    // Humeroulnar - flexion/extension
    if(Add_Joint_Safe2(BONE_R_HUMERUS,BONE_R_ULNA,JOINT_R_HUMEROULNAR,new ANGLE_JOINT<TV>())){
        TV HU_FE_origin_in_anatomical_humerus_frame((T).00081,(T).02518,-(T).30329);
        TV HU_FE_origin_in_anatomical_ulna_frame;
        ROTATION<TV> HU_FE_orientation_in_anatomical_ulna_frame=Rotation_From_Axis(TV((T).976409,0,(T)-.215927),TV(0,1,0));
        ROTATION<TV> HU_FE_orientation_in_anatomical_humerus_frame=Get_Joint_Orientation_In_Parent_Anatomical_Frame(BONE_R_HUMERUS,BONE_R_ULNA,(T)64.727,0,0,
            HU_FE_orientation_in_anatomical_ulna_frame);
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_R_HUMERUS)*FRAME<TV>(HU_FE_origin_in_anatomical_humerus_frame,HU_FE_orientation_in_anatomical_humerus_frame);
        joint(JOINT_R_HUMEROULNAR)->Set_Joint_To_Parent_Frame(parent_frame); //humerus
        child_frame=anatomical_frame_in_rigid_body_frame(BONE_R_ULNA)*FRAME<TV>(HU_FE_origin_in_anatomical_ulna_frame,HU_FE_orientation_in_anatomical_ulna_frame);
        joint(JOINT_R_HUMEROULNAR)->Set_Joint_To_Child_Frame(child_frame);} //ulna

    // Radioulnar - pronation/supination
    if(Add_Joint_Safe2(BONE_R_ULNA,BONE_R_RADIUS,JOINT_R_RADIOULNAR,new ANGLE_JOINT<TV>())){
        TV RU_PS_origin_in_anatomical_ulna_frame((T).02382,0,0);
        TV RU_PS_origin_in_anatomical_radius_frame;
        ROTATION<TV> RU_PS_orientation_in_anatomical_radius_frame=Rotation_From_Axis(TV((T).135789,(T)0.023796,(T).990452),TV(0,1,0));
        ROTATION<TV> RU_PS_orientation_in_anatomical_ulna_frame=Get_Joint_Orientation_In_Parent_Anatomical_Frame(BONE_R_ULNA,BONE_R_RADIUS,(T)92.096,0,0,
            RU_PS_orientation_in_anatomical_radius_frame);
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_R_ULNA)*FRAME<TV>(RU_PS_origin_in_anatomical_ulna_frame,RU_PS_orientation_in_anatomical_ulna_frame);
        joint(JOINT_R_RADIOULNAR)->Set_Joint_To_Parent_Frame(parent_frame); //ulna
        child_frame=anatomical_frame_in_rigid_body_frame(BONE_R_RADIUS)*FRAME<TV>(RU_PS_origin_in_anatomical_radius_frame,RU_PS_orientation_in_anatomical_radius_frame);
        joint(JOINT_R_RADIOULNAR)->Set_Joint_To_Child_Frame(child_frame);} //radius

    // Radiocarpal - flexion/extension
    if(Add_Joint_Safe2(BONE_R_RADIUS,BONE_R_WRIST,JOINT_R_RADIOCARPAL,new ANGLE_JOINT<TV>())){
        TV RC_FE_origin_in_anatomical_radius_frame((T)-.01834,(T)-.00034,(T)-.27811);
        TV RC_FE_axis_in_anatomical_radius_frame((T)0.966275,0,(T)-0.257513);
        FRAME<TV> RC_FE_in_anatomical_radius_frame(RC_FE_origin_in_anatomical_radius_frame,Rotation_From_Axis(RC_FE_axis_in_anatomical_radius_frame,TV(0,1,0)));
        parent_frame=anatomical_frame_in_rigid_body_frame(BONE_R_RADIUS)*RC_FE_in_anatomical_radius_frame;
        joint(JOINT_R_RADIOCARPAL)->Set_Joint_To_Parent_Frame(parent_frame); //radius
        //NOTE: we don't have anatomical frame for hand -- instead we assume the initial configuration of the hand is the reference configuration
        child_frame=bones(BONE_R_WRIST)->Frame().Inverse()*bones(BONE_R_RADIUS)->Frame()*parent_frame;
        joint(JOINT_R_RADIOCARPAL)->Set_Joint_To_Child_Frame(child_frame);} //wrist/hand

    // TODO: add radiocarpal radial/ulnar deviation between wrist and hand

    // Hip Joint -- ball and socket
    if(Add_Joint_Safe2(BONE_HIP,BONE_R_FEMUR,JOINT_R_HIP,new POINT_JOINT<TV>())){
        TV joint_center2((T)0.0263956,(T)-0.00229903,(T)0.235911);
        joint(JOINT_R_HIP)->Set_Joint_To_Parent_Frame(bones(BONE_HIP)->Frame().Inverse()*bones(BONE_R_FEMUR)->Frame()*FRAME<TV>(joint_center2));
        joint(JOINT_R_HIP)->Set_Joint_To_Child_Frame(FRAME<TV>(joint_center2));}

    // Knee Joint - flexion/extension
    if(Add_Joint_Safe2(BONE_R_FEMUR,BONE_R_TIBIA,JOINT_R_KNEE,new ANGLE_JOINT<TV>())){
        TV tibia_joint_center((T)0.00293271,(T)0.0138317,(T)-0.19082);
        TV tibia_axis_of_rotation((T)(-0.0413611-0.0426697),(T)(0.0110518-0.0158082),(T)(-0.189065+0.198578));
        FRAME<TV> knee_frame_in_parent(tibia_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),tibia_axis_of_rotation));
        joint(JOINT_R_KNEE)->Set_Joint_To_Parent_Frame(knee_frame_in_parent);
        joint(JOINT_R_KNEE)->Set_Joint_To_Child_Frame(bones(BONE_R_TIBIA)->Frame().Inverse()*bones(BONE_R_FEMUR)->Frame()*knee_frame_in_parent);}

    // Patella Joint - bend
    if(Add_Joint_Safe2(BONE_R_FEMUR,BONE_R_PATELLA,JOINT_R_PATELLA,new POINT_JOINT<TV>())){
        TV patella_joint_center;
        joint(JOINT_R_PATELLA)->Set_Joint_To_Parent_Frame(bones(BONE_R_FEMUR)->Frame().Inverse()*bones(BONE_R_PATELLA)->Frame()*FRAME<TV>(patella_joint_center));
        joint(JOINT_R_PATELLA)->Set_Joint_To_Child_Frame(FRAME<TV>(patella_joint_center));}

    // Ankle Joint - ball and socket -ish
    if(Add_Joint_Safe2(BONE_R_TIBIA,BONE_R_ANKLE,JOINT_R_ANKLE,new POINT_JOINT<TV>())){
        TV ankle_joint_center((T)0.00380192,(T)0.00355539,(T)-0.259706);
        joint(JOINT_R_ANKLE)->Set_Joint_To_Parent_Frame(FRAME<TV>(ankle_joint_center));
        joint(JOINT_R_ANKLE)->Set_Joint_To_Child_Frame(bones(BONE_R_ANKLE)->Frame().Inverse()*bones(BONE_R_TIBIA)->Frame()*FRAME<TV>(ankle_joint_center));}

    // Big toe
    if(Add_Joint_Safe2(BONE_R_ANKLE,BONE_R_TOE_1,JOINT_R_TOE_1,new ANGLE_JOINT<TV>())){
        TV toe_1_joint_center((T)0.0335615,(T)-0.0907241,(T)-0.00454778);
        TV toe_1_axis_of_rotation((T)(0.0446287-0.0240753),(T)(-0.0891431+0.0986293),(T)(-0.00612882+0.00296675));
        FRAME<TV> toe_1_frame_in_parent(toe_1_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),toe_1_axis_of_rotation));
        joint(JOINT_R_TOE_1)->Set_Joint_To_Parent_Frame(toe_1_frame_in_parent);
        joint(JOINT_R_TOE_1)->Set_Joint_To_Child_Frame(bones(BONE_R_TOE_1)->Frame().Inverse()*bones(BONE_R_ANKLE)->Frame()*toe_1_frame_in_parent);}

    // Pointer toe
    if(Add_Joint_Safe2(BONE_R_ANKLE,BONE_R_TOE_2,JOINT_R_TOE_2,new ANGLE_JOINT<TV>())){
        TV toe_2_joint_center((T)0.0130081,(T)-0.101791,(T)-0.014034);
        TV toe_2_axis_of_rotation((T)(0.00984601-0.0177512),(T)(-0.104953+0.103372),(T)(-0.0171961+0.0108719));
        FRAME<TV> toe_2_frame_in_parent(toe_2_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),toe_2_axis_of_rotation));
        joint(JOINT_R_TOE_2)->Set_Joint_To_Parent_Frame(toe_2_frame_in_parent);
        joint(JOINT_R_TOE_2)->Set_Joint_To_Child_Frame(bones(BONE_R_TOE_2)->Frame().Inverse()*bones(BONE_R_ANKLE)->Frame()*toe_2_frame_in_parent);}

    // Middle toe
    if(Add_Joint_Safe2(BONE_R_ANKLE,BONE_R_TOE_3,JOINT_R_TOE_3,new ANGLE_JOINT<TV>())){
        TV toe_3_joint_center((T)0.00194084,(T)-0.0986293,(T)-0.0219392);
        TV toe_3_axis_of_rotation((T)(-0.00122122-0.00668395),0,(T)(-0.0266823+0.0187771));
        FRAME<TV> toe_3_frame_in_parent(toe_3_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),toe_3_axis_of_rotation));
        joint(JOINT_R_TOE_3)->Set_Joint_To_Parent_Frame(toe_3_frame_in_parent);
        joint(JOINT_R_TOE_3)->Set_Joint_To_Child_Frame(bones(BONE_R_TOE_3)->Frame().Inverse()*bones(BONE_R_ANKLE)->Frame()*toe_3_frame_in_parent);}

    // Ring toe
    if(Add_Joint_Safe2(BONE_R_ANKLE,BONE_R_TOE_4,JOINT_R_TOE_4,new ANGLE_JOINT<TV>())){
        TV toe_4_joint_center((T)-0.00912639,(T)-0.0923052,(T)-0.0314254);
        TV toe_4_axis_of_rotation((T)(-0.0138695+0.00596432),(T)(-0.0907241+0.0923052),(T)(-0.0345874+0.0282633));
        FRAME<TV> toe_4_frame_in_parent(toe_4_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),toe_4_axis_of_rotation));
        joint(JOINT_R_TOE_4)->Set_Joint_To_Parent_Frame(toe_4_frame_in_parent);
        joint(JOINT_R_TOE_4)->Set_Joint_To_Child_Frame(bones(BONE_R_TOE_4)->Frame().Inverse()*bones(BONE_R_ANKLE)->Frame()*toe_4_frame_in_parent);}

    // Pinky toe
    if(Add_Joint_Safe2(BONE_R_ANKLE,BONE_R_TOE_5,JOINT_R_TOE_5,new ANGLE_JOINT<TV>())){
        TV toe_5_joint_center((T)-0.0249367,(T)-0.0812379,(T)-0.0393305);
        TV toe_5_axis_of_rotation((T)(-0.0296798+0.0186126),(T)(-0.0796569+0.0812379),(T)(-0.0424926+0.0377495));
        FRAME<TV> toe_5_frame_in_parent(toe_5_joint_center,ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),toe_5_axis_of_rotation));
        joint(JOINT_R_TOE_5)->Set_Joint_To_Parent_Frame(toe_5_frame_in_parent);
        joint(JOINT_R_TOE_5)->Set_Joint_To_Child_Frame(bones(BONE_R_TOE_5)->Frame().Inverse()*bones(BONE_R_ANKLE)->Frame()*toe_5_frame_in_parent);}

    for(int i=0;i<num_joints;i++) if(joint(i)) joint(i)->Set_Name(joint_names[i-1]);

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(data_directory+"/SIMM_Data/hand_joints",false);
    std::string next;
    std::string joint_name,joint_type,bone_child,bone_parent,joint_center_frame;
    int num_joints_to_read=0;
    TV joint_center,axis;
    FRAME<TV> joint_frame;
    *input>>num_joints_to_read;
    if(verbose) LOG::cout<<"Found "<<num_joints_to_read<<" to read\n";
    *input>>next;
    int current_joint_index=JOINT_R_5PROXPH;
    for(int j=0;j<num_joints_to_read;j++){
        if(next=="beginjoint"){
            JOINT<TV>* new_joint=0;int bone_parent_index=0,bone_child_index=0;
            *input>>next;
            while(next!="endjoint"){
                if(next=="joint_name"){*input>>joint_name;if(verbose) LOG::cout<<"Processing joint: "<<joint_name<<std::endl;}
                else if(next=="joint_type"){
                    *input>>joint_type;
                    if(use_only_point_joints || joint_type=="point") new_joint=new POINT_JOINT<TV>();
                    else if(joint_type=="bend") new_joint=new ANGLE_JOINT<TV>();
                    else std::cerr<<"Unrecognized joint_type: "<<joint_type<<std::endl;}
                else if(next=="bone_child"){
                    *input>>bone_child;
                    bone_child_index=Get_Bone_Index(bone_child);
                    if(!bone_child_index) std::cerr<<"Unrecognized bone_child: "<<bone_child<<std::endl;}
                else if(next=="bone_parent"){
                    *input>>bone_parent;
                    bone_parent_index=Get_Bone_Index(bone_parent);
                    if(!bone_parent_index) std::cerr<<"Unrecognized bone_parent: "<<bone_parent<<std::endl;}
                else if(next=="joint_center"){*input>>joint_center;if(verbose) LOG::cout<<"Joint_center: "<<joint_center<<std::endl;}
                else if(next=="joint_frame") *input>>joint_center_frame;
                else if(next=="axis") *input>>axis;
                else std::cerr<<"Unrecognized directive: "<<next<<std::endl;
                *input>>next;}

            if(Add_Joint_Safe(bone_parent_index,bone_child_index,current_joint_index++,new_joint)){
                RIGID_BODY<TV> *parent=bones(bone_parent_index),*child=bones(bone_child_index);
                joint_frame=FRAME<TV>(joint_center);
                if(joint_type=="bend" || use_only_point_joints) joint_frame.r=ROTATION<TV>::From_Rotated_Vector(TV(1,0,0),axis);

                if(joint_center_frame=="parent"){
                    new_joint->Set_Joint_To_Parent_Frame(joint_frame);
                    new_joint->Set_Joint_To_Child_Frame(child->Frame().Inverse()*parent->Frame()*joint_frame);}
                else if(joint_center_frame=="child"){
                    new_joint->Set_Joint_To_Parent_Frame(parent->Frame().Inverse()*child->Frame()*joint_frame);
                    new_joint->Set_Joint_To_Child_Frame(joint_frame);}
                else if(joint_center_frame=="world"){
                    new_joint->Set_Joint_To_Parent_Frame(parent->Frame().Inverse()*joint_frame);
                    new_joint->Set_Joint_To_Child_Frame(child->Frame().Inverse()*joint_frame);}
                else std::cerr<<"Unrecognized joint_center_frame: "<<joint_center_frame<<std::endl;

                new_joint->Set_Name(joint_name);}}
        *input>>next;}
    delete input;

    if(do_left) for(int i=medial_joints+Individual_Side_Joints()+1;i<JOINT_END;i++) Add_Reflected_Left_Joint(i);
}
//#####################################################################
// Function Get_Bone_Index
//#####################################################################
int Get_Bone_Index(const std::string& bone_name)
{
    int index=0;
    if(bone_names.Find(bone_name,index)) return index;
    std::string search_name;
    if(bone_name=="foot") search_name="foot_no_toes_right";
    else if(bone_name=="thorax") search_name="sternum_ribs_back";
    else if(bone_name=="pelvis") search_name="hip";
    else search_name=bone_name+"_right";
    if(bone_names.Find(search_name,index)) return index;
    return 0;
}
//#####################################################################
// Function Reflected_Constrained_Point
//#####################################################################
ATTACHMENT_POINT<TV>* Reflected_Constrained_Point(ATTACHMENT_POINT<TV>* constrained_point)
{
    int local_bone_index=0;
    int reflected_bone_index=Reflected_Bone(local_bone_index);
    RIGID_BODY<TV>* reflected_body=bones(reflected_bone_index);if(!reflected_body) return 0; // bone might not exist if we're filtering some out
    TV reflected_point_in_object_space=Reflected_Frame_Relative_To_Body(reflected_bone_index,FRAME<TV>(constrained_point->object_space_position)).t;
    return new ATTACHMENT_POINT<TV>(*reflected_body,reflected_point_in_object_space);
}
//#####################################################################
// Function Reflected_Muscle
//#####################################################################
// MUSCLE<T,GRID<TV> >* Reflected_Muscle(MUSCLE<T,GRID<TV> >* muscle)
MUSCLE<TV>* Reflected_Muscle(MUSCLE<TV>* muscle)
{
    // MUSCLE<T,GRID<TV> >* reflected_muscle=new MUSCLE<T,GRID<TV> >(muscle->force_curve);

    MUSCLE<TV>* reflected_muscle=new MUSCLE<TV>(muscle->force_curve);
    reflected_muscle->Set_Name(muscle->name+"_left");
    reflected_muscle->Set_Attachment_Point_1(Reflected_Constrained_Point(muscle->attachment_point_1));
    reflected_muscle->Set_Attachment_Point_2(Reflected_Constrained_Point(muscle->attachment_point_2));
    for(int i=1;i<=muscle->via_points.m;i++) reflected_muscle->Add_Via_Point(Reflected_Constrained_Point(muscle->via_points(i)));
    // if any of the attachment/via points returned null then we can't create the muscle (some body doesn't exist)
    int dummy_index;if(!reflected_muscle->attachment_point_1 || !reflected_muscle->attachment_point_2 || reflected_muscle->via_points.Find(0,dummy_index)){delete reflected_muscle;return 0;}
    if(verbose) LOG::cout<<"Reflecting "<<muscle->name<<" -> "<<reflected_muscle->name<<std::endl;
    if(muscle->optimal_length) reflected_muscle->Set_Optimal_Length(muscle->optimal_length);
    if(muscle->peak_force) reflected_muscle->Set_Peak_Force(muscle->peak_force);
    if(muscle->pennation_angle) reflected_muscle->Set_Pennation_Angle(muscle->pennation_angle);
    if(muscle->tendon_slack_length) reflected_muscle->Set_Tendon_Slack_Length(muscle->tendon_slack_length);
    if(muscle->max_shortening_velocity) reflected_muscle->Set_Max_Shortening_Velocity(muscle->max_shortening_velocity);
    return reflected_muscle;
}
//#####################################################################
// Function Make_Muscles
//#####################################################################
void Make_Muscles()
{
    std::string viapoint_file_directory=data_directory+"/SIMM_Data/muscle_viapoints/";
    std::string muscle_params_directory=data_directory+"/SIMM_Data/muscle_params/";
    std::string muscle_files[]={"add_brev","ECU","FDSM","LAT2","sar",
                                "add_long","EDCI","FDSR","LAT3","semimem",
                                "add_mag1","EDCL","flex_dig","lat_gas","semiten",
                                "add_mag2","EDCM","flex_hal","med_gas","soleus",
                                "add_mag3","EDCR","FPL","pat_lig","solstiff",
                                "ANC","EDM","gem","PECM1","SUBSC",
                                "APL","EIP","glut_max1","PECM2","SUPSP",
                                "BIClong","EPB","glut_max2","PECM3","SUP",
                                "BICshort","EPL","glut_max3","pect","tfl",
                                "bifemlh","ext_dig","glut_med1","per_brev","tib_ant",
                                "bifemsh","ext_hal","glut_med2","peri","tib_post",
                                "BRA","FCR","glut_med3","per_long","TMAJ",
                                "BRD","FCU","glut_min1","per_tert","TMIN",
                                "CORB","FDPI","glut_min2","PL","TRIlat",
                                "DELT1","FDPL","glut_min3","PQ","TRIlong",
                                "DELT2","FDPM","grac","psoas","TRImed",
                                "DELT3","FDPR","iliacus","PT","vas_int",
                                "ECRB","FDSI","INFSP","quad_fem","vas_lat",
                                "ECRL","FDSL","LAT1","rect_fem","vas_med",
                                "TRAP1","TRAP2","sternohyoid","omohyoid","PLAT1",
                                "PLAT2","PLAT3","SPLCAP1","SPLCAP2","supoblcap","infoblcap",
                                "rectcap","rectcapmaj","splcervicis","semisplcap","semisplcervicis",
                                "longcap","serisplcap","BKOBL1","BKOBL2","BKOBL3","FROBL1","FROBL2",
                                "transverse","lowerback"};
    int num_muscles=120;
    if(verbose) LOG::cout<<"Making muscles: ";

    // parameters for now
    ARRAY<T> parameters;
    parameters.Resize(5);for(int p=0;p<5;p++) parameters(p)=0;
    for(int m=0;m<num_muscles;m++){
        ARRAY<T_MUSCLE_SEGMENT_DATA>* segment_data=new ARRAY<T_MUSCLE_SEGMENT_DATA>();
        MUSCLE<TV>* muscle=new MUSCLE<TV>(arb->muscle_list->muscle_force_curve);
        muscle->Set_Name(muscle_files[m-1]);
        std::istream* input=FILE_UTILITIES::Safe_Open_Input(viapoint_file_directory+muscle_files[m-1]+".viapts",false);

        if(verbose) LOG::cout<<muscle_files[m-1]<<" ";
        int num_points;float x,y,z;std::string bone_name,attached,st;
        *input>>num_points;
        bool success=true;
        for(int i=0;i<num_points;i++){
            *input>>x>>y>>z>>bone_name>>attached;
            if (i!=1){
                *input>>st;
                if(st=="linear")
                    segment_data->Append(T_MUSCLE_SEGMENT_DATA(MUSCLE_SEGMENT<TV>::LINEAR_SEGMENT,ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_NONE,parameters));
                if(st=="analytic"){
                    *input>>parameters(1)>>parameters(2)>>parameters(4)>>parameters(5);
                    segment_data->Append(T_MUSCLE_SEGMENT_DATA(MUSCLE_SEGMENT<TV>::ANALYTIC_SURFACE_SEGMENT,ANALYTIC_SURFACE_MUSCLE_SEGMENT<T>::CURVE_COSINE,parameters));}}
            std::cout<<parameters(1)<<" "<<parameters(2)<<" "<<parameters(3)<<" "<<parameters(4)<<" "<<parameters(5)<<std::endl;
            std::cout<<bone_name<<std::endl;
            RIGID_BODY<TV>* rigid_body=bones(Get_Bone_Index(bone_name));
            if(rigid_body){ // body might not exist if we're filtering part of the skeleton
                ATTACHMENT_POINT<TV>* attachment_point=new ATTACHMENT_POINT<TV>(*rigid_body,TV(x,y,z));
                if(attached=="attach"){//attchment point
                    if(!muscle->attachment_point_1) muscle->attachment_point_1=attachment_point;
                    else muscle->attachment_point_2=attachment_point;}
                else muscle->Add_Via_Point(attachment_point);

//                if(previous_constrained_point) muscle->muscle_segments.Append(muscle->Create_Muscle_Segment(previous_constrained_point,attachment_point, segment_type,lots of parameteras));
            }
            else success=false;}
        delete input;
        if((!muscle->attachment_point_1 || !muscle->attachment_point_2) && !bones(BONE_HIP)) success=false;

        if(success){
            if(!muscle->attachment_point_1) muscle->attachment_point_1=new ATTACHMENT_POINT<TV>(*bones(BONE_HIP),TV());
            if(!muscle->attachment_point_2) muscle->attachment_point_2=new ATTACHMENT_POINT<TV>(*bones(BONE_HIP),TV());

            PARAMETER_LIST muscle_param;muscle_param.Read(muscle_params_directory+muscle_files[m-1]+".param");
            if(muscle_param.Is_Defined("optimal_fiber_length")) muscle->Set_Optimal_Length(muscle_param.Get_Parameter("optimal_fiber_length",(T)0));
            if(muscle_param.Is_Defined("peak_force")) muscle->Set_Peak_Force(muscle_param.Get_Parameter("peak_force",(T)0));
            if(muscle_param.Is_Defined("pennation_angle")) muscle->Set_Pennation_Angle(muscle_param.Get_Parameter("pennation_angle",(T)0));
            if(muscle_param.Is_Defined("tendon_slack_length")) muscle->Set_Tendon_Slack_Length(muscle_param.Get_Parameter("tendon_slack_length",(T)0));
            muscles.Append(muscle);segment_data_map.Set(muscle->name,segment_data);
            MUSCLE<TV>* reflected_muscle=Reflected_Muscle(muscle);
            if(reflected_muscle){muscles.Append(reflected_muscle);segment_data_map.Set(reflected_muscle->name,segment_data);}}
        else{delete muscle;muscle=0;}}
}
//#####################################################################
// Function Initialize Muscle_Segments
//#####################################################################
void Initialize_Muscle_Segments()
{
    for(int i=1;i<=arb->muscle_list->muscles.m;i++){
        ARRAY<T_MUSCLE_SEGMENT_DATA> *curr_segment_data=0;segment_data_map.Get(arb->muscle_list->muscles(i)->name,curr_segment_data);
        arb->muscle_list->muscles(i)->Initialize(*curr_segment_data);}
}
//#####################################################################
// Function Replace_Bones_With_Fused_Bone
//#####################################################################
void Replace_Bones_With_Fused_Bone(const std::string& merge_filename,RIGID_BODY<TV>* new_bone,RIGID_BODY<TV>* new_reflected_bone=0,int specific_side=0)
{
    std::ifstream input(merge_filename.c_str());
    if(!input){LOG::cerr<<"Can't open "<<merge_filename<<std::endl;exit(1);}
    int n;
    input>>n;
    ARRAY<RIGID_BODY<TV>*> replaced_bones(n);
    ARRAY<FRAME<TV> > original_to_fused_world_transform(n);
    FRAME<TV> subtract_frame,align_frame;
    for(int i=0;i<n;i++){
        std::string bone_name;T scale;
        input>>bone_name>>scale>>original_to_fused_world_transform(i);
        if(verbose) LOG::cout<<"BONE "<<bone_name<<" got frame "<<original_to_fused_world_transform(i)<<std::endl;
        bone_name=FILE_UTILITIES::Get_Basename(FILE_UTILITIES::Get_Short_Name(bone_name));
        int index=Get_Bone_Index(bone_name);
        if(!index) continue;
        if(index) replaced_bones(i)=bones(index);
        if(i==1){assert(replaced_bones(i));subtract_frame=original_to_fused_world_transform(i);align_frame=replaced_bones(i)->Frame();}
        original_to_fused_world_transform(i)=align_frame*subtract_frame.Inverse()*original_to_fused_world_transform(i)*bones(index)->Frame().Inverse();}
    Replace_Bones_With_Fused_Bone(replaced_bones,original_to_fused_world_transform,new_bone);

    if(new_reflected_bone){
        for(int i=0;i<n;i++){int index=0;
            bones.Find(replaced_bones(i),index);
            if(!index || !bones(Reflected_Bone(index))){LOG::cerr<<"could not replace fused reflected bone"<<std::endl;return;}
            replaced_bones(i)=bones(Reflected_Bone(index));
            MATRIX<T,4> reflection=transform.Matrix_4X4()*Reflection_Matrix(0,specific_side)*transform.Inverse().Matrix_4X4();
            MATRIX<T,4> reflected_transform=reflection.Inverse()*original_to_fused_world_transform(i).Matrix_4X4()*reflection;
            original_to_fused_world_transform(i)=FRAME<TV>(TV(),ROTATION<TV>(reflected_transform.Upper_3X3()));}
        Replace_Bones_With_Fused_Bone(replaced_bones,original_to_fused_world_transform,new_reflected_bone);}
}
//#####################################################################
// Function Replace_Bones_With_Fused_Bone
//#####################################################################
void Replace_Bones_With_Fused_Bone(const ARRAY<RIGID_BODY<TV>*>& replaced_bones,const ARRAY<FRAME<TV> >& original_to_fused_world_transform,RIGID_BODY<TV>* new_bone)
{
    for(int i=0;i<muscles.m;i++){
        bool add_muscle=false;int index=0;
        for(int j=0;j<=muscles(i)->via_points.m+1;j++){
            ATTACHMENT_POINT<TV>* constrained_point=(j==0)?muscles(i)->attachment_point_1:((j<=muscles(i)->via_points.m)?muscles(i)->via_points(j):muscles(i)->attachment_point_2);

            for(int b=0;b<replaced_bones.m;b++) if(replaced_bones(b)->particle_index==constrained_point->particle_index) index=b;
            if(constrained_point && index){
                TV old_object_space_position=(constrained_point)->object_space_position;
                TV new_object_space_position=new_bone->Frame().Inverse()*original_to_fused_world_transform(index)*replaced_bones(index)->Frame()*old_object_space_position;
                if(verbose)
                    LOG::cout<<"Moving attachment point of muscle "<<muscles(i)->name<<" from body "<<replaced_bones(index)->name<<" ("<<
                    old_object_space_position<<") to fused body "<<new_bone->name<<" ("<<new_object_space_position<<")"<<std::endl;
                delete constrained_point;
                constrained_point=new ATTACHMENT_POINT<TV>(arb->muscle_list->particles,0,new_bone->rigid_body_collection,new_bone->particle_index,new_object_space_position);
                add_muscle=true;
                if(j==0) muscles(i)->attachment_point_1=constrained_point;
                else if(j<=muscles(i)->via_points.m) muscles(i)->via_points(j)=constrained_point;
                else muscles(i)->attachment_point_2=constrained_point;}}
        if(add_muscle){
            assert(!arb->muscle_list->muscles.Find(muscles(i),index));
            if(verbose) LOG::cout<<"Adding muscle "<<muscles(i)->name<<" attached to fused bone"<<std::endl;
            arb->muscle_list->Add_Muscle(muscles(i));}}
}
//#####################################################################
// Function Make Muscles -- makes lots and lots of muscles
//#####################################################################
void Make_Wrapping_Objects()
{
    // Okay, this isn't going to do what we want it to, but we might use it later if we add in the wrapping surfaces so I'm leaving it here
    std::string next,object_name,wrap_type,bone_name;
    float x,y,z,radius,radiust,height,scale_factor;
    int id;
    TV translation,rotation,radius_vector;
    RIGID_BODY<TV>* rigid_body=0;
    MATRIX<T,4> bone_rotation_matrix;
    TV bone_translation;
    TV bone_scale;

    std::istream* input=FILE_UTILITIES::Safe_Open_Input(upper_body_file_jnt,false);
    while(!input->eof()){
        *input>>next;
        if(next=="beginwrapobject"){
            while(next!="endwrapobject"){
                *input>>next;
                if(next=="wraptype") *input>>wrap_type;
                else if(next=="segment") *input>>bone_name;
                else if(next=="visible") *input>>next;
                else if(next=="xyz_body_rotation"){*input>>x>>y>>z;rotation=TV(x,y,z);}
                else if(next=="translation"){*input>>x>>y>>z;translation=TV(x,y,z);}
                else if(next=="radius"){
                    if(wrap_type=="ellipsoid"){*input>>x>>y>>z;radius_vector=TV(x,y,z);}
                    else if(wrap_type=="sphere") *input>>radius;
                    else if(wrap_type=="cylinder") *input>>radius>>next>>height;
                    else if(wrap_type=="torus") *input>>radius>>radiust;
                    else std::cerr<<"Unrecognized wrap type "<<wrap_type<<std::endl;}
                else if(next=="quadrant") *input>>next; // don't know what to do with this information
                else std::cerr<<"Unrecognized token: "<<next<<std::endl;}}
            // now we make the wrap object and insert it into the rigid bodies list
            if(wrap_type=="ellipsoid"){}
            else if(wrap_type=="sphere"){}
            else if(wrap_type=="cylinder"){}
            else if(wrap_type=="torus"){}
            else std::cerr<<"Unrecognized wrap type "<<wrap_type<<std::endl;

            rigid_body->frame=FRAME<TV>(translation,ROTATION<TV>::From_Euler_Angles(rotation.x,rotation.y,rotation.z));

            // now get the transform for the bone that it's attached to
            std::istream* input_bone;
            if(bone_name=="thorax") input_bone=FILE_UTILITIES::Safe_Open_Input(simm_directory+"PhysBAM_to_SIMM_transforms/"+bone_name+".transform",false);
            else input_bone=FILE_UTILITIES::Safe_Open_Input(simm_directory+"PhysBAM_to_SIMM_transforms/"+bone_name+"_right.transform",false);
            input_bone>>bone_rotation_matrix>>bone_translation>>bone_scale;
            //mouse_handler->Update_Particles();

            //id=solids_parameters.rigid_body_parameters.list.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/"+bone_files[i-1],scale_factor,true,false,false);
            //apply bone transform

            rigid_body=bones(id);
            rigid_body->Set_Coefficient_Of_Restitution((T).5);
            rigid_body->Set_Coefficient_Of_Friction((T).5);
            rigid_body->Set_Name(object_name);

/* code for updating hte points of the triangulated surface . . may need to use to update the via_points on the wrapping surfaces
        for(int i=1;i<=triangulated_surface.particles.array_collection->Size();i++)
            triangulated_surface.particles.X(i)=Current_World_Space_Position(original_particle_positions(i));
        for(int j=0;j<points.m;j++)
        points(j)=Current_World_Space_Position(original_points_positions(j));*/
    }
    delete input;
}
//#####################################################################
// Function Read_Joint_Limits
//#####################################################################
static void Read_Joint_Limits(ARTICULATED_RIGID_BODY<TV>* arb,PARAMETER_LIST& parameter_list)
{
    for(int i=1;i<=arb->joint_mesh.joints.m;i++){JOINT<TV>* joint=arb->joint_mesh.joints(i);
        std::string prefix=joint->name+".limits.";
        VECTOR<T,2> x,y,z;
        JOINT<TV>* mirrored_joint=0;
        for(int j=1;j<=arb->joint_mesh.joints.m;j++) if(arb->joint_mesh.joints(j)->name==joint->name+"_left") mirrored_joint=arb->joint_mesh.joints(j);
        if(joint->joint_type==JOINT<TV>::TYPE_POINT_JOINT){
            if(parameter_list.Is_Defined(prefix+"x")){
                VECTOR<T,2> limits=(T)(pi/180)*parameter_list.Get_Parameter(prefix+"x",VECTOR<T,2>());
                ((POINT_JOINT<TV>*)joint)->Use_Twist_Constraint(limits.x,limits.y);
                if(mirrored_joint) ((POINT_JOINT<TV>*)mirrored_joint)->Use_Twist_Constraint(limits.x,limits.y);}
            if(parameter_list.Is_Defined(prefix+"y")){
                VECTOR<T,2> limits=(T)(pi/180)*parameter_list.Get_Parameter(prefix+"y",VECTOR<T,2>());
                ((POINT_JOINT<TV>*)joint)->Use_Phi_Constraint(limits.x,limits.y);
                if(mirrored_joint) ((POINT_JOINT<TV>*)mirrored_joint)->Use_Phi_Constraint(-limits.y,-limits.x);}
            if(parameter_list.Is_Defined(prefix+"z")){
                VECTOR<T,2> limits=(T)(pi/180)*parameter_list.Get_Parameter(prefix+"z",VECTOR<T,2>());
                ((POINT_JOINT<TV>*)joint)->Use_Theta_Constraint(limits.x,limits.y);
                if(mirrored_joint) ((POINT_JOINT<TV>*)mirrored_joint)->Use_Theta_Constraint(-limits.y,-limits.x);}}
        else if(joint->joint_type==JOINT<TV>::TYPE_ANGLE_JOINT){
            if(parameter_list.Is_Defined(prefix+"x")){
                VECTOR<T,2> limits=(T)(pi/180)*parameter_list.Get_Parameter(prefix+"x",VECTOR<T,2>());
                ((ANGLE_JOINT<TV>*)joint)->Set_Angle_Constraints(true,limits.x,limits.y);
                if(mirrored_joint) ((ANGLE_JOINT<TV>*)mirrored_joint)->Set_Angle_Constraints(true,limits.x,limits.y);}}}
}
//#####################################################################
// Function Load_Splines
//#####################################################################
// - call load splines
// - then call Read_Frame_Track_From_Spline and set the result to the joint's joint_function's track
static void Load_Splines(const std::string& filename,ARRAY<PAIR<std::string,BSPLINE_QUATERNION<T>* > >& motion_splines,const int spline_order,const bool quaternion_check=true)
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(filename,false);
    std::cout<<"Loading from "<<filename<<std::endl;

    ARRAY<T> keyframe_times;
    int number_of_keyframes,number_of_joints;
    *input>>number_of_keyframes;
    keyframe_times.Resize(number_of_keyframes);
    for(int i=0;i<keyframe_times.m;i++) (*input)>>keyframe_times(i);;

    *input>>number_of_joints;
    motion_splines.Resize(number_of_joints);
    for(int i=0;i<motion_splines.m;i++){
        (*input)>>motion_splines(i).x;
        ARRAY<ROTATION<TV> > control_points;control_points.Resize(number_of_keyframes);
        for(int j=0;j<control_points.m;j++){FRAME<TV> frame;(*input)>>frame;control_points(j)=frame.r;}
        motion_splines(i).y=new BSPLINE_QUATERNION<T>(keyframe_times,control_points,spline_order);
        if(quaternion_check) motion_splines(i).y->Quaternion_Check();
        motion_splines(i).y->Create_Closed_Points();}
    delete input;
}
//#####################################################################
// Function Print_Splines
//#####################################################################
static void Print_Splines(ARRAY<PAIR<std::string,BSPLINE_QUATERNION<T>* > >*& motion_splines)
{
    for(int i=0;i<motion_splines.m;i++){
        std::cout<<"Joint: "<<motion_splines(i).x<<std::endl;
        motion_splines(i).y->Print_Control_Points_And_Times();}
}
//#####################################################################
// Function Read_Frame_Track_From_Spline
//#####################################################################
static INTERPOLATION_CURVE<T,ROTATION<TV> >* Read_Frame_Track_From_Spline(ARRAY<PAIR<std::string,BSPLINE_QUATERNION<T>* > >& motion_splines,const std::string& jointname,const bool periodic=true,
    const int samples=1000,const T track_speedup_factor=1)
{
    PHYSBAM_FATAL_ERROR();
#if 0
    T tmin=0;int index=0;
    for(int i=0;i<motion_splines.m;i++) if(motion_splines(i).x==jointname) index=i;
    if(!index) return 0;
    FRAME_TRACK_3D<T>* frame_track=new FRAME_TRACK_3D<T>(samples,tmin,track_speedup_factor);
    T time_increment=(1/(T)(samples-1))*motion_splines(index).y->Range();
    for(int i=0;i<samples-1;i++) frame_track->trajectory(i+1)=FRAME<TV>(motion_splines(index).y->Evaluate(motion_splines(index).y->Start_Time()+i*time_increment));
    frame_track->trajectory(samples)=FRAME<TV>(motion_splines(index).y->Evaluate(motion_splines(index).y->Start_Time()));
    frame_track->periodic=periodic;
    return frame_track;
#endif
    return 0;
}
//#####################################################################
};
}
#endif
