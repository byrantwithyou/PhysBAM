//#####################################################################
// Copyright 2004-2006, Kevin Der, Ranjitha Kumar, Mike Rodgers, Rachel Weinstein.
// This file is part of PhysBAM whose distribution is governed by the license 
// contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class VISIBLE_HUMAN_EXAMPLE
//##################################################################### 
#ifndef __VISIBLE_HUMAN_EXAMPLE__
#define __VISIBLE_HUMAN_EXAMPLE__

#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/LEVELSET_MAKER_UNIFORM.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TETRAHEDRON_COLLISION_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Dynamics/Motion/MOTION_SEQUENCE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include "../ARB_PARAMETERS.h"
#include "../VISIBLE_HUMAN.h"

namespace PhysBAM{

template<class T>
class VISIBLE_HUMAN_EXAMPLE:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> >
{
    typedef VECTOR<T,3> TV;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
public:
    using BASE::first_frame;using BASE::last_frame;using BASE::frame_rate;using BASE::stream_type;
    using BASE::restart;using BASE::restart_frame;using BASE::output_directory;
    using BASE::solids_parameters;using BASE::write_last_frame;using BASE::data_directory;using BASE::fluids_parameters;
    
    enum MESH_TYPE{FULL_BODY=1,TORSO_AND_ARMS};
    static const MESH_TYPE mesh_file=TORSO_AND_ARMS;
    static const bool USE_TESTING_TOLERANECES=true;

    ARTICULATED_RIGID_BODY<TV>* arb;
    PARAMETER_LIST parameter_list;
    FRAME_TRACK_3D<T> *shoulder_track,*elbow_track; 
    VISIBLE_HUMAN<T>* da_man;
    FRAME<TV> skeleton_frame;
    ARRAY<RIGID_BODY_MASS<TV> > saved_masses;

    enum GEOMETRY_TYPE {TETRAHEDRALIZED_VOLUME_TYPE,TRIANGULATED_SURFACE_TYPE,EMBEDDED_TETRAHEDRALIZED_VOLUME_TYPE};
    int number_of_deformable_bodies;
    ARRAY<std::string> deformable_body_geometry_filenames;
    ARRAY<GEOMETRY_TYPE> deformable_body_geometry_types;
    ARRAY<RIGID_BODY_STATE<TV>* > deformable_body_initial_states;
    ARRAY<ARRAY<TV> > deformable_body_rest_positions;
    
    T tracks_start_time;
    int root,first_wrist_joint_id;
    bool add_ground,add_block,use_gravity,draw_flesh,simulate;
    bool tracks_added,use_tracks,test_frame_tracks;

    // Data pertaining to enslaving flesh to bone
    int num_bones_present;
    ARRAY<int> bone_ids;
    ARRAY<ARRAY<int> > enslaved_nodes; // indicies of enslaved nodes. 
    ARRAY<ARRAY<TV> > positions_relative_to_bone_frames; // corresponding positions
    
    // Callbacks
    //void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}

    VISIBLE_HUMAN_EXAMPLE(const STREAM_TYPE stream_type,const int simulate_input,std::string parameter_file="")
        :BASE(stream_type,0,FLUIDS_PARAMETERS_UNIFORM<GRID<TV> >::NONE),
        number_of_deformable_bodies(0),root(1),simulate(simulate_input>0),tracks_added(false)
    {
        fluids_parameters.simulate=false;
        std::cout<<"simulate: "<<simulate<<std::endl;
        //restart=true;
        //restart_frame=3;
        last_frame=600;
        frame_rate=24;
        write_last_frame=true;
        output_directory="Visible_Human/output";

        //#####################################################################
        // SET SOLID PARAMETERS
        solids_parameters.cfl=(T).5;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
        solids_parameters.implicit_solve_parameters.cg_iterations=900;
        solids_parameters.newton_iterations=10;
        solids_parameters.newton_tolerance=(T)1e-2;
        solids_parameters.use_partially_converged_result=true;
        solids_parameters.collision_tolerance=(T)1e-4;
        solids_parameters.perform_self_collision=false;
        solids_parameters.perform_collision_body_collisions=false;
        solids_parameters.gravity=(T)9.8;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=simulate;;
    }

    ~VISIBLE_HUMAN_EXAMPLE()
    {
        delete arb;
    }

//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    //#####################################################################
    // CREATE AND INITIALIZE ARTICULATED_RIGID_BODY
    arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb->Set_Iterative_Tolerance((T)1e-2/*8*/);
    arb->Set_Line_Search_Interval_Tolerance((T)1e-2/*8*/);
    arb->Set_Contact_Level_Iterations(5);
    arb->Set_Max_Line_Search_Iterations(20/*1000*/);
    arb->Set_Max_Iterations(/*1000*/20);
    //arb->Set_Shock_Propagation_Level_Iterations(100);
    //arb->Set_Poststabilization_Iterations(100);
    //arb->Set_Use_Shock_Propagation(false);
    arb->Set_Do_Final_Pass(false);

    //#####################################################################
    // LOAD PARAMETERS
    if(parameter_file.empty()) parameter_file="Visible_Human/example.param";
    ARB_PARAMETERS::Read_Common_Parameters(parameter_file,*this,parameter_list);
    add_ground=parameter_list.Get_Parameter("add_ground",true); 
    add_block=parameter_list.Get_Parameter("add_block",false);
    use_gravity=parameter_list.Get_Parameter("use_gravity",true);
    draw_flesh=parameter_list.Get_Parameter("draw_flesh",true);
    use_tracks=parameter_list.Get_Parameter("use_tracks",false);
    test_frame_tracks=parameter_list.Get_Parameter("test_frame_tracks",false);
    tracks_start_time=parameter_list.Get_Parameter("tracks_start_time",(T)0);
    if(test_frame_tracks) simulate=false;

    //#####################################################################
    // Turn down precision?
    if(USE_TESTING_TOLERANECES){
        solids_parameters.implicit_solve_parameters.cg_iterations=10;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.collision_tolerance=(T)1e-2;
        arb->Set_Iterative_Tolerance((T)1e-2);}

    //#####################################################################
    // MUSCLES AND TRACKS
    if(simulate) arb->Use_Muscle_Actuators();   
    if(use_tracks) Create_Push_Up_Tracks();
    //arb->Use_PD_Actuators();

    // Initialize deformable bodies
    Get_Initial_Data();
    if(draw_flesh){
        DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
        DEFORMABLE_PARTICLES<TV>& particles=deformable_object.particles;    
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        if(use_gravity) solid_body_collection.Add_Force(new GRAVITY<TV>(particles));
        solid_body_collection.Add_Force(Create_Quasistatic_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,true));
        solid_body_collection.Update_Fragments();}

    std::cout<<"done initializing example\n";

    if(!restart && tracks_start_time==0) Add_Tracks(0);
    if(!restart && parameter_list.Get_Parameter("update_joints",true)) Update_Joints(0);
    if(draw_flesh) Get_Constrained_Particle_Data();
    
    std::cout<<"SOLIDS_PARAMTER "<< solids_parameters.cfl<<std::endl;
}
//#####################################################################
// Function Initialize_Rigd_Bodies
//#####################################################################
void Initialize_Rigid_Bodies()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    // Initialize muscle list for future muscle creation
    arb->muscle_list->muscle_force_curve.Initialize(data_directory);
    //Default_Example2();
    Push_Up();
    if(add_block){
        int id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/box");
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->Frame().t=TV(0,1,(T)1.4);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("block");
        rigid_body->is_static=true;}

    if(add_ground){
        int id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/ground");
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("ground");
        rigid_body->is_static=true;
        rigid_body->add_to_spatial_partition=false;}

    if(parameter_list.Get_Parameter("add_hand_joints",false)) Add_Hand_Ground_Joints(1,2);
    if(parameter_list.Get_Parameter("add_foot_joints",false)) Add_Foot_Ground_Joints();
    if(parameter_list.Get_Parameter("static_hands",false)){
        arb->rigid_body_list.rigid_bodies(1)->is_static=true;
        arb->rigid_body_list.rigid_bodies(2)->is_static=true;}
    if(parameter_list.Get_Parameter("static_feet",false)){
        da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ANKLE)->is_static=true;
        da_man->bones(VISIBLE_HUMAN<T>::BONE_L_ANKLE)->is_static=true;}

    ARRAY<RIGID_BODY<TV>*>& rigid_bodies=rigid_body_list.rigid_bodies;
    if(use_gravity) for(int i=0;i<rigid_bodies.m;i++) if(!rigid_bodies(i)->is_static)
        rigid_bodies(i)->Add_Basic_Forces(solids_parameters.gravity,solids_parameters.gravity_direction,solids_parameters.rigid_body_evolution_parameters.rigid_body_ether_viscosity,0);
    solids_parameters.collision_body_list.Add_Bodies(rigid_body_particles);
    
    da_man->Initialize_Muscle_Segments();
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    draw_flesh = false;
    if(draw_flesh){
        TV flesh_translation(-(T).294107,(T).9355,(T).165); //determined via estimate-and-test
        flesh_translation+=TV(0,(T).05,0);
        QUATERNION<T> flesh_rotation=QUATERNION<T>((T)pi,TV(0,0,1));
        // TODO: update filename for flesh tet_vol.  (File is not currently checked into CVS.)
        Add_Deformable_Body(data_directory+"/VH_Flesh/single_body.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<TV>(FRAME<TV>(flesh_translation,flesh_rotation))); // High-res version
        // Add_Deformable_Body(data_directory+"/VH_Flesh/body_mesh_lowres.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<T,TV>(FRAME<TV>(flesh_translation,flesh_rotation))); // Low-res version
        /*TV flesh_translation; QUATERNION<T> flesh_rotation;*/
        switch(mesh_file){
          case FULL_BODY:
            flesh_translation=TV(-(T).294107,(T).9355,(T).165); // match skeleton default position (determined via estimate-and-test)
            flesh_translation+=TV(0,(T).05,0);            // skeleton translated up this amount
            flesh_rotation=QUATERNION<T>((T)pi,TV(0,0,1));
            // TODO: update filename for flesh tet_vol.  (File is not currently checked into CVS.)
            Add_Deformable_Body(data_directory+"/VH_Flesh/single_body.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<TV>(FRAME<TV>(flesh_translation,flesh_rotation))); // High-res version
            // Add_Deformable_Body(data_directory+"/VH_Flesh/body_mesh_lowres.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<T,TV>(FRAME<TV>(flesh_translation,flesh_rotation))); // Low-res version
            // Add_Deformable_Body(data_directory+"/Tetrahedralized_Volumes/fine_sphere.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<T,TV>(FRAME<TV>(flesh_translation,flesh_rotation))); // For testing only.
            break;
            
          case TORSO_AND_ARMS:
            flesh_translation=TV(-(T).294107,(T)1.3455,(T).157); // match skeleton default position (determined via estimate-and-test)
            flesh_translation+=TV(0,(T).05,0);             // skeleton translated up this amount
            flesh_rotation=QUATERNION<T>(-(T)pi/2,TV(1,0,0))*QUATERNION<T>((T)pi,TV(0,0,1));
            Add_Deformable_Body(data_directory+"/VH_Flesh/union_24.tet",TETRAHEDRALIZED_VOLUME_TYPE,new RIGID_BODY_STATE<TV>(FRAME<TV>(flesh_translation,flesh_rotation)));
            break;
        }

        // deformable bodies
        DEFORMABLE_OBJECT<TV>& deformable_object=solid_body_collection.deformable_object;
        deformable_body_rest_positions.Resize(number_of_deformable_bodies);

        for(int index=0;index<number_of_deformable_bodies;index++){
            // initialize geometry
            STRUCTURE<TV>* structure=Create_Tetrahedralized_Volume(index);

            // add to deformable_object
            deformable_object.Add_Structure(structure->Append_Particles_And_Create_Copy(deformable_object.particles));
            //Disable collisions - they aren't needed for quasistatic
            //deformable_object.collisions.collision_structures.Append(deformable_object.structures.Last());
            //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_object.structures.Last());
            delete structure;}

        // correct number nodes
        for(int i=0;i<deformable_object.structures.m;i++) deformable_object.structures(i)->Update_Number_Nodes();}

    // Rigid Bodies
    Initialize_Rigid_Bodies();

    if(draw_flesh){
        // Collisions [Disable collisions - they aren't needed for quasistatic]
        //solids_parameters.collision_body_list.Add_Bodies(rigid_body_particles);
        solids_parameters.perform_self_collision=false;}
        //Initialize_Tetrahedron_Collisions(1,solid_body_collection.deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>());}
}
//#####################################################################
// Write_Output_Files
//#####################################################################
/*void Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);

    std::ostream* output=FILE_UTILITIES::Safe_Open_Output("handtransform_"+STRING_UTILITIES::string_sprintf("%d",frame),false);
    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++){
        (*output)<<data_directory<<"/Rigid_Bodies/New_Visible_Human_Bones/";
        (*output)<<rigid_body_particles.Rigid_Body(i).name<<"\n"<<rigid_body_particles.Rigid_Body(i).frame<<std::endl;
    }
    delete output;
}*/
//#####################################################################
// Function Humerus_Ulna_Radius
//#####################################################################
void Humerus_Ulna_Radius()
{
    skeleton_frame=FRAME<TV>(TV(),QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T)0.5*(T)pi,TV(1,0,0)));
    da_man=new VISIBLE_HUMAN<T>(arb,data_directory,skeleton_frame);
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Humerus_And_Ulna_Filter(parameter_list.Get_Parameter("with_radius",true),false));
    da_man->bones(VISIBLE_HUMAN<T>::BONE_R_HUMERUS)->is_static=true;
}
//#####################################################################
// Function Skeleton_In_Flesh
//#####################################################################
void Skeleton_In_Flesh()
{
    skeleton_frame=FRAME<TV>(TV(0,(T).05,0),QUATERNION<T>((T)pi/2,TV(1,0,0))*QUATERNION<T>((T)pi,TV(0,1,0)));
    da_man=new VISIBLE_HUMAN<T>(arb,data_directory,skeleton_frame);
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Skeleton_In_Flesh_Filter());

    if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->id_number;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->id_number;

    T k_p=parameter_list.Get_Parameter("k_p",(T)10000);

    if(parameter_list.Get_Parameter("set_static",true)){
        da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->is_static=true;
    }

    Adjust_Joints_For_Skeleton_In_Flesh();

    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;

    // Set joint frames
    /*
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)&&da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL)){
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)->Set_Joint_Frame(shoulder_track->Frame(tracks_start_time));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL)->Set_Joint_Frame(da_man->Rotated_Track(*shoulder_track)->Frame(tracks_start_time));}    
    if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR)&&da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR)){
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR)->Set_Joint_Frame(elbow_track->Frame(tracks_start_time));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR)->Set_Joint_Frame(da_man->Rotated_Track(*elbow_track)->Frame(tracks_start_time));}*/
    
    arb->Update_With_Breadth_First_Directed_Graph(root);

    // make joint functions to try to keep this pose
    for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->muscle_control=true;
        joint_function->Set_k_p(k_p);
        joint_function->Set_Target_Angle(joint_function->Angle());
    }

    for(int i=0;i<da_man->muscles.m;i++) if(da_man->muscles(i)){
        T length=da_man->muscles(i)->Total_Length();
        T total_rest_length=da_man->muscles(i)->optimal_length + da_man->muscles(i)->tendon_slack_length;
        if(total_rest_length>(T)1.5*length){
            std::cout << "Muscle " << da_man->muscles(i)->name << ": (" << da_man->muscles(i)->optimal_length << "+" << da_man->muscles(i)->tendon_slack_length << ") " << total_rest_length << " vs " << length << std::endl;}}

    // Determine which joints will get PD as opposed to muscle actuation
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

    if(arb->use_muscle_actuators){
        LOG::cout << "Joints with only muscle control"<<std::endl;
        for(int i=0;i<arb->joint_mesh.joints.m;i++) if(arb->joint_mesh.joints(i)->joint_function && arb->joint_mesh.joints(i)->joint_function->muscle_control)
            LOG::cout << "\t" << arb->joint_mesh.joints(i)->name << std::endl;}
    for(int i=0;i<rigid_body_list.Number_Of_Elements();i++) rigid_body_particles.Rigid_Body(i).Set_Coefficient_Of_Friction(1);
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


    switch(mesh_file){
      case FULL_BODY:
        // Update all Left Joints to be symmetric
        for(int i=VISIBLE_HUMAN<T>::medial_joints+VISIBLE_HUMAN<T>::Individual_Side_Joints()+1;i<VISIBLE_HUMAN<T>::JOINT_END;i++) da_man->Update_Reflected_Left_Joint(i);
        
        // Tweak individual left joints where flesh-mesh is not semetric (still requires additional refinement).
        Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),(T).035*(T)pi);
        Augment_Child_Angle_Y(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),0); 
        Augment_Child_Angle_Z(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),0*(T)pi);
        break;

      case TORSO_AND_ARMS:
        Augment_Child_Angle_Y(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL),-(T).1049*(T)pi);
        Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL),-(T).0572*(T)pi);
        Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR),(T).051*(T)pi);
        Augment_Child_Angle_Y(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR),-(T).0257*(T)pi);
        for(int i=VISIBLE_HUMAN<T>::medial_joints+VISIBLE_HUMAN<T>::Individual_Side_Joints()+1;i<VISIBLE_HUMAN<T>::JOINT_END;i++) da_man->Update_Reflected_Left_Joint(i);
        
        Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL),-(T).0033*(T)pi);
        Augment_Child_Angle_X(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR),-(T).008*(T)pi);
        Augment_Child_Angle_Y(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR),(T).009*(T)pi);
        // TODO: da_man's elbows still protrude just a sliver in back.
        break;
    }
}
//#####################################################################
// Function Default_Example
//#####################################################################
void Default_Example2()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    skeleton_frame=FRAME<TV>(QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T).5*(T)pi,TV(1,0,0)));

    // real hands versus fake hands?
    bool upper_body=parameter_list.Get_Parameter("upper_body",true);
    bool lower_body=parameter_list.Get_Parameter("lower_body",true);
    bool head_and_neck=parameter_list.Get_Parameter("head_and_neck",true);
    bool hands=parameter_list.Get_Parameter("hands",true);
    bool real_hands=parameter_list.Get_Parameter("real_hands",false);
    bool feet=parameter_list.Get_Parameter("feet",true);

    if(hands){
        int id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_right",1,true,false,false);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_right");
        rigid_body->Frame()=skeleton_frame*rigid_body->Frame();

        id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_left",1,true,false,false);
        rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_left");
        rigid_body->Frame()=skeleton_frame*rigid_body->Frame();
    }

    bool verbose=parameter_list.Get_Parameter("verbose",false);

    da_man=new VISIBLE_HUMAN<T>(stream_type,arb,solid_body_collection.deformable_object.rigid_body_particles,data_directory,skeleton_frame,verbose);
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Basic_Filter(upper_body,lower_body,head_and_neck,real_hands,feet,true));
    da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Mass());
    da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Mass());
    // what does this line do exactly?
    if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index;

    if(hands){
        Add_Wrist_Joint(true);
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_transform",
            rigid_body_list(1),rigid_body_list(2),1);}
    
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);
    std::cout<<"THE K_P IS: "<<k_p<<std::endl;

    da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->is_static=true;
    /*
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)pi/(T)4.5,TV(0,0,1)))*FRAME<TV>(QUATERNION<T>((T)pi/6,TV(0,1,0))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T)pi/(T)4.5,TV(0,0,1)))*FRAME<TV>(QUATERNION<T>(-(T)pi/6,TV(0,1,0))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T)pi/(T)2.5,TV(1,0,0))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T)pi/(T)2.5,TV(1,0,0))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)pi/6,TV(1,0,0))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)pi/6,TV(1,0,0))));*/
    arb->Update_With_Breadth_First_Directed_Graph(root);
    
    // what are joint functions, should all joints have joint functions?
    for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i); JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->Set_k_p(k_p);
        joint_function->Set_Target_Angle(joint_function->Angle());
        joint_function->muscle_control=true;
    }

    // add wrist joint functions
    if(hands && parameter_list.Get_Parameter("wrist_joint_function",true)){
        JOINT_FUNCTION<TV>* jfunc=arb->Create_Joint_Function(first_wrist_joint_id);
        jfunc->muscle_control=true;
        jfunc->Set_Target_Angle(jfunc->Angle());jfunc->Set_k_p(k_p);
        jfunc=arb->Create_Joint_Function(first_wrist_joint_id+1);
        jfunc->muscle_control=true;
        jfunc->Set_Target_Angle(jfunc->Angle());jfunc->Set_k_p(k_p);}
    
}

//#####################################################################
// Function Default_Example
//#####################################################################
void Default_Example()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    skeleton_frame=FRAME<TV>(QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T).5*(T)pi,TV(1,0,0)));
    //skeleton_frame=FRAME<TV>(TV(0,(T)0.01,-(T)0.2))*FRAME<TV>(QUATERNION<T>((T)0.3,TV(1,0,0)))*FRAME<TV>(TV(0,-(T)0.12,0),QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T)0.5*(T)pi,TV(1,0,0)));
    //da_man=new VISIBLE_HUMAN<T>(arb,data_directory,skeleton_frame);
    //da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::All_But_Right_Hand_Filter());
    //da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ANKLE)->is_static=true;
    //da_man->bones(VISIBLE_HUMAN<T>::BONE_L_ANKLE)->is_static=true;
    bool upper_body=parameter_list.Get_Parameter("upper_body",true);
    bool lower_body=parameter_list.Get_Parameter("lower_body",true);
    bool head_and_neck=parameter_list.Get_Parameter("head_and_neck",true);
    bool hands=parameter_list.Get_Parameter("hands",true);
    bool real_hands=parameter_list.Get_Parameter("real_hands",false);
    bool feet=parameter_list.Get_Parameter("feet",true);

    if(hands){
        int id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_right",1,true,false,false);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_right");
        rigid_body->frame=skeleton_frame*rigid_body->frame;

        id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_left",1,true,false,false);
        rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_left");
        rigid_body->frame=skeleton_frame*rigid_body->frame;
    }
    
    
    bool verbose=parameter_list.Get_Parameter("verbose",false);

    da_man=new VISIBLE_HUMAN<T>(arb,data_directory,skeleton_frame,verbose);
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Basic_Filter(upper_body,lower_body,head_and_neck,real_hands,feet,true));
    
    da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->mass.mass);
    da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->mass.mass);

    if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number;
    //else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->id_number;
    //else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->id_number) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->id_number;

   
    
    if(hands){
        Add_Wrist_Joint(true);
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_transform",
            rigid_body_list(1),rigid_body_list(2),1);}
    
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);

    if(parameter_list.Get_Parameter("set_static",false)){
        //da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->is_static=true;
        //da_man->bones(VISIBLE_HUMAN<T>::BONE_R_HUMERUS)->is_static=true;
        //da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ULNA)->is_static=true;
        da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->is_static=true;
    }

    //RIGID_BODY_LIST<T,TV>& rigid_body_list=rigid_body_list;

    
    
    for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->Set_k_p(k_p);
//        if(i==VISIBLE_HUMAN<T>::JOINT_R_STERNOCLAVICULAR || i==VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR || i==VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)
//            joint_function->Set_Target_Angle(joint_function->Angle());
//        else
            joint_function->Set_Target_Angle(QUATERNION<T>());

//        joint_function->muscle_control=(i!=VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL &&
//                                        i!=VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR &&
//                                        i!=VISIBLE_HUMAN<T>::JOINT_R_RADIOULNAR);
//                                        //STERNOCLAVICULAR && 
//                                        //i!=VISIBLE_HUMAN<T>::JOINT_L_STERNOCLAVICULAR && 
//                                        //i!=VISIBLE_HUMAN<T>::JOINT_R_ACROMIOCLAVICULAR && 
//                                        //i!=VISIBLE_HUMAN<T>::JOINT_L_ACROMIOCLAVICULAR);
        joint_function->muscle_control=true;
    }
    
    
    //if(getenv("SET_ZERO")){
    //for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)) da_man->joint(i)->Set_Joint_Frame(FRAME<TV>());
      // arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number);
      //arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->id_number);//}
    }

void Set_Toes_Static()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    for(int i=0;i<arb->rigid_body_list.Number_Of_Elements();i++){
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list.rigid_bodies(i);
        if(rigid_body->name.substr(0,3)=="toe") rigid_body->is_static=true;
    }
}
void Set_Ankle_Joints()
{
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.730149,-(T)0.214719,(T)0.0977022))));
    da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.806459,(T)0.112275,-(T)0.133572))));
    arb->Update_With_Breadth_First_Directed_Graph(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->id_number);
}
//#####################################################################
// Function Add_Wrist_Joint
//#####################################################################
void Add_Wrist_Joint(const bool add_left=false)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    JOINT<TV>* joint=0;
    if(parameter_list.Get_Parameter("bend_joint_for_wrist",false)) joint=new ANGLE_JOINT<TV>();
    else joint=new POINT_JOINT<TV>();
    TV RC_FE_origin_in_anatomical_radius_frame(-(T).01834,-(T).00034,-(T).27811);
    TV RC_FE_axis_in_anatomical_radius_frame((T)0.966275,(T)0.0,-(T)0.257513);
    FRAME<TV> RC_FE_in_anatomical_radius_frame(RC_FE_origin_in_anatomical_radius_frame,da_man->Rotation_From_Axis(RC_FE_axis_in_anatomical_radius_frame,TV(0,1,0)));
    FRAME<TV> parent_frame=da_man->anatomical_frame_in_rigid_body_frame(VISIBLE_HUMAN<T>::BONE_R_RADIUS)*RC_FE_in_anatomical_radius_frame;
    joint->Set_Joint_To_Parent_Frame(parent_frame); //radius
    joint->name="wrist";
    //NOTE: we don't have anatomical frame for hand -- instead we assume the initial configuration of the hand is the reference configuration
    FRAME<TV> child_frame=arb->rigid_body_list.rigid_bodies(1)->Frame().Inverse()*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->Frame()*parent_frame;
    joint->Set_Joint_To_Child_Frame(child_frame); //wrist/hand
    first_wrist_joint_id=arb->joint_mesh.Add_Articulation(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_RADIUS)->particle_index,1,joint);

    if(add_left)  da_man->Add_Additional_Reflected_Left_Joint(arb->joint_mesh.joints.m,da_man->bones(VISIBLE_HUMAN<T>::BONE_L_RADIUS),arb->rigid_body_list(2),1);

    if(parameter_list.Is_Defined("wrist_rotation")){
        TV wrist_rotation=parameter_list.Get_Parameter("wrist_rotation",TV());
        QUATERNION<T> desired_orientation=QUATERNION<T>::From_Rotation_Vector(wrist_rotation);
        QUATERNION<T> desired_orientation_left=QUATERNION<T>::From_Rotation_Vector(parameter_list.Get_Parameter("wrist_rotation_left",wrist_rotation));
        if(add_left){
            arb->joint_mesh.joints(arb->joint_mesh.joints.m-1)->Set_Joint_Frame(FRAME<TV>(desired_orientation));
            arb->joint_mesh.joints(arb->joint_mesh.joints.m)->Set_Joint_Frame(FRAME<TV>(desired_orientation_left));}
        else arb->joint_mesh.joints(arb->joint_mesh.joints.m)->Set_Joint_Frame(FRAME<TV>(desired_orientation));}
}
//#####################################################################
// Function Push_Up
//#####################################################################
void Push_Up()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    T k_p=parameter_list.Get_Parameter("k_p",(T)100);
    add_block=false;

    T skeleton_tilt=parameter_list.Get_Parameter("skeleton_tilt",(T)0);
    TV skeleton_offset=parameter_list.Get_Parameter("skeleton_offset",TV());

    skeleton_frame=FRAME<TV>(skeleton_offset)*FRAME<TV>(QUATERNION<T>(-(T)1.28+skeleton_tilt,TV(1,0,0))*QUATERNION<T>((T)pi,TV(0,1,0))*QUATERNION<T>(-(T)0.5*(T)pi,TV(1,0,0)));

    bool upper_body=parameter_list.Get_Parameter("upper_body",true);
    bool lower_body=parameter_list.Get_Parameter("lower_body",true);
    bool head_and_neck=parameter_list.Get_Parameter("head_and_neck",true);
    bool hands=parameter_list.Get_Parameter("hands",true);
    bool real_hands=parameter_list.Get_Parameter("real_hands",false);
    bool feet=parameter_list.Get_Parameter("feet",true);

    if(hands){
        int id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_right",1,true,false,false);
        RIGID_BODY<TV>* rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_right");
        rigid_body->Frame()=skeleton_frame*rigid_body->Frame();
        //rigid_body->Set_Mass(200);

        id=rigid_body_particles.Add_Rigid_Body(stream_type,data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_left",1,true,false,false);
        rigid_body=arb->rigid_body_list(id);
        rigid_body->Set_Coefficient_Of_Restitution(0);
        rigid_body->Set_Coefficient_Of_Friction(1);
        rigid_body->Set_Name("flat_hand_left");
        rigid_body->Frame()=skeleton_frame*rigid_body->Frame();
        //rigid_body->Set_Mass(200);
    }
    
    bool verbose=parameter_list.Get_Parameter("verbose",false);

    da_man=new VISIBLE_HUMAN<T>(stream_type,arb,solid_body_collection.deformable_object.rigid_body_particles,data_directory,skeleton_frame,verbose);
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Basic_Filter(upper_body,lower_body,head_and_neck,real_hands,feet,true));
   
    da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Mass());
    da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Set_Mass(10*da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Mass());

    if(da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->particle_index;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_THORAX)->particle_index;
    else if(da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index) root=da_man->bones(VISIBLE_HUMAN<T>::BONE_HIP)->particle_index;

//    T elbow_angle=parameter_list.Get_Parameter("elbow_angle",3*(T)pi/4);
//    T shoulder_angle=parameter_list.Get_Parameter("shoulder_angle",(T)pi/4);
    RIGID_BODY_LIST<TV>& rigid_body_list=rigid_body_list;

    if(hands){
        Add_Wrist_Joint(true);
        da_man->Replace_Bones_With_Fused_Bone(data_directory+"/Rigid_Bodies/New_Visible_Human_Bones/flat_hand_transform",rigid_body_list(1),rigid_body_list(2),1);}

    // Set joint frames
    if(upper_body){
        // get these from the tracks
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)->Set_Joint_Frame(shoulder_track->Frame(tracks_start_time));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR)->Set_Joint_Frame(elbow_track->Frame(tracks_start_time));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL)->Set_Joint_Frame(da_man->Rotated_Track(*shoulder_track)->Frame(tracks_start_time));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR)->Set_Joint_Frame(da_man->Rotated_Track(*elbow_track)->Frame(tracks_start_time));

        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_RADIOULNAR)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T)pi+(T).2,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_RADIOULNAR)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T)pi+(T).2,TV(1,0,0))));
    }
    if(feet){
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_TOE_1)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-1,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_TOE_2)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)5.4,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_TOE_3)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)5.4,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_TOE_4)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T).6,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_TOE_5)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T).3,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.6,-(T)0.123,(T)0.0477))));
    //    da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.840725,-(T)0.160726,(T)0.359952))));
        
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_TOE_1)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-1,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_TOE_2)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)5.4,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_TOE_3)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>(-(T)5.4,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_TOE_4)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T).6,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_TOE_5)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>((T).3,TV(1,0,0))));
        da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_ANKLE)->Set_Joint_Frame(FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(TV(-(T)0.689,-(T)0.00996,-(T)0.0822813))));
    }

    arb->Update_With_Breadth_First_Directed_Graph(root);

    // make joint functions to try to keep this pose
    for(int i=0;i<da_man->joint.m;i++) if(da_man->joint(i)){
        da_man->Create_Joint_Function(i);JOINT_FUNCTION<TV>* joint_function=da_man->joint(i)->joint_function;
        joint_function->muscle_control=true;
        joint_function->Set_k_p(k_p);
        joint_function->Set_Target_Angle(joint_function->Angle());
    }

    // add wrist joint functions
    if(hands && parameter_list.Get_Parameter("wrist_joint_function",true)){
        JOINT_FUNCTION<TV>* jfunc=arb->Create_Joint_Function(first_wrist_joint_id);
        jfunc->muscle_control=true;
        jfunc->Set_Target_Angle(jfunc->Angle());jfunc->Set_k_p(k_p);
        jfunc=arb->Create_Joint_Function(first_wrist_joint_id+1);
        jfunc->muscle_control=true;
        jfunc->Set_Target_Angle(jfunc->Angle());jfunc->Set_k_p(k_p);}

    for(int i=0;i<da_man->muscles.m;i++) if(da_man->muscles(i) && !da_man->muscles(i)->peak_force && da_man->muscles(i)->name.find("_left")==std::string::npos){
        T length=da_man->muscles(i)->Total_Length();
        LOG::cout << "cat > " << da_man->muscles(i)->name << ".param << EOF" << std::endl;
        LOG::cout << "// values not available from SIMM model so we made them up by dividing length at rest state into half muscle half tendon" << std::endl;
        LOG::cout << "optimal_fiber_length = " << (T).5*length << std::endl;
        LOG::cout << "tendon_slack_length = " << (T).5*length << std::endl;
        LOG::cout << "peak_force = 1000" << std::endl;
        LOG::cout << "pennation_angle = " << 0 << std::endl;
        LOG::cout << "EOF" << std::endl;
    }

    for(int i=0;i<da_man->muscles.m;i++) if(da_man->muscles(i)){
        T length=da_man->muscles(i)->Total_Length();
        T total_rest_length=da_man->muscles(i)->optimal_length + da_man->muscles(i)->tendon_slack_length;
        if(total_rest_length>(T)1.5*length){
            std::cout << "Muscle " << da_man->muscles(i)->name << ": (" << da_man->muscles(i)->optimal_length << "+" << da_man->muscles(i)->tendon_slack_length << ") " << total_rest_length << " vs " << length << std::endl;}}

//    T peak_force=parameter_list.Get_Parameter("peak_force",(T)1);
//    for(int i=0;i<da_man->muscles.m;i++) if(da_man->muscles(i)) da_man->muscles(i)->Set_Peak_Force(peak_force);

    // Determine which joints will get PD as opposed to muscle actuation
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

    if(arb->use_muscle_actuators){
        LOG::cout << "Joints with only muscle control"<<std::endl;
        for(int i=0;i<arb->joint_mesh.joints.m;i++) if(arb->joint_mesh.joints(i)->joint_function && arb->joint_mesh.joints(i)->joint_function->muscle_control)
            LOG::cout << "\t" << arb->joint_mesh.joints(i)->name << std::endl;}

    for(int i=0;i<rigid_body_list.rigid_bodies.m;i++) rigid_body_list.rigid_bodies(i)->Set_Coefficient_Of_Friction(1);
}
//#####################################################################
// Function Read_Mocap_Data
//#####################################################################
void Read_Mocap_Data()
{
    std::istream* input=FILE_UTILITIES::Safe_Open_Input("../../Tools/c3d2motion/motion",false);
    MOTION_SEQUENCE<float,QUATERNION<float> > my_motion;
    Read_Binary<float>(*input,my_motion);
    skeleton_frame=FRAME<TV>(TV(0,(T).05,0),QUATERNION<T>((T)pi/2,TV(1,0,0))*QUATERNION<T>((T)pi,TV(0,1,0)));
    da_man=new VISIBLE_HUMAN<T>(arb,data_directory,skeleton_frame);
    root=VISIBLE_HUMAN<T>::BONE_CRANIUM;
    da_man->Initialize_Bodies(VISIBLE_HUMAN<T>::Basic_Filter(true,true,true,true,true,true));
    da_man->bones(VISIBLE_HUMAN<T>::BONE_CRANIUM)->is_static=false;
    /*da_man->bones(VISIBLE_HUMAN<T>::BONE_L_ANKLE)->is_static=true;
      da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ANKLE)->is_static=true;*/
    T k_p=parameter_list.Get_Parameter("k_p",(T)50);
    for(int i=0;i<da_man->joint.m;i++){
        std::cout<<"joint name is "<<da_man->joint(i)->name<<"\n";
        int track = my_motion.Track_Index(da_man->joint(i)->name);
        if (track != -1){
            da_man->Create_Joint_Function(i);
            da_man->joint(i)->joint_function->track=new FRAME_TRACK_3D<T>(my_motion.trajectories(1).m,my_motion.time_grid.xmin, my_motion.time_grid.xmax);
            da_man->joint(i)->joint_function->Set_k_p(k_p);
            for (int j = 1; j<= my_motion.trajectories(1).m; j++){
                if (da_man->joint(i)->name == "joint_r_knee" || da_man->joint(i)->name == "joint_r_knee_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)pi-(T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_r_humeroulnar" || da_man->joint(i)->name == "joint_r_humeroulnar_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)pi+(T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_r_radiocarpal" || da_man->joint(i)->name == "joint_r_radiocarpal_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)pi-(T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_r_hip" || da_man->joint(i)->name == "joint_r_hip_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)my_motion.trajectories(track)(j).Angle()-4*(T)pi/5,TV(my_motion.trajectories(track)(1).v)));
                if (da_man->joint(i)->name == "joint_r_ankle" || da_man->joint(i)->name == "joint_r_ankle_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>(5*(T)pi/8-(T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_r_glenohumeral" || da_man->joint(i)->name == "joint_r_glenohumeral_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_atlas")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>(3*(T)pi/4-(T)my_motion.trajectories(track)(j).Angle(),TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_l1")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)pi-(T)my_motion.trajectories(track)(j).Angle(), TV((T)my_motion.trajectories(track)(j).v.x, 0, 0)));
                if (da_man->joint(i)->name == "joint_l5_coccyx")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)my_motion.trajectories(track)(j).Angle()-(T)pi,TV((T)0.78,-(T)0.22,0)));
                if (da_man->joint(i)->name == "joint_r_sternoclavicular" || da_man->joint(i)->name == "joint_r_sternoclavicular_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)my_motion.trajectories(track)(j).Angle()-(T)pi/3,TV(my_motion.trajectories(track)(j).v)));
                if (da_man->joint(i)->name == "joint_r_acromioclavicular" || da_man->joint(i)->name == "joint_r_acromioclavicular_left")
                    da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>((T)my_motion.trajectories(track)(j).Angle()+3*(T)pi/2,TV(my_motion.trajectories(track)(j).v)));}}
        else{
            da_man->Create_Joint_Function(i);
            da_man->joint(i)->joint_function->track=new FRAME_TRACK_3D<T>(my_motion.trajectories(1).m,my_motion.time_grid.xmin, my_motion.time_grid.xmax);
            da_man->joint(i)->joint_function->Set_k_p(k_p);
            for (int j = 1; j<= my_motion.trajectories(1).m; j++) da_man->joint(i)->joint_function->track->trajectory(j)=FRAME<TV>(QUATERNION<T>(0,TV()));}}
}
//#####################################################################
// Function Create_Push_Up_Tracks
//#####################################################################
void Create_Push_Up_Tracks()
{
    int samples=1000;T period=5;
    shoulder_track=new FRAME_TRACK_3D<T>(samples,0,period);shoulder_track->periodic=true;
    elbow_track=new FRAME_TRACK_3D<T>(samples,0,period);elbow_track->periodic=true;
    T shoulder_start=parameter_list.Get_Parameter("shoulder_start",(T)0),shoulder_end=parameter_list.Get_Parameter("shoulder_end",(T)pi/4);
    T elbow_start=parameter_list.Get_Parameter("elbow_start",3*(T)pi/4),elbow_end=parameter_list.Get_Parameter("elbow_end",(T).2);
    TV shoulder_prerot=parameter_list.Get_Parameter("shoulder_prerotation",TV());
    TV shoulder_postrot=parameter_list.Get_Parameter("shoulder_postrotation",TV());
    for(int i=0;i<samples;i++){
        T t=(T)(i-1)/(samples-1) - tracks_start_time/period; // normally starts at zero but offset by time at which motion is starting
        shoulder_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(shoulder_prerot)*
                                                  QUATERNION<T>(shoulder_start+(shoulder_end-shoulder_start)*(T)0.5*(1-cos(2*(T)pi*t)),TV(1,0,0))*
                                                  QUATERNION<T>::From_Rotation_Vector(shoulder_postrot));
        elbow_track->trajectory(i)=FRAME<TV>(QUATERNION<T>(elbow_start+(elbow_end-elbow_start)*(T)0.5*(1-cos(2*(T)pi*t)),TV(1,0,0)));
    }
}
//#####################################################################
// Function Create_Skeleton_In_Flesh_Tracks
//#####################################################################
void Create_Skeleton_In_Flesh_Tracks()
{
    int samples=1000;T period=2;
    shoulder_track=new FRAME_TRACK_3D<T>(samples,0,period);shoulder_track->periodic=true;
    elbow_track=new FRAME_TRACK_3D<T>(samples,0,period);elbow_track->periodic=true;
    T shoulder_start=parameter_list.Get_Parameter("shoulder_start",(T)0),shoulder_end=parameter_list.Get_Parameter("shoulder_end",(T)pi/4);
    T elbow_start=parameter_list.Get_Parameter("elbow_start",(T)0),elbow_end=parameter_list.Get_Parameter("elbow_end",(T)pi/4);
    TV shoulder_prerot=parameter_list.Get_Parameter("shoulder_prerotation",TV());
    TV shoulder_postrot=parameter_list.Get_Parameter("shoulder_postrotation",TV());
    for(int i=0;i<samples;i++){
        T t=(T)(i-1)/(samples-1) - tracks_start_time/period; // normally starts at zero but offset by time at which motion is starting
        shoulder_track->trajectory(i)=FRAME<TV>(QUATERNION<T>::From_Rotation_Vector(shoulder_prerot)*
                                                  QUATERNION<T>(shoulder_start+(shoulder_end-shoulder_start)*(T)0.5*(1-cos(2*(T)pi*t)),TV(1,0,0))*
                                                  QUATERNION<T>::From_Rotation_Vector(shoulder_postrot));
        elbow_track->trajectory(i)=FRAME<TV>(QUATERNION<T>(elbow_start+(elbow_end-elbow_start)*(T)0.5*(1-cos(2*(T)pi*t)),TV(1,0,0)));
    }
}
//#####################################################################
// Function Add_Tracks
//#####################################################################
void Add_Tracks(const T time)
{
    LOG::cout << "Adding tracks (time " << time << ")" << std::endl;
    if(use_tracks){
        if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL))
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_GLENOHUMERAL)->joint_function->track=shoulder_track;
        if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR))
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_R_HUMEROULNAR)->joint_function->track=elbow_track;
        if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL))
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_GLENOHUMERAL)->joint_function->track=da_man->Rotated_Track(*shoulder_track);
        if(da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR))
            da_man->joint(VISIBLE_HUMAN<T>::JOINT_L_HUMEROULNAR)->joint_function->track=da_man->Rotated_Track(*elbow_track);}
    tracks_added=true;
}
//#####################################################################
// Function Add_Hand_Ground_Joints
//#####################################################################
void Add_Hand_Ground_Joints(const int right_hand_index,const int left_hand_index,const T offset=0)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    //T k_p=parameter_list.Get_Parameter("k_p",(T)100);
    FRAME<TV> vertical_offset(TV());
    RIGID_BODY<TV>* ground=arb->rigid_body_list.rigid_bodies(arb->rigid_body_list.rigid_bodies.m);assert(ground->name=="ground");
    
    JOINT<TV>* joint=new ANGLE_JOINT<TV>();
    //joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*(arb->rigid_body_list.rigid_bodies(right_hand_index)->Frame()+vertical_offset));
    joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*(FRAME<TV>((arb->rigid_body_list.rigid_bodies(right_hand_index)->Frame()).t+vertical_offset.t,(arb->rigid_body_list.rigid_bodies(right_hand_index)->Frame()).r)));
    arb->joint_mesh.Add_Articulation(ground->particle_index,arb->rigid_body_list.rigid_bodies(right_hand_index)->particle_index,joint);
    //JOINT_FUNCTION<TV>* joint_function=arb->Create_Joint_Function(arb->joint_mesh.joints.m);
    //joint_function->Set_k_p(k_p);
    //joint_function->Set_Target_Angle(joint_function->Angle());

    joint=new ANGLE_JOINT<TV>();
    //joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*(arb->rigid_body_list.rigid_bodies(left_hand_index)->Frame()+vertical_offset));
    joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*(FRAME<TV>((arb->rigid_body_list.rigid_bodies(left_hand_index)->Frame()).t+vertical_offset.t,(arb->rigid_body_list.rigid_bodies(left_hand_index)->Frame()).r)));
    arb->joint_mesh.Add_Articulation(ground->particle_index,arb->rigid_body_list.rigid_bodies(left_hand_index)->particle_index,joint);
    //joint_function=arb->Create_Joint_Function(arb->joint_mesh.joints.m);
    //joint_function->Set_k_p(k_p);
    //joint_function->Set_Target_Angle(joint_function->Angle());
}
//#####################################################################
// Function Add_Foot_Ground_Joints
//#####################################################################
void Add_Foot_Ground_Joints()
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    RIGID_BODY<TV>* ground=arb->rigid_body_list.rigid_bodies(arb->rigid_body_list.rigid_bodies.m);
    JOINT<TV>* joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*FRAME<TV>(da_man->bones(VISIBLE_HUMAN<T>::BONE_R_TOE_1)->Frame().t));
    arb->joint_mesh.Add_Articulation(ground->particle_index,da_man->bones(VISIBLE_HUMAN<T>::BONE_R_ANKLE)->particle_index,joint);
    
    joint=new POINT_JOINT<TV>();
    joint->Set_Joint_To_Parent_Frame(ground->Frame().Inverse()*FRAME<TV>(da_man->bones(VISIBLE_HUMAN<T>::BONE_L_TOE_1)->Frame().t));
    arb->joint_mesh.Add_Articulation(ground->particle_index,da_man->bones(VISIBLE_HUMAN<T>::BONE_L_ANKLE)->particle_index,joint);
}
//#####################################################################
// Function Update_Joints
//#####################################################################
void Update_Joints(const /*int frame*/T time)
{
    for(int i=0;i<arb->joint_mesh.joints.m;i++)
        if(arb->joint_mesh.joints(i)->joint_function && arb->joint_mesh.joints(i)->joint_function->track)
            arb->joint_mesh.joints(i)->Set_Joint_Frame(arb->joint_mesh.joints(i)->joint_function->track->Frame(time));
    arb->Update_With_Breadth_First_Directed_Graph(root);
}

//#####################################################################
// Function Preprocess_frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
   if (frame > 0 && !simulate)  Update_Joints((T)frame);
}

//#####################################################################
// Function Apply_Constraints
//#####################################################################
// used in test_frame_tracks case to dump the animated body without simulating
void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE
{
    BASE::Apply_Constraints(dt,time);
    if(test_frame_tracks) Update_Joints(time);
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
// for restart!!
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    LOG::cout << "*** IN Read_Output_Files_Solids" << std::endl;
    Save_Masses();
    BASE::Read_Output_Files_Solids(frame);
    Restore_Masses();
}
void Save_Masses()
{
    if(!da_man) return;
    LOG::cout << "--- saving masses" << std::endl;
    saved_masses.Resize(da_man->bones.m);
    // TODO: confirm this mass is world-space
    for(int i=0;i<da_man->bones.m;i++) if(da_man->bones(i))
        saved_masses(i)=da_man->bones(i)->Mass();
}
void Restore_Masses()
{
    if(!da_man) return;
    LOG::cout << "--- restoring masses" << std::endl;
    for(int i=0;i<da_man->bones.m;i++) if(da_man->bones(i)){
        da_man->bones(i)->Mass()=saved_masses(i);}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    if(use_tracks && !tracks_added && time>=tracks_start_time) Add_Tracks(time);
    if(da_man&&arb) da_man->Turn_Off_Collisions(arb,solids_evolution->rigid_body_collisions->collision_manager);
}
//#####################################################################
// Function Add_Deformable_Body
//#####################################################################
void Add_Deformable_Body(const std::string& filename,const GEOMETRY_TYPE type,RIGID_BODY_STATE<TV>* initial_state)
{
    number_of_deformable_bodies++;
    deformable_body_geometry_filenames.Append(filename);
    deformable_body_geometry_types.Append(type);
    deformable_body_initial_states.Append(initial_state);
}
//#####################################################################
// Function Initialize_Tetrahedron_Collisions
//#####################################################################
void Initialize_Tetrahedron_Collisions(const int id_number,TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,TRIANGULATED_SURFACE<T>* triangulated_surface=0)
{
    TETRAHEDRON_COLLISION_BODY<T>& tetrahedron_collision_body=*(new TETRAHEDRON_COLLISION_BODY<T>(tetrahedralized_volume,triangulated_surface));
    DEFORMABLE_PARTICLES<TV>& undeformed_particles=*(new DEFORMABLE_PARTICLES<TV>(tetrahedralized_volume.particles));
    TRIANGULATED_SURFACE<T>& undeformed_triangulated_surface=*(new TRIANGULATED_SURFACE<T>(tetrahedron_collision_body.triangulated_surface->mesh,undeformed_particles));
    undeformed_triangulated_surface.Update_Triangle_List();undeformed_triangulated_surface.Initialize_Triangle_Hierarchy();

    // undeformed levelset
    LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*LEVELSET_IMPLICIT_OBJECT<TV>::Create();
    std::string levelset_filename=output_directory+STRING_UTILITIES::string_sprintf("/test_%d_deformable_body_undeformed_levelset_%d.phi",6,id_number);
    if(FILE_UTILITIES::File_Exists(levelset_filename)) FILE_UTILITIES::Read_From_File(stream_type,levelset_filename,undeformed_levelset);
    else{
        undeformed_triangulated_surface.Update_Bounding_Box();BOX_3D<T>& box=*undeformed_triangulated_surface.bounding_box;
        GRID<TV>& grid=undeformed_levelset.levelset.grid;ARRAY<T,VECTOR<int,3> >& phi=undeformed_levelset.levelset.phi;
        grid=GRID_3D<T>::Create_Grid_Given_Cell_Size(box,(T)1e-2*box.Edge_Lengths().Max(),false,5);
        phi.Resize(grid);
        LEVELSET_MAKER_UNIFORM<T> levelset_maker;
        levelset_maker.Verbose_Mode();
        levelset_maker.Set_Surface_Padding_For_Flood_Fill((T)1e-3);
        levelset_maker.Use_Fast_Marching_Method(true,0);
        levelset_maker.Compute_Level_Set(undeformed_triangulated_surface,grid,phi);
        FILE_UTILITIES::Create_Directory(output_directory);FILE_UTILITIES::Write_To_File(stream_type,levelset_filename,undeformed_levelset);} 
    undeformed_levelset.Update_Box();
    tetrahedron_collision_body.Set_Implicit_Geometry(&undeformed_levelset);
    tetrahedron_collision_body.Set_Undeformed_Triangulated_Surface(&undeformed_triangulated_surface);
    solids_parameters.collision_body_list.Add_Body(&tetrahedron_collision_body);
}
//#####################################################################
// Function Create_Tetrahedralized_Volume
//#####################################################################
STRUCTURE<TV>* Create_Tetrahedralized_Volume(int index)
{
    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=*TETRAHEDRALIZED_VOLUME<T>::Create();
    DEFORMABLE_PARTICLES<TV>& particles=tetrahedralized_volume.particles;
    assert(deformable_body_geometry_filenames(index)!="");
    FILE_UTILITIES::Read_From_File(stream_type,deformable_body_geometry_filenames(index),tetrahedralized_volume);
    LOG::cout<<"Deformable body "<<index<<" - Total Tetrahedra : "<<tetrahedralized_volume.mesh.elements.m<<std::endl;
    particles.Store_Velocity();particles.Store_Mass();
    tetrahedralized_volume.Set_Density(1000);tetrahedralized_volume.Set_Mass_Of_Particles(true);
    Set_Initial_Particle_Configuration(particles,index);
    return &tetrahedralized_volume;
}
//#####################################################################
// Function Set_Initial_Particle_Configuration
//#####################################################################
void Set_Initial_Particle_Configuration(DEFORMABLE_PARTICLES<TV>& particles,const int index)
{
    if(deformable_body_initial_states(index)){
        LOG::cout<<"Deformable body "<<index<<" - Total Particles : "<<particles.Size()<<std::endl;
        BOX_3D<T> bounding_box(particles.X(1));for(int i=2;i<=particles.Size();i++) bounding_box.Enlarge_To_Include_Point(particles.X(i));TV center=bounding_box.Center();
        RIGID_BODY_STATE<TV>& state=*deformable_body_initial_states(index);
        for(int p=0;p<particles.Size();p++){
            particles.X(p)=state.frame*(particles.X(p)-center);
            particles.V(p)=state.twist.linear+TV::Cross_Product(state.twist.angular,particles.X(p)-state.frame.t);}}
}
//#####################################################################
// Function Get_Constrained_Particle_Data
//#####################################################################
void Get_Constrained_Particle_Data()
{
    assert(draw_flesh);

    bone_ids.Resize(num_bones_present=0);
    for(int i=0;i<VISIBLE_HUMAN<T>::num_bones;i++)
        if(da_man->bones(i)->particle_index) {bone_ids.Append(da_man->bones(i)->particle_index); num_bones_present++;}

    DEFORMABLE_PARTICLES<TV>& particles=(solid_body_collection.deformable_object.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>()).particles;
 
    enslaved_nodes.Resize(num_bones_present);
    positions_relative_to_bone_frames.Resize(num_bones_present);

    for(int i=0;i<particles.Size();i++)
        for(int p=0;p<num_bones_present;p++){
            RIGID_BODY<TV>& bone_rigid_body=*arb->rigid_body_list.rigid_bodies(bone_ids(p));
            if (!bone_rigid_body.Implicit_Geometry_Lazy_Outside(particles.X(i))){
                enslaved_nodes(p).Append(i);
                positions_relative_to_bone_frames(p).Append(bone_rigid_body.Frame().Inverse_Times(particles.X(i)));}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    assert(fragment_id==FRAGMENT_ID(1));
    for(int p=0;p<num_bones_present;p++)for(int i=0;i<enslaved_nodes(p).m;i++){
        X(enslaved_nodes(p)(i))=TV();} 
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.deformable_object.rigid_body_particles;
    assert(fragment_id==FRAGMENT_ID(1));
    for(int p=0;p<num_bones_present;p++){
        RIGID_BODY<TV>& bone_rigid_body=*arb->rigid_body_list.rigid_bodies(bone_ids(p));
        FRAME<TV> inverted_frame=bone_rigid_body.Frame().Inverse();
        for(int i=0;i<enslaved_nodes(p).m;i++){
            X(enslaved_nodes(p)(i))=inverted_frame.Inverse_Times(positions_relative_to_bone_frames(p)(i));}}
}
//#####################################################################
// Function Untweak_Joint
//#####################################################################
void Untweak_Joint(RIGID_BODY<TV>* parent,RIGID_BODY<TV>* child,JOINT<TV>* joint) PHYSBAM_OVERRIDE
{
    joint->Set_Joint_To_Child_Frame(child->Frame().Inverse()*parent->Frame()*joint->F_pj());
}
//#####################################################################
// Function Augment_Angle - NOTE: NOT COMUTATIVE
//#####################################################################
void Augment_Child_Angle_X(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_twist=QUATERNION<T>(angle,TV(1,0,0));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_twist));
}
void Augment_Child_Angle_Y(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_rot=QUATERNION<T>(angle,TV(0,1,0));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_rot));
}
void Augment_Child_Angle_Z(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_bend=QUATERNION<T>(angle,TV(0,0,1));
    FRAME<TV> child_frame=joint->F_cj(); 
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(child_frame.t,child_frame.r*pre_bend));
}
void Augment_Parent_Angle_X(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_twist=QUATERNION<T>(angle,TV(1,0,0));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
void Augment_Parent_Angle_Y(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_twist=QUATERNION<T>(angle,TV(0,1,0));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
void Augment_Parent_Angle_Z(JOINT<TV>* joint,T angle)
{
    QUATERNION<T> pre_twist=QUATERNION<T>(angle,TV(0,0,1));
    FRAME<TV> parent_frame=joint->F_pj(); 
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(parent_frame.t,parent_frame.r*pre_twist));
}
//#####################################################################
};
}
#endif
