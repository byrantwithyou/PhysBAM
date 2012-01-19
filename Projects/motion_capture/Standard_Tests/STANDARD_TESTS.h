//#####################################################################
// Copyright 2008, Michael Lentine
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Two bodies with a bend joint and gravity
//   2. Two bodies with a bend joint and constant wind
//   3. Two bodies with a bend joint and varying wind
//   4. Two bodies with a point joint and gravity
//   5. Two bodies with a point joint and constant wind
//   6. Two bodies with a point joint and varying wind
//   7. Human Legs and a spine (1 contolled dof)
//   8. Human Legs and a spine (multiple controlled dof)
//   9. Human with no spine control
//   10. Human
//   11. Cluster Human Top
//   11. Cluster Human Spine
//   14. Octosquid
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Data_Structures/READ_WRITE_TRIPLE.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Geometry/Grids_Uniform_Computations/DUALCONTOUR_3D.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Read_Write/Geometry/READ_WRITE_TRIANGULATED_SURFACE.h>
#include <PhysBAM_Geometry/Read_Write/Implicit_Objects_Uniform/READ_WRITE_LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Motion/BODY_MOTION_SEQUENCE.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_Evolution/SEARCH_CONTROLLER.h>
#include <string>

#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>

#define ANGLE_JOINT_TYPE 1
#define POINT_JOINT_TYPE 2

namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef GRID<TV> T_GRID;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_POLICY<GRID<TV> >::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename GRID<TV>::CELL_ITERATOR CELL_ITERATOR;
    typedef typename RIGID_BODY_POLICY<TV>::INERTIA_TENSOR T_INERTIA_TENSOR;
    typedef VECTOR<int,3> TV_INT;
    typedef typename TV::SPIN T_SPIN;

    typedef typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER T_CLUSTER;
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
public:
    using BASE::Add_Volumetric_Body_To_Fluid_Simulation;
    SOLIDS_FLUIDS_DRIVER<TV>* driver;
    SOLIDS_STANDARD_TESTS<TV> tests;
    ARRAY<int>* referenced_rigid_particles;
    RIGID_BODY_COLLISION_MANAGER_HASH* collision_manager;
    ROTATION<TV> rotation;
    SEARCH_CONTROLLER<T_GRID>* controller;
    BODY_MOTION_SEQUENCE<T> body_motion;
    ARRAY<ARRAY<int> > bone_hierarchy;
    ARRAY<std::string> controlled_bones;
    ARRAY<std::string> uncontrolled_bones;
    ARRAY<ARRAY<std::string> > cluster_bones;
    FRAME<TV> rigid_base_transform;
    ARRAY<int,int> id_to_index;
    ARRAY<int> rigid_body_ids;
    int motion_frame_rate;
    bool kinematic_motion,use_clustering,use_limits,use_motion_capture,use_deformable,use_embedding;
    RIGID_BODY<TV>* octosquid_body;
    GRAVITY<TV> *solids_source;
    //spatially varying wind
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<TV,TV_INT> v_array;
    ARRAY<T,TV_INT> p_array,d_array;
    GRID<TV> grid;
    std::string input_directory;

    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::frame_rate;using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::parse_args;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),driver(0),tests(*this,solid_body_collection),referenced_rigid_particles(0),collision_manager(0),kinematic_motion(true),
        use_clustering(false),use_limits(false),use_motion_capture(true),use_deformable(false),use_embedding(false),octosquid_body(0),solids_source(0)
    {
    }

    ~STANDARD_TESTS()
    {}
    
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {
        if(controller && controller->hypothetical_step) dt=controller->dt_hyp;
        else if(test_number==14) dt=1/frame_rate+(T)1e-4;}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_String_Argument("-mtn","Specify the physbam motion file to use");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.use_rigid_deformable_contact=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
    solids_parameters.cfl=1;
    solids_parameters.implicit_solve_parameters.cg_iterations=1000;

    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    input_directory="/n/field/disk2/mlentine/PhysBAM/Projects/fluids_3d/Standard_Tests_Smoke/Test_1__Resolution_50_50_75/";
    frame_rate=30;motion_frame_rate=60;

    if(test_number<6||test_number>13){use_motion_capture=false;use_deformable=true;use_embedding=true;}
    if(test_number==13){use_embedding=false;use_deformable=true;}

    if(parse_args->Is_Value_Set("-mtn")) FILE_UTILITIES::Read_From_File(stream_type,parse_args->Get_String_Value("-mtn"),body_motion);
    else if((test_number>6&&test_number<14)||use_motion_capture) PHYSBAM_FATAL_ERROR("You must supply a mtn file to simulate a human.");
        
    switch(test_number){
        case 1: last_frame=360;break;
        case 2: last_frame=360;break;
        case 3: last_frame=360;break;
        case 4: last_frame=360;break;
        case 5: last_frame=360;break;
        case 6: last_frame=360;break;
        case 7: last_frame=620;break;
        case 8: last_frame=620;break;
        case 9: last_frame=620;break;
        case 10: last_frame=620;break;
        case 11: last_frame=620;break;
        case 12: last_frame=620;break;
        case 13: last_frame=620;break;
        case 14: last_frame=620;break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
void Set_Driver(SOLIDS_FLUIDS_DRIVER<TV>* driver_input)
{
    driver=driver_input;
}
void Setup_Spatially_Varying_Wind(const int frame)
{
    std::string f=STRING_UTILITIES::string_sprintf("%d",13);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/mac_velocities."+f,face_velocities);
    v_array.Resize(face_velocities.Component(1).Domain_Indices(),false,false);
    p_array.Resize(face_velocities.Component(1).Domain_Indices(),false,false);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/"+f+"/pressure",p_array);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/"+f+"/density.",d_array);
    if(frame==0){
        FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/common/grid",grid);
        grid.Initialize(face_velocities.Component(1).counts.x,face_velocities.Component(1).counts.y,face_velocities.Component(1).counts.z,-5.0,5.0,0.0,10.0,-5.0,5.0);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/common/grid",grid);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        v_array(cell_index)=TV((face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x+1,cell_index.y,cell_index.z)))/(T)2,
            (face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x,cell_index.y+1,cell_index.z)))/(T)2,
            (face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x,cell_index.y,cell_index.z+1)))/(T)2);
    }
    for(int i=1;WIND_DRAG_3D<T>* drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>*>(i);i++){
        drag->Use_Spatially_Varying_Wind((T)0,grid,v_array);drag->Set_Wind_Pressure(p_array);}
}
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{  
    if(test_number==14){
        T cycle_time=time*4; //.5 second cycle
        T force_magnitude=3;
        T alpha=(cycle_time-7*(((int)cycle_time)/7))/4;
        if(((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==2 || ((int)cycle_time)%7==3) force_magnitude*=-1*sin(alpha*(T)pi);
        else force_magnitude=0;
        solids_source->Set_Gravity(TV(0,force_magnitude,0));}
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    if(test_number==14 && controller && !controller->hypothetical_step){
        T cycle_time=4*time;
        controller->drag_direction=octosquid_body->Rotation().Rotate(TV(0,-1,0));
        bool old_min=controller->minimize;
        controller->minimize=((((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==6)?false:true);
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(((int)cycle_time)%7==0 || ((int)cycle_time)%7==1 || ((int)cycle_time)%7==6) joint.joint_function->Set_k_p(750);
            else joint.joint_function->Set_k_p(1500);}
        if(old_min ^ controller->minimize){
            LOG::cout<<"Switching, so filling dF_array_multipliers with PAIR(0,1)"<<std::endl;
            controller->dF_array_multipliers.Fill(PAIR<T,T>((T)0,(T)1));}}
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    BASE::Write_Output_Files(frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    LOG::cout<<"Writing "<<output_directory+"/mac_velocities."+f<<std::endl;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/density.",d_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/mac_velocities",face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/"+f+"/pressure",p_array);
    controller->Write(stream_type,output_directory,frame);
}
//#####################################################################
// Function Set_PD_Targets
//#####################################################################
void Set_PD_Targets(const T dt,const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    //RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    int frame=(int)(time*motion_frame_rate)+1;
    T alpha=time*motion_frame_rate-frame+1;
    if(use_motion_capture && !kinematic_motion){
       //motion capture input
       for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            bool controlled=false;for(int j=0;j<T_SPIN::dimension;j++) if(joint.control_dof(j)) controlled=true;
            if(!joint.joint_function || !controlled) continue;
            RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
            int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
            FRAME<TV> rigid_base_transform_parent1=rigid_base_transform,rigid_base_transform_child1=rigid_base_transform;
            FRAME<TV> rigid_base_transform_parent2=rigid_base_transform,rigid_base_transform_child2=rigid_base_transform;
            rigid_base_transform_parent1.t*=body_motion.trajectories(parent_index)(frame).length;
            rigid_base_transform_child1.t*=body_motion.trajectories(child_index)(frame).length;
            rigid_base_transform_parent2.t*=body_motion.trajectories(parent_index)(frame+1).length;
            rigid_base_transform_child2.t*=body_motion.trajectories(child_index)(frame+1).length;
            FRAME<TV> parent_frame=FRAME<TV>::Interpolation(body_motion.trajectories(parent_index)(frame).targeted_transform*rigid_base_transform_parent1,
                body_motion.trajectories(parent_index)(frame+1).targeted_transform*rigid_base_transform_parent2,alpha);
            FRAME<TV> child_frame=FRAME<TV>::Interpolation(body_motion.trajectories(child_index)(frame).targeted_transform*rigid_base_transform_child1,
                body_motion.trajectories(child_index)(frame+1).targeted_transform*rigid_base_transform_child2,alpha);
            joint.joint_function->Set_Target_Angle((joint.frame_jp*parent_frame.Inverse()*child_frame*joint.frame_cj).r);}}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(test_number!=3 && test_number!=6) return;
    Setup_Spatially_Varying_Wind(frame);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    //Write_Output_Files(frame);
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(controller && !controller->hypothetical_step) controller->Update_Position_Based_State(fluid_collection.incompressible_fluid_collection.face_velocities,dt,time);
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    //DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    //deformable_body_collection.particles.V*=0;solid_body_collection.rigid_body_collection.rigid_body_particle.twist*=0;solid_body_collection.rigid_body_collection.rigid_body_particle.angular_momentum*=0;
}
//#####################################################################
// Function Post_Initialization
//#####################################################################
void Post_Initialization() PHYSBAM_OVERRIDE
{
    solids_evolution->rigid_body_collisions->collision_manager=collision_manager;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    if(use_motion_capture && kinematic_motion && id_to_index(id)<=body_motion.trajectories.m){
        int frame_num=(int)(time*motion_frame_rate)+1;
        T alpha=time*motion_frame_rate-frame_num+1;
        RIGID_BODY<TV>* parent=&arb.rigid_body_collection.Rigid_Body(id);
        int parent_index=id_to_index(parent->particle_index);
        FRAME<TV> rigid_base_transform_parent1=rigid_base_transform,rigid_base_transform_child1=rigid_base_transform;
        FRAME<TV> rigid_base_transform_parent2=rigid_base_transform,rigid_base_transform_child2=rigid_base_transform;
        rigid_base_transform_parent1.t*=body_motion.trajectories(parent_index)(frame_num).length;
        rigid_base_transform_parent2.t*=body_motion.trajectories(parent_index)(frame_num+1).length;
        frame=FRAME<TV>::Interpolation(body_motion.trajectories(parent_index)(frame_num).targeted_transform*rigid_base_transform_parent1,
            body_motion.trajectories(parent_index)(frame_num+1).targeted_transform*rigid_base_transform_parent2,alpha);}
    else{
        RIGID_BODY<TV>* kinematic_body=&solid_body_collection.rigid_body_collection.Rigid_Body(id);
        frame=kinematic_body->Frame();}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    controller=new SEARCH_CONTROLLER<T_GRID>(solid_body_collection,driver);
    controller->solve_minimization=false;
    controller->max_iterations=1;
    controller->minimize=true;
    controller->use_projection=false;
    controller->drag_direction=TV(0,0,-1);
    controller->dt_per_search_step=(T).1;
    controller->dt_hyp=(T).05;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(10);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(10);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=10;
    arb.constrain_pd_directions=true;
    
    T height=0;
    
    //set up rigid_bodies and external forces
    switch(test_number){
        case 1: 
        case 2: 
        case 3: 
        case 4: 
        case 5: 
        case 6: Cylinder_And_Block();break;
        case 7:
        case 8:
        case 9:
        case 10: Human(1); break;
        case 11: 
        case 12: 
        case 13: 
            if(test_number==11){
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(1).Append("Spine");cluster_bones(1).Append("Neck");cluster_bones(1).Append("Head");
                cluster_bones(1).Append("Clavicle");cluster_bones(1).Append("Arm");cluster_bones(1).Append("Hand");
                controlled_bones.Append("Spine");}
            else if(test_number==12){
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(1).Append("LUpper_Arm");cluster_bones(1).Append("LLower_Arm");cluster_bones(1).Append("LHand");
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(2).Append("RUpper_Arm");cluster_bones(2).Append("RLower_Arm");cluster_bones(2).Append("RHand");
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(3).Append("Spine");cluster_bones(3).Append("Neck");cluster_bones(3).Append("Head");cluster_bones(3).Append("Clavicle");
                controlled_bones.Append("Arm");controlled_bones.Append("Spine");/*controlled_bones.Append("Head");controlled_bones.Append("Neck");controlled_bones.Append("Clavicle");*/}
            else if(test_number==13){
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(1).Append("LUpper_Arm");cluster_bones(1).Append("LLower_Arm");cluster_bones(1).Append("LHand");
                cluster_bones.Append(ARRAY<std::string>());
                cluster_bones(2).Append("RUpper_Arm");cluster_bones(2).Append("RLower_Arm");cluster_bones(2).Append("RHand");
                controlled_bones.Append("Arm");/*controlled_bones.Append("Head");controlled_bones.Append("Neck");controlled_bones.Append("Clavicle");*/}
            Human_Cluster(1); 
            break;
        case 14:
            controller->drag_direction=TV(0,(T)-1,0);
            controller->real_dx=(T)20*(T)pi/180;
            controller->dx=(T)10*(T)pi/180;
            use_deformable=true;use_embedding=false;
            controller->dt_per_search_step=(T).1;
            controller->dt_hyp=(T).0333333;
            controller->use_projection=false;
            Octosquid();
            break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(use_deformable){
        bool build_tri=false,build_tet=false;
        TETRAHEDRALIZED_VOLUME<T> *volume_original=TETRAHEDRALIZED_VOLUME<T>::Create();
        TRIANGULATED_SURFACE<T> *surface_original=0,*surface=0;
        if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number)) && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number))){
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);
            if(use_embedding){
                surface_original=TRIANGULATED_SURFACE<T>::Create();
                FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}}
        else if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number))){
            if(use_embedding) build_tri=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);}
        else if(use_embedding && FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number))){
            surface_original=TRIANGULATED_SURFACE<T>::Create();
            build_tet=true;FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}
        else{if(use_embedding) build_tri=true;build_tet=true;}

        if(surface_original) for(int i=1;i<=surface_original->particles.array_collection->Size();++i) surface_original->particles.X(i)+=TV(0,(T)15,0);
        for(int i=1;i<=volume_original->particles.array_collection->Size();++i) volume_original->particles.X(i)+=TV(0,(T)15,0);                
        PHYSBAM_ASSERT(!build_tri); // Make baby jesus cry

        if(build_tri || build_tet){
            BOX<TV> grid_domain;
            T thickness=(T).2;
            for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                grid_domain.Enlarge_To_Include_Box(solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Box());

            bool read_in_phi_from_file=false;
            GRID<TV> new_grid((T).02,grid_domain.Thickened(2*thickness));
            ARRAY<T,TV_INT> new_phi(new_grid.Domain_Indices());new_phi.Fill(1e10);
            LEVELSET_IMPLICIT_OBJECT<TV>* implicit=new LEVELSET_IMPLICIT_OBJECT<TV>(new_grid,new_phi);
            
            if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number))){
                read_in_phi_from_file=true;
                FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number),*implicit);}
            else{
                for(CELL_ITERATOR iterator(new_grid);iterator.Valid();iterator.Next()){const TV_INT &cell_index=iterator.Cell_Index();
                    //new_phi(cell_index)=-1;
                    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                        new_phi(cell_index)=min(new_phi(cell_index),solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Extended_Phi(iterator.Location())-thickness);}
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number),*implicit);}
            if(build_tet){
                int x_edge=40;
                TV edges=new_grid.Domain().Edge_Lengths();
                TV edge_cells_float=TV((T)x_edge,(T)x_edge*edges.y/edges.x,(T)x_edge*edges.z/edges.x);
                TV_INT edge_cells(edge_cells_float);
                volume_original->Initialize_Cube_Mesh_And_Particles(GRID<TV>(edge_cells,new_grid.Domain()));
                volume_original->Discard_Tetrahedrons_Outside_Implicit_Surface(*implicit);
                volume_original->Discard_Valence_Zero_Particles_And_Renumber();
                volume_original->Update_Number_Nodes();
                if(read_in_phi_from_file) for(int i=1;i<=volume_original->particles.array_collection->Size();++i) volume_original->particles.X(i)+=TV(0,(T)15,0);                
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),*volume_original);}
            if(build_tri){
                surface_original=DUALCONTOUR_3D<T>::Create_Triangulated_Surface_From_Levelset(implicit->levelset);  
                surface_original=dynamic_cast<TRIANGULATED_SURFACE<T>*>(surface_original->Append_Particles_And_Create_Copy(particles));
                surface_original->Update_Number_Nodes();
                FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tri",test_number),*surface_original);}
            delete implicit;}

        if(surface_original){
            surface=(TRIANGULATED_SURFACE<T>*)surface_original->Append_Particles_And_Create_Copy(particles);
            surface->Update_Number_Nodes();}
        TETRAHEDRALIZED_VOLUME<T>* volume=(TETRAHEDRALIZED_VOLUME<T>*)volume_original->Append_Particles_And_Create_Copy(particles);
        volume->Update_Number_Nodes();
        volume->Initialize_Hierarchy();
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*volume,1,true);
        particles.Store_Velocity();
        if(surface) deformable_body_collection.deformable_geometry.Add_Structure(surface);
        deformable_body_collection.deformable_geometry.Add_Structure(volume);

        if(use_embedding){
            assert(surface);
            ARRAY<TRIPLE<int,int,TV> > bindings; 
            if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number))) FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number),bindings);
            else{
                ARRAY<int> tets;const T tolerance=(T)1e-4;
                ARRAY<int> surface_particles;surface->mesh.elements.Flattened().Get_Unique(surface_particles);
                for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
                    tets.Remove_All();volume->hierarchy->Intersection_List(particles.X(p),tets,tolerance);bool got_bind=false;
                    for(int tt=0;tt<tets.m;tt++){int t=tets(tt);
                        TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(particles.X(p),particles.X.Subset(volume->mesh.elements(t)));
                        if(bary.x>-tolerance && bary.y>-tolerance && bary.z>-tolerance && bary.x+bary.y+bary.z<(T)1+tolerance){bindings.Append(TRIPLE<int,int,TV>(p,t,bary));got_bind=true;break;}}
                    if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(TRIPLE<int,int,TV>(p,0,TV(0,0,0)));}}
                FILE_UTILITIES::Write_To_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/bindings_%d",test_number),bindings);}
            for(int i=0;i<bindings.m;i++){
                if(bindings(i).y==0) continue;
                VECTOR<int,4> nodes=volume->mesh.elements(bindings(i).y);
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,bindings(i).x,nodes,bindings(i).z));
                int soft_bound_particle=particles.array_collection->Add_Element_From_Deletion_List();
                soft_bindings.Add_Binding(VECTOR<int,2>(soft_bound_particle,bindings(i).x),true);}
            surface->Update_Number_Nodes();
            volume->Update_Number_Nodes();}

        //T stiffness=(T)2e5;
        //T damping=(T).01;
        //deformable_object.Add_Force(Create_Finite_Volume(*volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));
    
        T stiffness=(T)1e6;
        T damping=(T)10;
        T restlength_clamp=(T)1e-4;
        T cfl_strain_rate=(T).1;
        bool use_implicit=true,strain_limit=false;

        for(int i=0;i<particles.array_collection->Size();i++){particles.mass(i)=fluids_parameters.density/particles.array_collection->Size()/10;}
        particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
        T linear_stiffness=stiffness,linear_damping=damping;
        LINEAR_SPRINGS<TV>* edge_springs;
        edge_springs=Create_Edge_Springs(*volume,linear_stiffness,linear_damping,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
        edge_springs->Clamp_Restlength(restlength_clamp); 
        solid_body_collection.Add_Force(edge_springs);
        T tet_stiffness=stiffness,tet_damping=damping;
        LINEAR_TET_SPRINGS<T>* tet_springs;
        tet_springs=Create_Tet_Springs(*volume,tet_stiffness,tet_damping,false,(T).1,strain_limit,cfl_strain_rate,true,(T)0,true,use_implicit);
        tet_springs->Clamp_Restlength(restlength_clamp); 
        solid_body_collection.Add_Force(tet_springs);

        // binding the deformable particles to the rigid bodies
        ARRAY<int> particle_array;volume->mesh.elements.Flattened().Get_Unique(particle_array);
        for(int p=0;p<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Unbound_Particles_In_Rigid_Body(solid_body_collection.rigid_body_collection.Rigid_Body(p),particle_array);
    }

    //set up joints
    use_limits=false;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    if(test_number>6){
        arb.poststabilization_projection_iterations=100;
        solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
        height=-8.5;}
    switch(test_number){
        case 1: 
        case 2: 
        case 3: Joints_From_List(ANGLE_JOINT_TYPE);break;
        case 4:
        case 5:
        case 6: Joints_From_List(POINT_JOINT_TYPE);break;
        case 7: 
        case 8:
        case 9:
        case 10:
            Create_Joints_From_Hierarchy(ANGLE_JOINT_TYPE,body_motion.name_to_track_index.Get("Root"));
            if(test_number==7){
                controlled_bones.Append("Spine1");
                uncontrolled_bones.Append("Spine");uncontrolled_bones.Append("Arm");uncontrolled_bones.Append("Clavicle");
                uncontrolled_bones.Append("Hand");uncontrolled_bones.Append("Neck");uncontrolled_bones.Append("Head");}
            else if(test_number==8){
                controlled_bones.Append("Spine");controlled_bones.Append("Neck");controlled_bones.Append("Head");
                uncontrolled_bones.Append("Hand");uncontrolled_bones.Append("Arm");uncontrolled_bones.Append("Clavicle");}
            else if(test_number==9){
                //controlled_bones.Append("Hand");controlled_bones.Append("Arm");controlled_bones.Append("Clavicle");
                uncontrolled_bones.Append("Hand");uncontrolled_bones.Append("Arm");controlled_bones.Append("RUpper_Arm");
                /*uncontrolled_bones.Append("Spine");uncontrolled_bones.Append("Neck");uncontrolled_bones.Append("Head");*/}
            else if(test_number==10){
                controlled_bones.Append("Hand");controlled_bones.Append("Arm");controlled_bones.Append("Clavicle");
                controlled_bones.Append("Spine");controlled_bones.Append("Neck");controlled_bones.Append("Head");}
            if(use_limits) Incorporate_Joint_Limits();
            Initialize_Bone_Hierarchy_For_Human();
            Prune_Joints();
            break;
        case 11:
        case 12:
        case 13:{
            use_clustering=true;
            ARRAY<int> active_clusters;
            PHYSBAM_FATAL_ERROR("need to change the next line to use the actual fluid collision body list");
            //rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,&solids_parameters.collision_body_list);
            Create_Joints_From_Hierarchy(POINT_JOINT_TYPE,body_motion.name_to_track_index.Get("Root"));
            Initialize_Bone_Hierarchy_For_Cluster_Human();
            rigid_bindings.Update_All_Joint_Structures();
            rigid_bindings.Reactivate_Bindings(active_clusters);
            if(use_limits) Incorporate_Joint_Limits();
            Prune_Joints();
            break;}
        case 14: break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up forces
    switch(test_number){
        case 1: 
        case 4: Sideways_Gravity();break;
        case 2: 
        case 5: Constant_Wind();break;
        case 3: 
        case 6: Spatially_Varying_Wind();break;
        case 7:
        case 8:
        case 9: 
        case 10: break;
        case 11:
        case 12:
        case 13: 
        case 14: /*Sideways_Gravity();*/Constant_Wind();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    T friction=(T).8;
    tests.Add_Ground(friction,height,0);
    GRAVITY<TV>* gravity=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,referenced_rigid_particles);
    gravity->Set_Gravity(TV(0,(T)-10,0));
    solid_body_collection.Add_Force(gravity);
    controller->bone_hierarchy=bone_hierarchy;
    
    if(test_number>=7 && test_number<14) controller->root_particle_index=rigid_body_ids(body_motion.name_to_track_index.Get("Root"));
    if((test_number>=7 && test_number<=10) || test_number==12 || test_number==13){
        ARRAY<int> active_clusters;
        PHYSBAM_FATAL_ERROR("need to change the next line to use the actual fluid collision body list");
        //rigid_bindings.Deactivate_And_Return_Clusters(active_clusters,&solids_parameters.collision_body_list);
        controller->Create_All_Clusters(collision_manager);
        rigid_bindings.Reactivate_Bindings(active_clusters);
        RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            LOG::cout<<"Joint ID "<<joint.id_number<<std::endl;
            if(controller->joint_clusters(joint.id_number).x!=0){
                T_CLUSTER& bindings_of_parent_1=*rigid_bindings.reverse_bindings.Get(controller->joint_clusters(joint.id_number).x);        
                if(controller->joint_clusters(joint.id_number).y!=0){
                    T_CLUSTER& bindings_of_parent_2=*rigid_bindings.reverse_bindings.Get(controller->joint_clusters(joint.id_number).y);
                    LOG::cout<<"Cluster "<<controller->joint_clusters(joint.id_number).x<<" and "<<controller->joint_clusters(joint.id_number).y<<std::endl;
                    LOG::cout<<"M: Both Clusters are ["<<bindings_of_parent_1.children<<"]"<<std::endl;
                    LOG::cout<<"M: and ["<<bindings_of_parent_2.children<<"]"<<std::endl;}
                else{
                    LOG::cout<<"Cluster "<<controller->joint_clusters(joint.id_number).x<<std::endl;
                    LOG::cout<<"M: Top Clusters are ["<<bindings_of_parent_1.children<<"]"<<std::endl;}}
            else{
                if(controller->joint_clusters(joint.id_number).y!=0){
                    T_CLUSTER& bindings_of_parent_2=*rigid_bindings.reverse_bindings.Get(controller->joint_clusters(joint.id_number).y);        
                    LOG::cout<<"Cluster "<<controller->joint_clusters(joint.id_number).y<<std::endl;
                    LOG::cout<<"M: Bottom Clusters are ["<<bindings_of_parent_2.children<<"]"<<std::endl;}
                else{
                    LOG::cout<<"M: Clusters are []"<<std::endl;}}
        }
        for(int i=0;i<controller->global_clusters.m;i++){
            T_CLUSTER& bindings_of_parent=*rigid_bindings.reverse_bindings.Get(controller->global_clusters(i));        
            LOG::cout<<"M: Global Cluster ["<<bindings_of_parent.children<<"]";
        }
    }
    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    //deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.template Find_Structure);
    //solids_parameters.collision_body_list.Add_Bodies(solid_body_collection.rigid_body_collection);
    //solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.template Find_Structure);
    
    solids_source=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,new ARRAY<int>());
    solids_source->Set_Gravity(TV(0,0,0));solid_body_collection.Add_Force(solids_source);

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Joint_Between_Parallel
//#####################################################################
void Initialize_Joint_Between_For_Parallel(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,const RIGID_BODY<TV>& parent_cluster,const RIGID_BODY<TV>& child_cluster,ROTATION<TV> rotation,T ratio)
{
    FRAME<TV> J((ratio*child.X()+(1-ratio)*parent.X()),rotation);
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(parent_cluster.Frame()));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(child_cluster.Frame()));
}
//#####################################################################
// Function Initialize_Joint_Between_For_Lines
//#####################################################################
void Initialize_Joint_Between_For_Lines(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,const RIGID_BODY<TV>& parent_cluster,
    const RIGID_BODY<TV>& child_cluster,ROTATION<TV> rotation,FRAME<TV>& inverse_parent,FRAME<TV>& inverse_child)
{
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    TV px1=pf*inverse_parent*TV(0,0,0),px2=pf*inverse_parent*TV(1,0,0);
    TV cx1=cf*inverse_child*TV(0,0,0),cx2=cf*inverse_child*TV(1,0,0);
    SEGMENT_3D<T> parent_segment(px1,px2),child_segment(cx1,cx2);
    VECTOR<T,2> weights;parent_segment.Shortest_Vector_Between_Lines(child_segment,weights);
    FRAME<TV> J((T).5*((px1+weights.x*(px2-px1))+(cx1+weights.y*(cx2-cx1))),rotation);
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(parent_cluster.Frame()));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(child_cluster.Frame()));
}
//#####################################################################
// Function Initialize_Joint_Between
//#####################################################################
void Initialize_Joint_Between(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,const RIGID_BODY<TV>& parent_cluster,const RIGID_BODY<TV>& child_cluster,
    TV up,FRAME<TV> inverse_parent,FRAME<TV> inverse_child,T ratio)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.joint_mesh.Add_Articulation(parent_cluster.particle_index,child_cluster.particle_index,joint);
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    //TV x=(cf.t-pf.t).Normalized();
    //TV y=up.Projected_Orthogonal_To_Unit_Direction(x).Normalized();
    //TV z=TV::Cross_Product(x,y);
    TV x=TV(1,0,0);
    TV y=TV(0,1,0);
    TV z=TV(0,0,1);
    TV parent_vector=(pf.r*inverse_parent.r).Rotate(TV(1,0,0));
    TV child_vector=(cf.r*inverse_child.r).Rotate(TV(1,0,0));
    ROTATION<TV> rotation(MATRIX<T,3>(x,y,z));
    if(body_motion.names(id_to_index(child.particle_index))=="Spine1"){
        joint->impulse_accumulator->weight=1;
        controller->objective.Resize(joint->id_number);
        controller->objective(joint->id_number)=DRAG;
        controller->strain_joints.Append(joint->id_number);}
    else if(body_motion.names(id_to_index(child.particle_index))=="RUpper_Arm" || body_motion.names(id_to_index(child.particle_index))=="LUpper_Arm"){
        joint->impulse_accumulator->weight=1;
        rotation=ROTATION<TV>(MATRIX<T,3>(y,z,x));
        controller->objective.Resize(joint->id_number);
        controller->objective(joint->id_number)=DRAG;}
    else{
        controller->objective.Resize(joint->id_number);
        controller->objective(joint->id_number)=NOFUNCTION;}
    if(abs(TV::Dot_Product(parent_vector.Normalized(),child_vector.Normalized()))>.9) Initialize_Joint_Between_For_Parallel(joint,parent,child,parent_cluster,child_cluster,rotation,ratio);
    else Initialize_Joint_Between_For_Lines(joint,parent,child,parent_cluster,child_cluster,rotation,inverse_parent,inverse_child);
}
//#####################################################################
// Function Human
//#####################################################################
void Human(int frame)
{
    //adding bones
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    referenced_rigid_particles=new ARRAY<int>;
    T density=1000;collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    id_to_index.Resize(int(body_motion.trajectories.m));
    rigid_body_ids.Resize(body_motion.trajectories.m);
    for(int i=0;i<body_motion.trajectories.m;i++){
        std::string body_name=body_motion.names(i);
        T scale=(T).9,length=body_motion.trajectories(i)(1).length,height=scale*length,radius=(T).18;
        RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Cylinder(height,radius,32,32);
        int id=rigid_body.particle_index;assert(id);
        if(id_to_index.Size()<id) id_to_index.Resize(id);
        id_to_index(id)=i;rigid_body_ids(i)=id;
        T volume=(T)pi*radius*radius*height;
        rigid_body.Set_Mass(density*volume);
        rigid_body.Update_Bounding_Box();
        referenced_rigid_particles->Append(rigid_body.particle_index);
        if(i==1){
            T_INERTIA_TENSOR& inertia_tensor=arb.rigid_body_collection.rigid_body_particle.inertia_tensor(rigid_body_ids(i));
            T max=inertia_tensor(1),translation=rigid_body.axis_aligned_bounding_box.min_corner.x;int axis=1;
            if(max<inertia_tensor(2)){max=inertia_tensor(2);axis=2;}
            if(max<inertia_tensor(3)){max=inertia_tensor(3);axis=3;}
            axis=3; //manually set for now
            if(axis==2){translation=rigid_body.axis_aligned_bounding_box.min_corner.y;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,0,-1)));}
            if(axis==3){translation=rigid_body.axis_aligned_bounding_box.min_corner.z;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0)));}
            rigid_base_transform=FRAME<TV>(TV(-1*translation,0,0),ROTATION<TV>())*rigid_base_transform;
            rigid_base_transform.t/=body_motion.trajectories(1)(1).length;}
        FRAME<TV> rigid_base_transform_i=rigid_base_transform;
        rigid_base_transform_i.t*=length;
        rigid_body.Set_Frame(body_motion.trajectories(i)(frame).targeted_transform*rigid_base_transform_i);}
}
//#####################################################################
// Function Human_Cluster
//#####################################################################
void Human_Cluster(int frame)
{
    //adding bones
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    ARRAY<ARRAY<int,RIGID_CLUSTER_CONSTITUENT_ID> > children(cluster_bones.m);
    referenced_rigid_particles=new ARRAY<int>;
    T density=1000;collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    id_to_index.Resize(int(body_motion.trajectories.m+cluster_bones.m));
    rigid_body_ids.Resize(body_motion.trajectories.m+cluster_bones.m);
    for(int i=0;i<body_motion.trajectories.m;i++){
        std::string body_name=body_motion.names(i);
        T scale=(T).9,length=body_motion.trajectories(i)(1).length,height=scale*length,radius=(T).18;
        RIGID_BODY<TV>& rigid_body=tests.Add_Analytic_Cylinder(height,radius);
        int id=rigid_body.particle_index;assert(id);
        if(id_to_index.Size()<id) id_to_index.Resize(id);
        id_to_index(id)=i;rigid_body_ids(i)=id;
        T volume=(T)pi*radius*radius*height;
        rigid_body.Set_Mass(density*volume);
        rigid_body.Update_Bounding_Box();
        //cluster code
        int done=false;
        for(int j=0;j<cluster_bones.m;j++) for(int k=1;k<=cluster_bones(j).m;k++) if(body_name.find(cluster_bones(j)(k))!=std::string::npos){children(j).Append(rigid_body.particle_index);done=true;}
        if(!done) referenced_rigid_particles->Append(rigid_body.particle_index);
        if(i==1){
            T_INERTIA_TENSOR& inertia_tensor=arb.rigid_body_collection.rigid_body_particle.inertia_tensor(rigid_body_ids(i));
            T max=inertia_tensor(1),translation=rigid_body.axis_aligned_bounding_box.min_corner.x;int axis=1;
            if(max<inertia_tensor(2)){max=inertia_tensor(2);axis=2;}
            if(max<inertia_tensor(3)){max=inertia_tensor(3);axis=3;}
            axis=3; //manually set for now
            if(axis==2){translation=rigid_body.axis_aligned_bounding_box.min_corner.y;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,0,-1)));}
            if(axis==3){translation=rigid_body.axis_aligned_bounding_box.min_corner.z;rigid_base_transform=FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(0,1,0)));}
            rigid_base_transform=FRAME<TV>(TV(-1*translation,0,0),ROTATION<TV>())*rigid_base_transform;
            rigid_base_transform.t/=body_motion.trajectories(1)(1).length;}
        FRAME<TV> rigid_base_transform_i=rigid_base_transform;
        rigid_base_transform_i.t*=length;
        rigid_body.Set_Frame(body_motion.trajectories(i)(frame).targeted_transform*rigid_base_transform_i);}
    for(int i=0;i<children.m;i++){
        int cluster_particle=rigid_bindings.Add_Binding(children(i));
        LOG::cout<<" Cluster particle "<<cluster_particle<<" has children "<<children(i)<<std::endl;
        RIGID_BODY<TV>* rigid_body_cluster=&solid_body_collection.rigid_body_collection.Rigid_Body(cluster_particle);
        rigid_body_cluster->Set_Name(cluster_bones(i)(1));
        if(0 && (test_number==11 || test_number==12 || test_number==13)){
            TV rotation_point=solid_body_collection.rigid_body_collection.rigid_body_particle.X(rigid_body_ids(body_motion.name_to_track_index.Get("Root")))*.5+
                solid_body_collection.rigid_body_collection.rigid_body_particle.X(rigid_body_ids(body_motion.name_to_track_index.Get("Spine1")))*.5;
            rigid_body_cluster->Set_Frame(FRAME<TV>(rotation_point,ROTATION<TV>())*FRAME<TV>(TV(),ROTATION<TV>::From_Euler_Angles(45*(T)pi/180,0,0))
                *FRAME<TV>(rotation_point*-1,ROTATION<TV>())*rigid_body_cluster->Frame());}            
        rigid_bindings.Clamp_Particles_To_Embedded_Positions();
        referenced_rigid_particles->Append(rigid_body_cluster->particle_index);
        //update structures
        PHYSBAM_ASSERT(cluster_particle==body_motion.trajectories.m+i);PHYSBAM_ASSERT(rigid_body_cluster->particle_index==int(body_motion.trajectories.m+i));
        rigid_body_ids(cluster_particle)=rigid_body_cluster->particle_index;
        id_to_index(rigid_body_cluster->particle_index)=cluster_particle;
        body_motion.names.Resize(cluster_particle);body_motion.names(cluster_particle)=cluster_bones(i)(1);}
}
//#####################################################################
// Function Cylinder_And_Block
//#####################################################################
void Cylinder_And_Block()
{
    //T friction=(T).8;
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("subdivided_box",1,friction);
    kinematic_body.Is_Kinematic()=1;
    kinematic_body.X()=TV(0,(T)1,0);
    kinematic_body.Update_Bounding_Box();

    //RIGID_BODY<TV>& dynamic_body=tests.Add_Analytic_Box(TV(2,6,.1));
    //RIGID_BODY<TV>& dynamic_body=tests.Add_Analytic_Cylinder(6,1);
    //dynamic_body.X()=TV(0,(T)5,0);
    RIGID_BODY<TV>& dynamic_body=tests.Add_Rigid_Body("cyllink",1,friction);
    dynamic_body.X()=TV(0,(T)5,0);
    dynamic_body.Update_Bounding_Box();
    dynamic_body.Set_Mass(50);

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    collision_manager->hash.Insert(PAIR<int,int>(kinematic_body.particle_index,dynamic_body.particle_index));
    collision_manager->hash.Insert(PAIR<int,int>(dynamic_body.particle_index,kinematic_body.particle_index));
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);
    referenced_rigid_particles->Append(dynamic_body.particle_index);
}
//#####################################################################
// Function Create_Joints_From_Cluster_Hierarchy
//#####################################################################
void Create_Joints_From_Cluster_Hierarchy(int joint_type,int parent)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    for(int i=1;i<=body_motion.bone_hierarchy(parent).m;i++){
        int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(parent)(i));
        int parent_index=rigid_bindings.Get_Parent_Index(rigid_body_ids(parent));
        int child_index=rigid_bindings.Get_Parent_Index(rigid_body_ids(child));
        if(parent_index==0) parent_index=rigid_body_ids(parent);
        if(child_index==0) child_index=rigid_body_ids(child);
        if(parent_index!=child_index){
            JOINT<TV>* joint=0;
            if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
            if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
            joint->impulse_accumulator=arb_impulse_accumulator;
            RIGID_BODY<TV>& parent_cluster_body=solid_body_collection.rigid_body_collection.Rigid_Body(parent_index);
            RIGID_BODY<TV>& child_cluster_body=solid_body_collection.rigid_body_collection.Rigid_Body(child_index);
            RIGID_BODY<TV>& parent_body=arb.rigid_body_list.Rigid_Body(rigid_body_ids(parent));
            RIGID_BODY<TV>& child_body=arb.rigid_body_list.Rigid_Body(rigid_body_ids(child));
            TV up=parent_body.Frame()*TV(0,0,1);
            T child_length=body_motion.trajectories(child)(1).length,parent_length=body_motion.trajectories(parent)(1).length;
            FRAME<TV> rigid_base_transform_child=rigid_base_transform;rigid_base_transform_child.t*=child_length;
            FRAME<TV> rigid_base_transform_parent=rigid_base_transform;rigid_base_transform_parent.t*=parent_length;
            Initialize_Joint_Between(joint,parent_body,child_body,parent_cluster_body,child_cluster_body,up,rigid_base_transform_parent.Inverse(),rigid_base_transform_child.Inverse(),parent_length/(child_length+parent_length));
            collision_manager->hash.Insert(PAIR<int,int>(child_body.particle_index,parent_cluster_body.particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(parent_body.particle_index,child_cluster_body.particle_index));
            JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
            joint_function->Set_k_p(1000);joint_function->Set_Target_Angle(ROTATION<TV>());}
        Create_Joints_From_Cluster_Hierarchy(joint_type,child);}
}
//#####################################################################
// Function Create_Joints_From_Hierarchy
//#####################################################################
void Create_Joints_From_Hierarchy(int joint_type,int parent)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=1;i<=body_motion.bone_hierarchy(parent).m;i++){
        std::string body_name=body_motion.names(parent);
        int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(parent)(i));
        JOINT<TV>* joint=0;
        if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
        if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
        joint->impulse_accumulator=arb_impulse_accumulator;
        RIGID_BODY<TV> &parent_body=arb.rigid_body_collection.Rigid_Body(rigid_body_ids(parent)),&child_body=arb.rigid_body_collection.Rigid_Body(rigid_body_ids(child));
        TV up=parent_body.Frame()*TV(0,0,1);
        T child_length=body_motion.trajectories(child)(1).length,parent_length=body_motion.trajectories(parent)(1).length;
        FRAME<TV> rigid_base_transform_child=rigid_base_transform;rigid_base_transform_child.t*=child_length;
        FRAME<TV> rigid_base_transform_parent=rigid_base_transform;rigid_base_transform_parent.t*=parent_length;
        Initialize_Joint_Between(joint,parent_body,child_body,parent_body,child_body,up,rigid_base_transform_parent.Inverse(),rigid_base_transform_child.Inverse(),parent_length/(child_length+parent_length));
        collision_manager->hash.Insert(PAIR<int,int>(rigid_body_ids(child),rigid_body_ids(parent)));
        collision_manager->hash.Insert(PAIR<int,int>(rigid_body_ids(parent),rigid_body_ids(child)));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(1000);joint_function->Set_Target_Angle(ROTATION<TV>());
        Create_Joints_From_Hierarchy(joint_type,child);}
}
//#####################################################################
// Function Initialize_Bone_Hierarchy_For_Human
//#####################################################################
void Initialize_Bone_Hierarchy_For_Human()
{
    bone_hierarchy.Resize(body_motion.bone_hierarchy.m);
    for(int i=0;i<body_motion.bone_hierarchy.m;i++){
        std::string body_name=body_motion.names(i);
        int parent=i,parent_index=rigid_body_ids(parent);
        for(int j=1;j<=body_motion.bone_hierarchy(i).m;j++){
            int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(i)(j)),child_index=rigid_body_ids(child);
            bone_hierarchy(parent_index).Append(child_index);}}
}
//#####################################################################
// Function Initialize_Bone_Hierarchy_For_Cluster_Human
//#####################################################################
void Initialize_Bone_Hierarchy_For_Cluster_Human()
{
    RIGID_BODY_CLUSTER_BINDINGS<TV>& rigid_bindings=solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings;
    bone_hierarchy.Resize(body_motion.bone_hierarchy.m+cluster_bones.m);
    for(int i=0;i<body_motion.bone_hierarchy.m;i++){
        int parent=i,parent_index=rigid_body_ids(parent);
        int cluster_parent_index=rigid_bindings.Get_Parent_Index(parent_index);
        if(!cluster_parent_index) cluster_parent_index=parent_index; 
        for(int j=1;j<=body_motion.bone_hierarchy(i).m;j++){
            int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(i)(j)),child_index=rigid_body_ids(child);
            int cluster_child_index=rigid_bindings.Get_Parent_Index(child_index);
            if(!cluster_child_index) cluster_child_index=child_index; 
            if(cluster_parent_index!=cluster_child_index) bone_hierarchy(cluster_parent_index).Append(cluster_child_index);}}
}
//#####################################################################
// Function Prune_Joints
//#####################################################################
void Prune_Joints()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
        bool done=false;RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        for(int i=0;i<controlled_bones.m;i++) if(child_name.find(controlled_bones(i))!=std::string::npos){
            joint.joint_function->active=false;joint.global_post_stabilization=false;
            for(int i=0;i<T_SPIN::dimension;i++) joint.control_dof(i)=true;
            done=true;break;}
        if(done) continue;
        for(int i=0;i<uncontrolled_bones.m;i++) if(child_name.find(uncontrolled_bones(i))!=std::string::npos){
            /*delete joint.joint_function;joint.joint_function=0;*/joint.joint_function->active=true;
            done=true;break;}
        if(done) continue;
        if(kinematic_motion){delete joint.joint_function;joint.joint_function=0;child->Is_Kinematic()=true;parent->Is_Kinematic()=true;}
        else joint.joint_function->active=true;}
    //for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
    //    if(joint.id_number>=JOINT_ID(7)&&joint.id_number<=JOINT_ID(14)){delete joint.joint_function;arb.joint_mesh.Remove_Articulation(joint.id_number);}}
}
//#####################################################################
// Function Incorporate_Joint_Limits
//#####################################################################
void Incorporate_Joint_Limits()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){
        POINT_JOINT<TV>* joint=dynamic_cast<POINT_JOINT<TV>*>(arb.joint_mesh.joints(i));
        if(!joint) continue;
        RIGID_BODY<TV>* parent=arb.Parent(joint->id_number),*child=arb.Child(joint->id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        //spine joints
        if(parent_name.find("Spine")!=std::string::npos && parent_name.find("Spine1")==std::string::npos && child_name.find("Spine")!=std::string::npos){
            joint->Use_Twist_Constraint(-15*(T)pi/180,15*(T)pi/180);joint->Use_Phi_Constraint(-15*(T)pi/180,15*(T)pi/180);joint->Use_Theta_Constraint(-15*(T)pi/180,15*(T)pi/180);}
        //hip joints
        if(parent_name.find("Spine1")!=std::string::npos){
            joint->Use_Twist_Constraint(-45*(T)pi/180,45*(T)pi/180);joint->Use_Phi_Constraint(-45*(T)pi/180,45*(T)pi/180);joint->Use_Theta_Constraint(-45*(T)pi/180,90*(T)pi/180);}
        //neck joints
        if(child_name.find("Neck")!=std::string::npos || child_name.find("Head")!=std::string::npos){
            joint->Use_Twist_Constraint(-45*(T)pi/180,45*(T)pi/180);joint->Use_Phi_Constraint(-45*(T)pi/180,45*(T)pi/180);joint->Use_Theta_Constraint(-45*(T)pi/180,45*(T)pi/180);}
        //clavicle
        if(child_name.find("Clavicle")!=std::string::npos){
            joint->Use_Twist_Constraint(-15*(T)pi/180,15*(T)pi/180);joint->Use_Phi_Constraint(-30*(T)pi/180,30*(T)pi/180);joint->Use_Theta_Constraint(-30*(T)pi/180,30*(T)pi/180);}
        //arms
        if(child_name.find("Arm")!=std::string::npos && child_name.find("Upper")!=std::string::npos){
            joint->Use_Twist_Constraint(0*(T)pi/180,0*(T)pi/180);joint->Use_Phi_Constraint(-60*(T)pi/180,150*(T)pi/180);joint->Use_Theta_Constraint(-60*(T)pi/180,150*(T)pi/180);}
        if(child_name.find("Arm")!=std::string::npos && child_name.find("Lower")!=std::string::npos){
            joint->Use_Twist_Constraint(0*(T)pi/180,90*(T)pi/180);joint->Use_Phi_Constraint(0*(T)pi/180,150*(T)pi/180);joint->Use_Theta_Constraint(0*(T)pi/180,0*(T)pi/180);}
        //hands
        if(child_name.find("Hand")!=std::string::npos){
            joint->Use_Twist_Constraint(0*(T)pi/180,0*(T)pi/180);joint->Use_Phi_Constraint(-30*(T)pi/180,30*(T)pi/180);joint->Use_Theta_Constraint(-90*(T)pi/180,90*(T)pi/180);}
        //TODO: Should we limit lower body?
    }
}
//#####################################################################
// Function Root_Joint
//#####################################################################
void Root_Joint()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
        RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        if(child_name.find("Spine1")!=std::string::npos){
            joint.use_joint_function=false;joint.global_post_stabilization=false;
            for(int i=0;i<T_SPIN::dimension;i++) joint.control_dof(i)=true;}
        else if(child_name.find("Spine")!=std::string::npos || child_name.find("Arm")!=std::string::npos || child_name.find("Clavicle")!=std::string::npos ||
            child_name.find("Hand")!=std::string::npos || child_name.find("Neck")!=std::string::npos || child_name.find("Head")!=std::string::npos){
            /*delete joint.joint_function;joint.joint_function=0;*/joint.use_pd=true;}
        else if(kinematic_motion){
            delete joint.joint_function;joint.joint_function=0;child->Is_Kinematic()=true;parent->Is_Kinematic()=true;}
        else joint.use_pd=true;}
}
//#####################################################################
// Function Spine_Joints
//#####################################################################
void Spine_Joints()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
        RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        if(child_name.find("Spine")!=std::string::npos || child_name.find("Neck")!=std::string::npos || child_name.find("Head")!=std::string::npos){
            joint.use_joint_function=false;joint.global_post_stabilization=false;
            for(int i=0;i<T_SPIN::dimension;i++) joint.control_dof(i)=true;}
        else if(child_name.find("Arm")!=std::string::npos || child_name.find("Clavicle")!=std::string::npos || child_name.find("Hand")!=std::string::npos){
            delete joint.joint_function;joint.joint_function=0;}
        else if(kinematic_motion){delete joint.joint_function;joint.joint_function=0;child->Is_Kinematic()=true;parent->Is_Kinematic()=true;}
        else joint.use_pd=true;}
}
//#####################################################################
// Function Arm_Joints
//#####################################################################
void Arm_Joints()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
        RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        if(child_name.find("Arm")!=std::string::npos || child_name.find("Clavicle")!=std::string::npos || child_name.find("Hand")!=std::string::npos){
            joint.use_joint_function=false;joint.global_post_stabilization=false;
            for(int i=0;i<T_SPIN::dimension;i++) joint.control_dof(i)=true;}
        else if(child_name.find("Spine")!=std::string::npos || child_name.find("Neck")!=std::string::npos || child_name.find("Head")!=std::string::npos){
            delete joint.joint_function;joint.joint_function=0;}
        else if(kinematic_motion){delete joint.joint_function;joint.joint_function=0;child->Is_Kinematic()=true;parent->Is_Kinematic()=true;}
        else joint.use_pd=true;}
}
//#####################################################################
// Function Upper_Joints
//#####################################################################
void Upper_Joints()
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
        RIGID_BODY<TV>* parent=arb.Parent(joint.id_number),*child=arb.Child(joint.id_number);
        int parent_index=id_to_index(parent->particle_index),child_index=id_to_index(child->particle_index);
        std::string parent_name=body_motion.names(parent_index),child_name=body_motion.names(child_index);
        if(child_name.find("Spine")!=std::string::npos || child_name.find("Arm")!=std::string::npos || child_name.find("Clavicle")!=std::string::npos ||
            child_name.find("Hand")!=std::string::npos || child_name.find("Neck")!=std::string::npos || child_name.find("Head")!=std::string::npos){
            joint.use_joint_function=false;joint.global_post_stabilization=false;
            for(int i=0;i<T_SPIN::dimension;i++) joint.control_dof(i)=true;}
        else if(kinematic_motion){delete joint.joint_function;joint.joint_function=0;child->Is_Kinematic()=true;parent->Is_Kinematic()=true;}
        else joint.use_pd=true;}
}
//#####################################################################
// Function Joints_From_List
//#####################################################################
void Joints_From_List(int joint_type)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    bone_hierarchy.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size());
    controller->root_particle_index=1;
    controller->objective.Resize(JOINT_ID(Value(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size())-1));
    for(int id(1);id<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        JOINT<TV>* joint=0;
        if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
        if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
        for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
        joint->impulse_accumulator=arb_impulse_accumulator;
        joint->global_post_stabilization=false;
        bone_hierarchy(id).Append(id+1);
        arb.joint_mesh.Add_Articulation(id,id+1,joint);
        controller->objective(joint->id_number)=DRAG;
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>());
        //joint->Set_Joint_To_Child_Frame(FRAME<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.X(id).t-solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Frame(),solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Rotation().Inverse()));
        joint->Set_Joint_To_Child_Frame(solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Frame().Inverse()*solid_body_collection.rigid_body_collection.Rigid_Body(id).Frame());
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->active=false;joint_function->Set_k_p(1000);}
}
//#####################################################################
// Function Octosquid
//#####################################################################
void Octosquid()
{
    int num_legs=5;
    int num_segments_per_leg=1;
    bone_hierarchy.Resize(solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size()+num_legs*(num_segments_per_leg*2)-num_legs+1);
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.Set_Iterative_Tolerance((T)1e-4);
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    arb.Set_Poststabilization_Iterations(5);
    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    arb.constrain_pd_directions=true;

    T radius=(T).9,length=(T)2.5,width=(T).25;
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("sphere",radius,friction);
    kinematic_body.Set_Mass((T)40*radius*radius*fluids_parameters.density);
    kinematic_body.X()=TV(0,(T)15.75,0);
    kinematic_body.Update_Bounding_Box();
    octosquid_body=&kinematic_body;
    if(!use_deformable) Add_Volumetric_Body_To_Fluid_Simulation(kinematic_body);
    controller->not_affected_by_fluid.Set(kinematic_body.particle_index);    
    controller->root_particle_index=kinematic_body.particle_index;

    collision_manager=new RIGID_BODY_COLLISION_MANAGER_HASH;
    referenced_rigid_particles=new ARRAY<int>;
    referenced_rigid_particles->Append(kinematic_body.particle_index);

    for(int i=0;i<num_legs;i++){
        T offset=radius-(T).1767767;
        //offset+=.2;
        offset-=(T).45;
        RIGID_BODY<TV>* prev_link=&kinematic_body;
        TV direction=ROTATION<TV>::From_Euler_Angles(0,2*(T)pi/num_legs*(i-1),0).Rotate(TV(0,0,1));
        for(int j=0;j<num_segments_per_leg;j++){
            RIGID_BODY<TV>& tail_link=tests.Add_Analytic_Cylinder(length,width);
            tail_link.Set_Frame(kinematic_body.Frame());tail_link.X().y-=radius/2+(T).3;tail_link.X()+=(offset+length)*direction;offset+=length+(T).1767767;
            tail_link.Rotation()=ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction);
            FRAME<TV> J;
            RIGID_BODY<TV>* joint_cover=0;
            if(j==1){
                tail_link.X()-=direction*(T).1;
                FRAME<TV> parent_frame=prev_link->Frame();parent_frame.t.y-=(radius/2+(T).3);
                J=FRAME<TV>((tail_link.X()*(T).4+parent_frame.t*(T).6),ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction));}
            else{
                J=FRAME<TV>((tail_link.X()*(T).5+prev_link->X()*(T).5),ROTATION<TV>::From_Rotated_Vector(TV(0,0,1),direction));
                if(!use_deformable){
                    joint_cover=&tests.Add_Rigid_Body("sphere",(T).1767767,friction); // sqrt(2)*.25/2
                    joint_cover->X()=J.t;}}

            static const T volume=(T).15*length;
            tail_link.Set_Mass((T).4*fluids_parameters.density*volume);
            if(!use_deformable){
                Add_Volumetric_Body_To_Fluid_Simulation(tail_link);
                controller->not_affected_by_fluid.Set(tail_link.particle_index);
                if(j!=1){Add_Volumetric_Body_To_Fluid_Simulation(*joint_cover);controller->not_affected_by_fluid.Set(joint_cover->particle_index);}}
            bone_hierarchy(prev_link->particle_index).Append(tail_link.particle_index);
            
            referenced_rigid_particles->Append(tail_link.particle_index);
            collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,tail_link.particle_index));
            collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,prev_link->particle_index));
            if(j!=1 && !use_deformable){
                collision_manager->hash.Insert(PAIR<int,int>(prev_link->particle_index,joint_cover->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,prev_link->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(tail_link.particle_index,joint_cover->particle_index));
                collision_manager->hash.Insert(PAIR<int,int>(joint_cover->particle_index,tail_link.particle_index));}

            ANGLE_JOINT<TV> *joint=new ANGLE_JOINT<TV>;
            { // Create the PD-controlled joint that drives the up-down motion
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
                joint->impulse_accumulator=arb_impulse_accumulator;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,tail_link.particle_index,joint);
                joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                joint->Set_Child_To_Joint_Frame(J.Inverse_Times(tail_link.Frame()));
                controller->objective.Resize(joint->id_number);
                controller->objective(joint->id_number)=DRAG;
                for(int i=0;i<T_SPIN::dimension;i++) joint->control_dof(i)=true;
                JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
                joint_function->Set_k_p(750);joint_function->Set_Target_Angle(ROTATION<TV>());
                joint->global_post_stabilization=false;joint->joint_function->active=false;
                T angle=ROTATION<TV>::From_Rotated_Vector(direction,TV(0,-1,0)).Angle();
                //joint->Use_Rotation_Constraint(-abs(angle),abs(angle));
                if(angle<0) joint->Set_Angle_Constraints(true,angle,(T)0);
                else joint->Set_Angle_Constraints(true,(T)0,angle);
    //            if(constrain_joints) driven_joint->Use_Rotation_Constraint((T)-pi/4,(T)pi/4); // TODO(jontg): broken in 2D...
            }
            if(j!=1 && !use_deformable) { // Create the joint that keeps our protective sphere in place
                POINT_JOINT<TV> *cover_joint=new POINT_JOINT<TV>;
                ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*cover_joint,arb);
                cover_joint->impulse_accumulator=arb_impulse_accumulator;cover_joint->joint_function->active=true;
                arb.joint_mesh.Add_Articulation(prev_link->particle_index,joint_cover->particle_index,cover_joint);
                cover_joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(prev_link->Frame()));
                cover_joint->Set_Child_To_Joint_Frame(J.Inverse_Times(joint_cover->Frame()));
                JOINT_FUNCTION<TV>* cover_joint_function=arb.Create_Joint_Function(cover_joint->id_number);
                cover_joint_function->Set_k_p(750);cover_joint_function->Set_Target_Angle(ROTATION<TV>());
            }
            prev_link=&tail_link;
        }
    }
}
//#####################################################################
// Function Sideways_Gravity
//#####################################################################
void Sideways_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* gravity=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,referenced_rigid_particles);
    gravity->Set_Gravity(TV(0,0,(T)-10));
    solid_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Constant_Wind
//#####################################################################
void Constant_Wind()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(use_deformable){
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(tetrahedralized_volume,solid_body_collection.rigid_body_collection);
        if(test_number==14){drag->Use_Constant_Wind(0,TV(0,0,0));drag->Set_Wind_Density((T).001);}
        else{drag->Use_Constant_Wind(0,TV(0,0,(T)-1000));drag->Set_Wind_Density((T).001);}
        solid_body_collection.Add_Force(drag);
        return;}
    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        //if(use_clustering && Value(id)<=body_motion.trajectories.m) continue;
        if(body_motion.trajectories.m && Value(id)>body_motion.trajectories.m) continue;
        WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(solid_body_collection.rigid_body_collection.Rigid_Body(id),
            static_cast<PARTICLES<TV>&>(solid_body_collection.deformable_body_collection.particles));
        //drag->Use_Constant_Wind(0,TV(0,0,(T)-50));drag->Set_Wind_Density((T)1);
        if(test_number==14){drag->Use_Constant_Wind(0,TV(0,0,0));drag->Set_Wind_Density((T).1);}
        else{drag->Use_Constant_Wind(0,TV(0,0,(T)-100));drag->Set_Wind_Density((T).1);}
        //drag->Use_Constant_Wind(0,TV(0,0,(T)-1000));drag->Set_Wind_Density((T).001);
        solid_body_collection.Add_Force(drag);}
}
//#####################################################################
// Function Spatially_Varying_Wind
//#####################################################################
void Spatially_Varying_Wind()
{
    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(solid_body_collection.rigid_body_collection.Rigid_Body(id),
            static_cast<PARTICLES<TV>&>(solid_body_collection.deformable_body_collection.particles));
        drag->Set_Wind_Density((T)-10);solid_body_collection.Add_Force(drag);}
    Setup_Spatially_Varying_Wind(0);
}
//#####################################################################
};
}
#endif
