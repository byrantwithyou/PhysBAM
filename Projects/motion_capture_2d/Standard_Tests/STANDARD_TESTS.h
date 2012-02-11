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
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform/READ_WRITE_GRID.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_ARRAYS.h>
#include <PhysBAM_Tools/Read_Write/Grids_Uniform_Arrays/READ_WRITE_FACE_ARRAYS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISION_MANAGER_HASH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Fragments/PARTICLE_CONNECTIVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Fluids/PhysBAM_Incompressible/Standard_Tests/SMOKE_STANDARD_TESTS_2D.h>
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
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,2> > >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef GRID<TV> T_GRID;
    typedef typename INTERPOLATION_POLICY<T_GRID>::FACE_LOOKUP T_FACE_LOOKUP;
    typedef typename INTERPOLATION_POLICY<T_GRID>::LINEAR_INTERPOLATION_SCALAR T_LINEAR_INTERPOLATION_SCALAR;
    typedef typename T_GRID::CELL_ITERATOR CELL_ITERATOR;
    typedef VECTOR<int,2> TV_INT;
    typedef typename TV::SPIN T_SPIN;

    typedef typename RIGID_BODY_CLUSTER_BINDINGS<TV>::CLUSTER T_CLUSTER;
public:
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
    bool kinematic_motion,use_clustering,use_limits,use_deformable;
    T initial_angle;
    //spatially varying wind
    ARRAY<T,FACE_INDEX<TV::dimension> > face_velocities;
    ARRAY<TV,VECTOR<int,2> > v_array;
    ARRAY<T,VECTOR<int,2> > p_array,d_array;
    T_GRID grid;
    std::string input_directory;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::fluid_collection;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::frame_rate;using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),driver(0),tests(*this,solid_body_collection),referenced_rigid_particles(0),
        collision_manager(0),kinematic_motion(true),use_clustering(false),use_limits(false),use_deformable(false),initial_angle(0)
    {
    }

    ~STANDARD_TESTS()
    {}
    
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_PD_Targets(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return false;}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}

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
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.rigid_body_evolution_parameters.rigid_geometry_evolution_parameters.use_kinematic_keyframes=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
    solids_parameters.use_rigid_deformable_contact=true;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    //solids_parameters.perform_collision_body_collisions=true;
    solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
    solids_parameters.cfl=5;

    initial_angle=parse_args->Get_Double_Value("-initial_angle");
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    input_directory="/n/field/disk2/mlentine/PhysBAM/Projects/fluids_3d/Standard_Tests_Smoke/Test_1__Resolution_50_50_75/";
    frame_rate=30;

    switch(test_number){
        case 1: last_frame=360;break;
        case 2: last_frame=360;break;
        case 3: last_frame=360;break;
        case 4: last_frame=360;break;
        case 5: last_frame=360;break;
        case 6: last_frame=360;break;
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
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/pressure."+f,p_array);
    FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/density."+f,d_array);
    if(frame==0){
        FILE_UTILITIES::Read_From_File(stream_type,input_directory+"/grid",grid);
        grid.Initialize(face_velocities.Component(1).counts.x,face_velocities.Component(1).counts.y,-5.0,5.0,0.0,10.0);
        FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/grid",grid);}
    for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){TV_INT cell_index=iterator.Cell_Index();
        v_array(cell_index)=TV((face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x+1,cell_index.y)))/(T)2,
            (face_velocities(1,cell_index)+face_velocities(1,TV_INT(cell_index.x,cell_index.y+1)))/(T)2);
    }
    for(int i=0;WIND_DRAG<TV>* drag=solid_body_collection.template Find_Force<WIND_DRAG<TV>*>(i);i++){
        drag->Use_Spatially_Varying_Wind((T)0,grid,v_array);drag->Set_Wind_Pressure(p_array);}
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
void Write_Output_Files(const int frame) const PHYSBAM_OVERRIDE
{
    BASE::Write_Output_Files(frame);
    std::string f=STRING_UTILITIES::string_sprintf("%d",frame);
    LOG::cout<<"Writing "<<output_directory+"/mac_velocities."+f<<std::endl;
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/density."+f,d_array);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/mac_velocities."+f,face_velocities);
    FILE_UTILITIES::Write_To_File(stream_type,output_directory+"/pressure."+f,p_array);
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
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(controller && !controller->hypothetical_step) controller->Update_Position_Based_State(fluid_collection.incompressible_fluid_collection.face_velocities,dt,time);
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
    RIGID_BODY<TV>* kinematic_body=&solid_body_collection.rigid_body_collection.Rigid_Body(id);
    frame=kinematic_body->Frame();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    controller=new SEARCH_CONTROLLER<T_GRID>(solid_body_collection,driver);
    controller->solve_minimization=false;
    controller->max_iterations=1;
    controller->minimize=true;
    controller->use_projection=true;
    controller->drag_direction=TV(1,0);
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
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
#if 0
    if(use_deformable){
        TETRAHEDRALIZED_VOLUME<T>& volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);        
        if(FILE_UTILITIES::File_Exists(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number)))
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),volume);
        else{
            RANGE<TV> grid_domain;
            for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++) grid_domain.Enlarge_To_Include_Box(solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Box());
            GRID<TV> new_grid((T).02,grid_domain);
            ARRAY<T,VECTOR<int,3> > new_phi(new_grid);new_phi.Fill(FLT_MAX);
            T thickness=(T).1;
            for(CELL_ITERATOR iterator(new_grid);iterator.Valid();iterator.Next()){const TV_INT &cell_index=iterator.Cell_Index();
                //new_phi(cell_index)=-1;
                for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++)
                    new_phi(cell_index)=min(new_phi(cell_index),solid_body_collection.rigid_body_collection.Rigid_Body(id).implicit_object->Extended_Phi(iterator.Location())-thickness);}
            
            LEVELSET_IMPLICIT_OBJECT<TV> implicit(new_grid,new_phi);
            FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.phi",test_number),implicit);
            int x_edge=5;
            TV edges=new_grid.Domain().Edge_Lengths();
            TV edge_cells_float=TV((T)x_edge,(T)x_edge*edges.y/edges.x,(T)x_edge*edges.z/edges.x);
            TV_INT edge_cells(edge_cells_float);
            volume.Initialize_Cube_Mesh_And_Particles(GRID_3D<T>(edge_cells,new_grid.Domain()));
            volume.Discard_Tetrahedrons_Outside_Implicit_Surface(implicit);
            volume.Discard_Valence_Zero_Particles_And_Renumber();
            volume.Update_Number_Nodes();
            FILE_UTILITIES::Write_To_File<float>(data_directory+STRING_UTILITIES::string_sprintf("/joint_levelset_%d.tet",test_number),volume);}

        T stiffness=(T)2e4;
        T damping=(T).01;
        volume.Initialize_Hierarchy();
        particles.Store_Velocity();
        volume.Set_Density(1);
        volume.Set_Mass_Of_Particles(true,true);
        deformable_body_collection.deformable_geometry.Add_Structure(&volume);
        solid_body_collection.Add_Force(Create_Finite_Volume(volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));
        
        // binding the deformable particles to the rigid bodies
        for(int p=0;p<solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();p++) tests.Bind_Particles_In_Rigid_Body(solid_body_collection.rigid_body_collection.Rigid_Body(p));
    }
#endif

    //set up joints
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    switch(test_number){
        case 1: 
        case 2: 
        case 3: Joints_From_List(ANGLE_JOINT_TYPE);break;
        case 4:
        case 5:
        case 6: Joints_From_List(POINT_JOINT_TYPE);break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    //set up forces
    switch(test_number){
        case 1: 
        case 4: Sideways_Gravity();break;
        case 2: 
        case 5: Constant_Wind();break;
        case 3: 
        case 6: Spatially_Varying_Wind();break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    T friction=(T).8;
    tests.Add_Ground(friction,height,0);
    GRAVITY<TV>* gravity=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,referenced_rigid_particles);
    gravity->Set_Gravity(TV(0,(T)-10));
    solid_body_collection.Add_Force(gravity);
    controller->bone_hierarchy=bone_hierarchy;
    
    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

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
    SEGMENT_2D<T> parent_segment(px1,px2),child_segment(cx1,cx2);
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
    }//    controller->strain_joints.Append(joint->id_number);}
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
// Function Cylinder_And_Block
//#####################################################################
void Cylinder_And_Block()
{
    T friction=(T)0;
    RIGID_BODY<TV>& kinematic_body=tests.Add_Rigid_Body("square_refined",(T).95,friction);
    kinematic_body.Is_Kinematic()=1;
    kinematic_body.X()=TV(0,(T)1);
    kinematic_body.Update_Bounding_Box();

    T width=0;
    //switch(test_number){
        //case 1: width=(T)2; break;
        //case 2: width=(T).25; break;
        //default: width=(T)2;}
    width=.25;
    RIGID_BODY<TV>& dynamic_body=tests.Add_Analytic_Box(TV(width,(T)4));
    dynamic_body.Rotation()=ROTATION<TV>::From_Angle((T)pi*initial_angle/180);
    dynamic_body.X()=TV(0,(T)1)+(T)3*dynamic_body.Rotation().Rotate(TV(0,(T)1));
    dynamic_body.Update_Bounding_Box();

    static const T volume=width*(T)4;
    dynamic_body.Set_Mass(fluids_parameters.density*volume);

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
    for(int i=0;i<body_motion.bone_hierarchy(parent).m;i++){
        int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(parent)(i));
        int parent_index=rigid_bindings.Get_Parent_Index(arb.rigid_body_list(rigid_body_ids(parent))->particle_index);
        int child_index=rigid_bindings.Get_Parent_Index(arb.rigid_body_list(rigid_body_ids(child))->particle_index);
        if(parent_index==0) parent_index=arb.rigid_body_list(rigid_body_ids(parent))->particle_index;
        if(child_index==0) child_index=arb.rigid_body_list(rigid_body_ids(child))->particle_index;
        if(parent_index!=child_index){
            JOINT<TV>* joint=0;
            if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
            if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
            ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
            joint->impulse_accumulator=arb_impulse_accumulator;
            RIGID_BODY<TV>& parent_cluster_body=solid_body_collection.rigid_body_collection.Rigid_Body(parent_index);
            RIGID_BODY<TV>& child_cluster_body=solid_body_collection.rigid_body_collection.Rigid_Body(child_index);
            RIGID_BODY<TV>& parent_body=*arb.rigid_body_list(rigid_body_ids(parent));
            RIGID_BODY<TV>& child_body=*arb.rigid_body_list(rigid_body_ids(child));
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
    for(int i=0;i<body_motion.bone_hierarchy(parent).m;i++){
        std::string body_name=body_motion.names(parent);
        int child=body_motion.name_to_track_index.Get(body_motion.bone_hierarchy(parent)(i));
        JOINT<TV>* joint=0;
        if(joint_type==ANGLE_JOINT_TYPE) joint=new ANGLE_JOINT<TV>;
        if(joint_type==POINT_JOINT_TYPE) joint=new POINT_JOINT<TV>;
        ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>* arb_impulse_accumulator=new ARTICULATED_RIGID_BODY_IMPULSE_ACCUMULATOR<TV>(*joint,arb);
        joint->impulse_accumulator=arb_impulse_accumulator;
        RIGID_BODY<TV> &parent_body=*arb.rigid_body_list(rigid_body_ids(parent)),&child_body=*arb.rigid_body_list(rigid_body_ids(child));
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
        bone_hierarchy(id).Append(id+1);
        arb.joint_mesh.Add_Articulation(id,id+1,joint);
        controller->objective(joint->id_number)=DRAG;
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>());
        //joint->Set_Joint_To_Child_Frame(FRAME<TV>(solid_body_collection.rigid_body_collection.rigid_body_particle.X(id).t-solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Frame(),solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Rotation().Inverse()));
        joint->Set_Joint_To_Child_Frame(solid_body_collection.rigid_body_collection.Rigid_Body(id+1).Frame().Inverse()*solid_body_collection.rigid_body_collection.Rigid_Body(id).Frame());
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(1000);}
}
//#####################################################################
// Function Sideways_Gravity
//#####################################################################
void Sideways_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* gravity=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,(ARRAY<int>*)0,referenced_rigid_particles);
    gravity->Set_Gravity(TV((T)10,0));
    solid_body_collection.Add_Force(gravity);
}
//#####################################################################
// Function Constant_Wind
//#####################################################################
void Constant_Wind()
{
    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        //if(use_clustering && Value(id)<=body_motion.trajectories.m) continue;
        if(body_motion.trajectories.m && Value(id)>body_motion.trajectories.m) continue;
        WIND_DRAG<TV>* drag=new WIND_DRAG<TV>(solid_body_collection.rigid_body_collection.Rigid_Body(id),
            static_cast<PARTICLES<TV>&>(solid_body_collection.deformable_body_collection.particles));
        drag->Use_Constant_Wind(0,TV((T)100,0));drag->Set_Wind_Density((T).1);
        solid_body_collection.Add_Force(drag);}
}
//#####################################################################
// Function Spatially_Varying_Wind
//#####################################################################
void Spatially_Varying_Wind()
{
    for(int id(1);id<=solid_body_collection.rigid_body_collection.rigid_body_particle.array_collection->Size();id++){
        WIND_DRAG<TV>* drag=new WIND_DRAG<TV>(solid_body_collection.rigid_body_collection.Rigid_Body(id),
            static_cast<PARTICLES<TV>&>(solid_body_collection.deformable_body_collection.particles));
        drag->Set_Wind_Density((T)-10);solid_body_collection.Add_Force(drag);}
    Setup_Spatially_Varying_Wind(0);
}
//#####################################################################
};
}
#endif
