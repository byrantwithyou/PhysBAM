//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Jon Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//  1. Kinematic cube suspending a rigid cube by a linear spring
//  2. Kinematic particle attached to a rigid body
//  3. Rigid Spring
//  4. Cloth panel created from rigid bodies with springs
//  5. Two springs around a single particle
//  6. Impulse down a rigid body chain
//  7. Cloth panel
//  8. Impulse down a particle chain
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Arrays_Computations/ARRAY_COPY.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_EXAMPLE<VECTOR<T_input,3> >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    int parameter;
    int kinematic_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;
    INTERPOLATION_CURVE<T,TV> particle_curve;
    int grid_m;
    int grid_n;

    T ground_angle_rad,mu;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    T aspect_ratio;
    int number_side_panels;
    T stiffness_multiplier;
    T damping_multiplier;
    T bending_stiffness_multiplier;
    T bending_damping_multiplier;
    ROTATION<TV> torus_rotation;
    TV torus_min,torus_max;
    T period,start_time;
    RIGID_BODY<TV>* ground;
    int number_of_joints;
    T wavelength,height_amplitude,bend_amplitude,joint_separation;
    ROTATION<TV> initial_orientation;
    ARRAY<PAIR<TETRAHEDRALIZED_VOLUME<T>*,T> > structure_clamp_time;
    ARRAY<PAIR<RIGID_BODY<TV>*,T> > rigid_body_clamp_time;
    T maximum_fall_speed;
    ARRAY<VECTOR<int,3> > surface_elements;
    ARRAY<int> surface_particles;
    int subsamples;
    RANDOM_NUMBERS<T> random_numbers;
    ARRAY<ARRAY<int> > triangle_free_particles;
    T refinement_distance;
    bool dynamic_subsampling,temporarily_disable_dynamic_subsampling;
    ARRAY<T> particle_distances;
    int old_number_particles;
    T ring_mass;
    T num_objects_multiplier;
    int number_of_maggots;
    ARRAY<T> periods,phases;
    bool fish_mattress;  // low res mattress for testing
    bool fully_implicit;
    T arg_kd;
    T arg_ks;
    bool use_be,use_tr;
    bool use_free_particles;


    typedef SOLIDS_EXAMPLE<TV> BASE;
    using BASE::solids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type),tests(*this,solid_body_collection),parameter(1),grid_m(10),grid_n(20),
        rigid_body_collection(solid_body_collection.rigid_body_collection),number_of_joints(2),subsamples(8),refinement_distance((T).2),dynamic_subsampling(false),
        temporarily_disable_dynamic_subsampling(false),old_number_particles(0),ring_mass(10000),num_objects_multiplier((T)1),fish_mattress(false),fully_implicit(false),
        use_be(false),use_tr(false),use_free_particles(false)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-parameter",1,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-fully_implicit","use fully implicit forces");
    parse_args->Add_Option_Argument("-fully_explicit","use fully explicit forces");
    parse_args->Add_Option_Argument("-half_fully_implicit","use fully implicit forces for position update");
    parse_args->Add_Option_Argument("-rigid_collisions","use rigid collisions");
    parse_args->Add_Option_Argument("-rigid_contact","use rigid contact");
    parse_args->Add_Double_Argument("-kd",(T)1,"overdamping fraction");
    parse_args->Add_Double_Argument("-ks",(T)100,"spring constant (default: 100)");
    parse_args->Add_Integer_Argument("-grid_m",10,"mesh size");
    parse_args->Add_Integer_Argument("-grid_n",20,"mesh size");
    parse_args->Add_Option_Argument("-use_be","use backward Euler");
    parse_args->Add_Option_Argument("-use_tr","use trapezoid rule");
    parse_args->Add_Option_Argument("-use_free_particles","use free particles to visualize");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    random_numbers.Set_Seed(1234);

    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);

    parameter=parse_args->Get_Integer_Value("-parameter");

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.verbose_dt=true;
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    fully_implicit=parse_args->Is_Value_Set("-fully_implicit");
    solids_parameters.implicit_solve_parameters.use_half_fully_implicit=parse_args->Is_Value_Set("-half_fully_implicit");
    solids_parameters.rigid_body_evolution_parameters.correct_evolution_energy=true;
    solids_parameters.rigid_body_collision_parameters.perform_collisions=parse_args->Is_Value_Set("-rigid_collisions");
    solids_parameters.rigid_body_collision_parameters.perform_contact=parse_args->Is_Value_Set("-rigid_contact");
    arg_ks=parse_args->Get_Double_Value("-ks");
    arg_kd=parse_args->Get_Double_Value("-kd");
    if(parse_args->Is_Value_Set("-fully_explicit")){
        solids_parameters.implicit_solve_parameters.cg_iterations=0;
        solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;}

    grid_m=parse_args->Get_Integer_Value("-grid_m");
    grid_n=parse_args->Get_Integer_Value("-grid_n");

    if(parse_args->Is_Value_Set("-use_tr")) use_tr=true;
    if(parse_args->Is_Value_Set("-use_be")) use_be=true;
    use_free_particles=parse_args->Is_Value_Set("-use_free_particles");

    switch(test_number){
        case 1: last_frame=240;break;
        case 2: last_frame=240;break;
        case 3: last_frame=240;break;
        case 4: last_frame=240;break;
        case 5: last_frame=120;break;
        case 6: last_frame=120;break;
        case 7: last_frame=240;break;
        case 8: last_frame=120;break;
        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}
    output_directory+=STRING_UTILITIES::string_sprintf("_ks_%d_kd_%d",(int)arg_ks,(int)arg_kd);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==2) X(2)=TV((time+(T)1.1),(time+(T)1.1),(time+(T)1.1));
    if(test_number==8){X(1)=particle_curve.Value(time);X(grid_m)=TV(grid_m*2-2,0,0);}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==2) V(2)=TV((T)1,(T)1,(T)1);
    if(test_number==8){V(1)=particle_curve.Derivative(velocity_time);V(grid_m)=TV();}
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==2) V(2)=TV();
    if(test_number==8){V(1)=TV();V(grid_m)=TV();}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1: Hanging_Sphere();break;
        case 2: Dragging_Cube();break;
        case 3: Rigid_Spring();break;
        case 4: Rigid_Spring_Cloth();break;
        case 5: Single_Particle();break;
        case 6: Impulse_Chain();break;
        case 7: Spring_Cloth();break;
        case 8: Particle_Impulse_Chain();break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(use_tr) solids_parameters.use_trapezoidal_rule_for_velocities=true;
    else if(use_be) solids_parameters.use_trapezoidal_rule_for_velocities=false;

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();

    T stiffness=(T)2e6;
    T damping=(T).01;
    for(int i=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        if(test_number==12 || test_number==21) solid_body_collection.Add_Force(Create_Tet_Springs(*tetrahedralized_volume,(T)stiffness/(1+sqrt((T)2)),(T)3));
        else solid_body_collection.Add_Force(Create_Finite_Volume(*tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));}

//     for(int i=1;SEGMENTED_CURVE<TV>* segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>*>(i);i++){
//         solid_body_collection.Add_Force(Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)3));}

    for(int i=1;TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
        solid_body_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping));
        T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
        solid_body_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));
        PHYSBAM_DEBUG_PRINT("Spring stiffnesses",linear_stiffness,linear_damping,bending_stiffness,bending_damping);}

    // disable strain rate CFL for all forces
    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(!dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.deformable_geometry.structures(i))){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && (!dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.deformable_geometry.structures(i))))
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}

    if(use_free_particles){
        FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();
        deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
        free_particles.nodes.Append_Elements(IDENTITY_ARRAY<>(particles.array_collection->Size()));}

    // correct mass
    //binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    if(fully_implicit) for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=0;i<rigid_body_collection.rigids_forces.m;i++) rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==6 && id==1) frame=curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==6 && id==1){twist=curve.Derivative(time);return false;}
    return true;
}
//#####################################################################
// Function Hanging_Sphere
//#####################################################################
void Hanging_Sphere()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;

    RIGID_BODY<TV>& box=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
    box.is_static=true;
    box.X()=TV(0,3,0);
    RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",1,(T)0);
    sphere.X()=TV();

    LOG::cout<<"we currently have particles.array_collection->Size()= "<<particles.array_collection->Size()<<std::endl;
    particles.array_collection->Add_Elements(2);particles.mass(1)=particles.mass(2)=(T)1;
    RIGID_BODY_BINDING<TV>* binding1=new RIGID_BODY_BINDING<TV>(particles,1,rigid_body_collection,box.particle_index,TV(0,-1,0));
    RIGID_BODY_BINDING<TV>* binding2=new RIGID_BODY_BINDING<TV>(particles,2,rigid_body_collection,sphere.particle_index,TV(1,0,0));
    binding_list.Add_Binding(binding1);
    binding_list.Add_Binding(binding2);
    binding1->Clamp_To_Embedded_Position();
    binding1->Clamp_To_Embedded_Velocity();
    binding2->Clamp_To_Embedded_Position();
    binding2->Clamp_To_Embedded_Velocity();
    SEGMENTED_CURVE<TV>* segmented_curve=new SEGMENTED_CURVE<TV>(*new SEGMENT_MESH,particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->mesh.Set_Number_Nodes(2);
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
    tests.Add_Gravity();
}
//#####################################################################
// Function Dragging_Cube
//#####################################################################
void Dragging_Cube()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;

    RIGID_BODY<TV>& box=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
    particles.array_collection->Add_Elements(2);particles.mass(1)=particles.mass(2)=(T)1;
    particles.X(2)=TV((T)1.1,(T)1.1,(T)1.1);
    RIGID_BODY_BINDING<TV>* binding=new RIGID_BODY_BINDING<TV>(particles,1,rigid_body_collection,box.particle_index,TV(0,0,0));
    binding_list.Add_Binding(binding);
    binding->Clamp_To_Embedded_Position();
    binding->Clamp_To_Embedded_Velocity();
    SEGMENTED_CURVE<TV>* segmented_curve=new SEGMENTED_CURVE<TV>(*new SEGMENT_MESH,particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->mesh.Set_Number_Nodes(2);
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
    tests.Add_Gravity();
}
//#####################################################################
// Function Rigid_Spring
//#####################################################################
void Rigid_Spring()
{
    RIGID_BODY<TV>& body_a=tests.Add_Rigid_Body("sphere",1,(T)0);
    RIGID_BODY<TV>& body_b=tests.Add_Rigid_Body("sphere",1,(T)0);
    body_a.Set_Frame(FRAME<TV>(TV(0,0,0)));
    body_b.Set_Frame(FRAME<TV>(TV(3,0,0)));

    TV at[6]={TV(),TV(),TV(1,0,0),TV(0,1,0),TV(0,0,1),TV(1,1,1)};

    RIGID_LINEAR_SPRINGS<TV>* spring=new RIGID_LINEAR_SPRINGS<TV>(rigid_body_collection);
    spring->Add_Spring(body_a.particle_index,body_b.particle_index,at[parameter/10],at[parameter%10]);
    spring->Set_Restlength(0,5);
    spring->Set_Stiffness(0,arg_ks);
    spring->Set_Overdamping_Fraction(0,arg_kd);
    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
}
//#####################################################################
// Function Rigid_Spring_Cloth
//#####################################################################
void Rigid_Spring_Cloth()
{
    ARRAY<int,VECTOR<int,2> > body_indices(1,grid_m,1,grid_n);
    for(int i=0;i<grid_m;i++)for(int j=0;j<grid_n;j++){
        RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",1,(T)0);
        body.Set_Frame(FRAME<TV>(TV(2*i,0,2*j)));
        body_indices(i,j)=body.particle_index;}

    RIGID_LINEAR_SPRINGS<TV>* spring=new RIGID_LINEAR_SPRINGS<TV>(rigid_body_collection);
    int spring_index=0;
    for(int i=0;i<grid_m;i++)for(int j=0;j<grid_n;j++){
        if(i<grid_m){
            spring->Add_Spring(body_indices(i,j),body_indices(i+1,j),TV(),TV());
            spring_index++;
            spring->Set_Restlength(spring_index,5);
            spring->Set_Stiffness(spring_index,arg_ks);
            spring->Set_Overdamping_Fraction(spring_index,arg_kd);}
        if(j<grid_n){
            spring->Add_Spring(body_indices(i,j),body_indices(i,j+1),TV(),TV());
            spring_index++;
            spring->Set_Restlength(spring_index,5);
            spring->Set_Stiffness(spring_index,arg_ks);
            spring->Set_Overdamping_Fraction(spring_index,arg_kd);}}
    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
//    solids_parameters.use_post_cg_constraints=false;
}
//#####################################################################
// Function Single_Particle
//#####################################################################
void Single_Particle()
{
    RIGID_BODY<TV>& body_a=tests.Add_Rigid_Body("sphere",1,(T)0);body_a.Set_Mass(FLT_MAX);
    RIGID_BODY<TV>& body_b=tests.Add_Rigid_Body("sphere",1,(T)0);
    RIGID_BODY<TV>& body_c=tests.Add_Rigid_Body("sphere",1,(T)0);body_c.Set_Mass(FLT_MAX);
    body_a.Set_Frame(FRAME<TV>(TV(-3,0,0)));
    body_b.Set_Frame(FRAME<TV>(TV(.5,0,0)));
    body_c.Set_Frame(FRAME<TV>(TV(3,0,0)));

    RIGID_LINEAR_SPRINGS<TV>* spring=new RIGID_LINEAR_SPRINGS<TV>(rigid_body_collection);
    spring->Add_Spring(body_a.particle_index,body_b.particle_index,TV(),TV());
    spring->Set_Restlength(0,5);
    spring->Set_Stiffness(0,arg_ks);
    spring->Set_Overdamping_Fraction(0,arg_kd);
    spring->Add_Spring(body_c.particle_index,body_b.particle_index,TV(),TV());
    spring->Set_Restlength(1,5);
    spring->Set_Stiffness(1,arg_ks);
    spring->Set_Overdamping_Fraction(1,arg_kd);

    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
}
//#####################################################################
// Function Impulse_Chain
//#####################################################################
void Impulse_Chain()
{
    ARRAY<int> body_indices;
    for(int i=0;i<grid_m;i++){
        RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T).5,(T)0);
        body.Set_Frame(FRAME<TV>(TV(2*i-2,0,0)));
        body_indices.Append(body.particle_index);}

    solid_body_collection.rigid_body_collection.Rigid_Body(body_indices(0)).Is_Kinematic()=true;
    solid_body_collection.rigid_body_collection.Rigid_Body(body_indices.Last()).is_static=true;

    RIGID_LINEAR_SPRINGS<TV>* spring=new RIGID_LINEAR_SPRINGS<TV>(rigid_body_collection);
    for(int i=0;i<grid_m-1;i++){
        spring->Add_Spring(body_indices(i),body_indices(i+1),TV(),TV());
        spring->Set_Restlength(i,2);
        spring->Set_Stiffness(i,arg_ks);
        spring->Set_Overdamping_Fraction(i,arg_kd);}
    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;

    curve.Add_Control_Point(0,FRAME<TV>());
    if(parameter==2) curve.Add_Control_Point((T)10,FRAME<TV>(TV(-20,0,0)));
    else{
        curve.Add_Control_Point((T).25,FRAME<TV>(TV(-2,0,0)));
        curve.Add_Control_Point((T).5,FRAME<TV>());}
}
//#####################################################################
// Function Spring_Cloth
//#####################################################################
void Spring_Cloth()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    particles.array_collection->Add_Elements(grid_m*grid_n);

    SEGMENTED_CURVE<TV> *segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);    

    for(int i=0;i<grid_m;i++) for(int j=0;j<grid_n;j++){
        particles.X((i-1)*grid_n+j)=TV(2*i-2,0,2*j-2);
        if(i<grid_m) segmented_curve->mesh.elements.Append(VECTOR<int,2>((i-1)*grid_n+j,i*grid_n+j));
        if(j<grid_n) segmented_curve->mesh.elements.Append(VECTOR<int,2>((i-1)*grid_n+j,(i-1)*grid_n+j+1));}

    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    particles.mass.Fill((T)1);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    ARRAY<T> restlength(grid_m*(grid_n-1)+(grid_m-1)*grid_n);
    restlength.Fill((T)5);
    LINEAR_SPRINGS<TV>* spring=new LINEAR_SPRINGS<TV>(particles,segmented_curve->mesh,fully_implicit);
    spring->Set_Restlength(restlength);
    spring->Set_Stiffness(arg_ks);
    spring->Set_Overdamping_Fraction(arg_kd);

    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;
//    solids_parameters.use_post_cg_constraints=false;
}
//#####################################################################
// Function Particle_Impulse_Chain
//#####################################################################
void Particle_Impulse_Chain()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    particles.array_collection->Add_Elements(grid_m);

    SEGMENTED_CURVE<TV> *segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);    

    for(int i=0;i<grid_m;i++){
        particles.X(i)=TV(2*i-2,0,0);
        if(i<grid_m-1) segmented_curve->mesh.elements.Append(VECTOR<int,2>(i,i+1));}

    particles.mass.Fill((T)1);
    particles.mass(0)=FLT_MAX;
    particles.mass(grid_m-1)=FLT_MAX;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();

    ARRAY<T> restlength(grid_m-1);
    restlength.Fill((T)2);
    LINEAR_SPRINGS<TV>* spring=new LINEAR_SPRINGS<TV>(particles,segmented_curve->mesh,fully_implicit);
    spring->Set_Restlength(restlength);
    spring->Set_Stiffness(arg_ks);
    spring->Set_Overdamping_Fraction(arg_kd);
    solid_body_collection.Add_Force(spring);
    solids_parameters.use_trapezoidal_rule_for_velocities=false;

    particle_curve.Add_Control_Point(0,TV());
    if(parameter==2) particle_curve.Add_Control_Point((T)10,TV(-20,0,0));
    else{
        particle_curve.Add_Control_Point((T).25,TV(-2,0,0));
        particle_curve.Add_Control_Point((T).5,TV());}
}
//#####################################################################
};
}
#endif
