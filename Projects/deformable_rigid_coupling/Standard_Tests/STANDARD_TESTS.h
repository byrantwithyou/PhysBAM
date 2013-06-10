//#####################################################################
// Copyright 2007-2008, Craig Schroeder, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. example using only push out
//   2. stack with rigid and deformable blocks
//   3. trampoline
//   4. two cubes and two particles slide along a plane
//   5. Point joint between 2 blocks connected to a deformable body
//   6. Point joint between 2 blocks connected to a deformable body with one rigid body suspended
//   7. Deformable torsion spring connected to kinematic rigid body and articulated rigid body
//   8. Floppy fish
//   9. Rigid/deformable/articulated hanging chain
//  10. Rigid/deformable/articulated twisting chain
//  11. Brush and wheel
//  12. Ring Drop
//  13. rigid and deformable cubes slide down incline plane with friction
//  14. Sidewinding
//  15. Row of spheres
//  16. Magget
//  17. Maggets in a bowl
//  18. Maggets and blocks in a bowl
//  19. Trapped Magget
//  20. Push out test
//  21. Sphere fall
//  22. Deformable sphere on rigid box
//  23. Rigid box on deformable sphere
//  24. Two-way coupling between rigid and embedded tori
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_MESH.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_FORCE_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    int parameter;
    int kinematic_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;

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
    bool use_forces_for_drift;
    bool project_nullspace;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;
    using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::test_number;using BASE::parse_args;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),parameter(1),
        rigid_body_collection(solid_body_collection.rigid_body_collection),number_of_joints(2),subsamples(8),refinement_distance((T).2),dynamic_subsampling(false),
        temporarily_disable_dynamic_subsampling(false),old_number_particles(0),ring_mass(10000),num_objects_multiplier((T)1),fish_mattress(false),fully_implicit(false),
        use_forces_for_drift(false),project_nullspace(false)
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
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_Deformable_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Set_Rigid_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

ROTATION<TV> Upright_Orientation(const TV& x,const TV& y)
{
    TV u=(y-x).Normalized(),j(0,1,0);
    MATRIX<T,3> rotation(u,j.Projected_Orthogonal_To_Unit_Direction(u),TV::Cross_Product(u,j));
    return ROTATION<TV>(rotation);
}
TV Sidewinding_Position(const T time,const T segment)
{
    const T k=(T)two_pi/wavelength,w=(T)two_pi/period,theta=segment*k+time*w;
    T y=height_amplitude*cos(theta);
    T x=bend_amplitude*sin(theta);
    return TV(x,y,segment);
}
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add("-parameter",&parameter,"value","parameter used by multiple tests to change the parameters of the test");
    parse_args->Add("-subsamples",&subsamples,"samples","number of particle subsamples per triangle");
    parse_args->Add("-n",&num_objects_multiplier,"scale","multiplier for number of objects in drop test");
    parse_args->Add("-mass",&ring_mass,"mass","mass of objects in drop test");
    parse_args->Add("-sp",&solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,"Use spatial partition for collision body collisions");
    parse_args->Add("-print_energy",&solid_body_collection.solid_force_collection.print_energy,"print energy statistics");
    parse_args->Add("-fully_implicit",&fully_implicit,"use fully implicit forces");
    parse_args->Add("-project_nullspace",&project_nullspace,"project out nullspace");
    parse_args->Add("-binding_springs",&use_forces_for_drift,"use binding springs for drift particles");
    parse_args->Add("-test_system",&solids_parameters.implicit_solve_parameters.test_system,"Test Krylov system");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    if(project_nullspace) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Self_Collisions_Begin_Callback
//#####################################################################
void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(dynamic_subsampling && solids_parameters.triangle_collision_parameters.perform_self_collision) Update_Subsamples();
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE
{
    if(dynamic_subsampling && !solids_parameters.triangle_collision_parameters.perform_self_collision) Update_Subsamples();
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY_COLLISIONS<TV>& collisions=*solids_evolution->rigid_body_collisions;
    if(test_number==1){
        solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
        solids_parameters.rigid_body_collision_parameters.collision_iterations=0;
        solids_parameters.rigid_body_collision_parameters.contact_iterations=0;
        collisions.Set_Shock_Propagation_Iterations(0);
        collisions.Set_Push_Out_Iterations(1);
        collisions.Set_Push_Out_Level_Iterations(1);
        collisions.Set_Push_Out_Pair_Iterations(1);
        collisions.Use_Freezing_With_Push_Out(false);}
    if(test_number==2){collisions.Set_Shock_Propagation_Iterations(0);collisions.Use_Freezing_With_Push_Out(false);}
    if(test_number==8){
        T desired_x=(T)two_pi/16;
        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x*sin(4*time),TV(0,1,0));
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(desired_rotation);}}
    if(test_number==14){
        T sim_time=max((T)0,time-start_time);
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            TV previous=Sidewinding_Position(sim_time,(i-1)*joint_separation);
            TV current=Sidewinding_Position(sim_time,i*joint_separation);
            TV next=Sidewinding_Position(sim_time,(i+1)*joint_separation);
            ROTATION<TV> previous_orientation=Upright_Orientation(previous,current); // parent
            ROTATION<TV> next_orientation=Upright_Orientation(current,next); // child
            ROTATION<TV> joint_frame=joint.frame_jp.r*initial_orientation.Inverse()*previous_orientation.Inverse()*next_orientation*initial_orientation*joint.frame_cj.r;
            if(time<start_time) joint_frame=joint_frame.Scale_Angle(time/start_time);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(joint_frame);}}
    if(test_number==16 || test_number==19){
        T sim_time=max((T)0,time-start_time);
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            ROTATION<TV> previous_orientation;
            ROTATION<TV> next_orientation(sin((T)two_pi*sim_time/period),TV(0,1,0));
            ROTATION<TV> joint_frame=joint.frame_jp.r*initial_orientation.Inverse()*previous_orientation.Inverse()*next_orientation*initial_orientation*joint.frame_cj.r;
            if(time<start_time) joint_frame=joint_frame.Scale_Angle(time/start_time);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(joint_frame);}}
    if(test_number==17 || test_number==18){
        T sim_time=max((T)0,time-start_time);
        for(int i=0;i<arb.joint_mesh.joints.m;i++){JOINT<TV>& joint=*arb.joint_mesh.joints(i);
            int maggot=(i+1)/2;
            ROTATION<TV> previous_orientation;
            ROTATION<TV> next_orientation(sin((T)two_pi*sim_time/periods(maggot)+phases(maggot)),TV(0,1,0));
            ROTATION<TV> joint_frame=joint.frame_jp.r*initial_orientation.Inverse()*previous_orientation.Inverse()*next_orientation*initial_orientation*joint.frame_cj.r;
            if(time<start_time) joint_frame=joint_frame.Scale_Angle(time/start_time);
            if(joint.joint_function) joint.joint_function->Set_Target_Angle(joint_frame);}}
    if(test_number==12){
        const bool do_collisions=time>=(T)1.45;
        solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=do_collisions;
        solids_parameters.triangle_collision_parameters.perform_self_collision=do_collisions;
        temporarily_disable_dynamic_subsampling=!do_collisions;
        solids_parameters.write_static_variables_every_frame=true;}
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    random_numbers.Set_Seed(1234);
    fluids_parameters.simulate=false;
    solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
    solids_parameters.cfl=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solids_parameters.verbose_dt=true;

    switch(test_number){
        case 1: Push_Out_Test();break;
        case 2: Coupled_Stack();break;
        case 3: Trampoline();break;
        case 4: Sliding_Test();break;
        case 5: Point_Joint_Coupled_Drop(false);break;
        case 6: Point_Joint_Coupled_Drop(true);break;
        case 7: Torsion_Spring();break;
        case 8: Floppy_Fish();break;
        case 9: Twisting_Chain(false);break;
        case 10: Twisting_Chain(true);break;
        case 11: Brush_And_Wheel();break;
        case 12: Ring_Drop();break;
        case 13: Cubes_Friction();break;
        case 14: Sidewinding();break;
        case 15: Row_Of_Spheres();break;
        case 16: Maggot();break;
        case 17: Bowl_Of_Maggots(false);break;
        case 18: Bowl_Of_Maggots(true);break;
        case 19: Trapped_Maggot();break;
        case 20: Push_Out_Test2();break;
        case 21: Sphere_Fall();break;
        case 22: Sphere_Block(true);break;
        case 23: Sphere_Block(false);break;
        case 24: Two_Way_Tori();break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    if(dynamic_subsampling) Initialize_Dynamic_Subsampling();

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    T stiffness=(T)2e6;
    if(test_number==3) stiffness=(T)5e3;
    else if(test_number==8) stiffness=(T)4e5;
    else if(test_number==12) stiffness=(T)1e5;
    else if(test_number==13) stiffness=(T)2e6;
    else if(test_number==14) stiffness=(T)2e5;
    else if(test_number==16) stiffness=(T)2e5;
    else if(test_number==17 || test_number==18) stiffness=(T)2e5;
    else if(test_number==19) stiffness=(T)2e5;
    else if(test_number==21) stiffness=(T)1e5;
    else if(test_number==39 || test_number==40) stiffness=(T)2e5;
    T damping=(T).01;
    if(test_number==8) damping=(T).03;
    if(test_number==19) damping=(T).1;
    for(int i=0;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        if(test_number==12 || test_number==21) solid_body_collection.solid_force_collection.Add_Force(Create_Tet_Springs(*tetrahedralized_volume,(T)stiffness/(1+sqrt((T)2)),(T)3));
        else
            if(test_number==8 || test_number==14){
                FINITE_VOLUME<TV,3>* fv=Create_Finite_Volume(*tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1,
                    true,false,true,(PLASTICITY_MODEL<T,3>*)0,true);
                fv->density=(T)1000;
                solid_body_collection.solid_force_collection.Add_Force(fv);}
            else solid_body_collection.solid_force_collection.Add_Force(Create_Finite_Volume(*tetrahedralized_volume,new NEO_HOOKEAN<T,3>(stiffness,(T).45,damping,(T).25),true,(T).1));}

    for(int i=0;SEGMENTED_CURVE<TV>* segmented_curve=deformable_body_collection.template Find_Structure<SEGMENTED_CURVE<TV>*>(i);i++){
        solid_body_collection.solid_force_collection.Add_Force(Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)3));}
    
    for(int i=0;TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
        solid_body_collection.solid_force_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping));
        T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
        solid_body_collection.solid_force_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));
        PHYSBAM_DEBUG_PRINT("Spring stiffnesses",linear_stiffness,linear_damping,bending_stiffness,bending_damping);}

    for(int i=0;i<deformable_body_collection.structures.m;i++) if(!dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.structures(i))){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision && (!dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.structures(i))))
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    if(test_number==8){
        // collide structures with the ground only
        deformable_body_collection.collisions.Use_Structure_Collide_Collision_Body(true);
        for(int s=0;s<deformable_body_collection.structures.m;s++)
            deformable_body_collection.collisions.structure_collide_collision_body(s).Set(rigid_body_collection.collision_body_list->geometry_id_to_collision_geometry_id.Get(ground->particle_index));}

    // correct mass
    //binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // add forces
    if(test_number!=1) tests.Add_Gravity();

    // disable strain rate CFL for all forces
    for(int i=0;i<rigid_body_collection.rigid_force_collection.rigids_forces.m;i++) rigid_body_collection.rigid_force_collection.rigids_forces(i)->limit_time_step_by_strain_rate=false;
    for(int i=0;i<deformable_body_collection.deformable_force_collection.deformables_forces.m;i++) deformable_body_collection.deformable_force_collection.deformables_forces(i)->limit_time_step_by_strain_rate=false;
    for(int i=0;i<solid_body_collection.solid_force_collection.solids_forces.m;i++) solid_body_collection.solid_force_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;

    for(int i=0;i<rigid_body_collection.rigid_force_collection.rigids_forces.m;i++) rigid_body_collection.rigid_force_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
    for(int i=0;i<deformable_body_collection.deformable_force_collection.deformables_forces.m;i++) deformable_body_collection.deformable_force_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
    for(int i=0;i<solid_body_collection.solid_force_collection.solids_forces.m;i++) solid_body_collection.solid_force_collection.solids_forces(i)->use_implicit_velocity_independent_forces=fully_implicit;
}
//#####################################################################
// Function Advance_One_Time_Step_End_Callback
//#####################################################################
void Advance_One_Time_Step_End_Callback(const T dt,const T time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(test_number==12){
        for(int i=0;i<rigid_body_clamp_time.m;i++) if(rigid_body_clamp_time(i).y>=time){RIGID_BODY<TV>& rigid_body=*rigid_body_clamp_time(i).x;
            T speed=rigid_body.Twist().linear.Magnitude();
            if(speed>maximum_fall_speed) rigid_body.Twist().linear*=maximum_fall_speed/speed;}

        for(int i=0;i<structure_clamp_time.m;i++) if(structure_clamp_time(i).y>=time){TETRAHEDRALIZED_VOLUME<T>& volume=*structure_clamp_time(i).x;
            for(int j=0;j<volume.mesh.elements.m;j++){const VECTOR<int,4>& e=volume.mesh.elements(j);
                for(int k=0;k<e.Size();k++){
                    T speed=deformable_body_collection.particles.V(e(k)).Magnitude();
                    if(speed>maximum_fall_speed) deformable_body_collection.particles.V(e(k))*=maximum_fall_speed/speed;}}}}
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==kinematic_id){
        frame.r=torus_rotation;
        if(time<start_time || time>start_time+period) frame.t=torus_min;
        else{TV half=(T).5*(torus_max-torus_min);frame.t=torus_min+half-half*cos((T)two_pi/period*(time-start_time));}}
    if((test_number==4 || test_number==13) && id==kinematic_id){
        T magnitude=max((T)0,(T).5*sqr(time)*(T)9.8*(sin(ground_angle_rad)-mu*cos(ground_angle_rad)));
        TV direction(-cos(ground_angle_rad),-sin(ground_angle_rad),0);
        frame.t=magnitude*direction-TV(0,0,(T).25);}
    if((test_number==7 || test_number==10 || test_number==11) && id==kinematic_id) frame=curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==kinematic_id){
        twist.angular=TV();
        if(time<start_time || time>start_time+period) twist.linear=TV();
        else{T w=(T)two_pi/period;TV half=(T).5*(torus_max-torus_min);twist.linear=w*half*sin(w*(time-start_time));}}
    if((test_number==4 || test_number==13) && id==kinematic_id){
        T magnitude=max((T)0,time*(T)9.8*(sin(ground_angle_rad)-mu*cos(ground_angle_rad)));
        TV direction(-cos(ground_angle_rad),-sin(ground_angle_rad),0);
        twist.linear=magnitude*direction;}
    if((test_number==7 || test_number==10 || test_number==11) && id==kinematic_id) twist=curve.Derivative(time);
    return true;
}
//#####################################################################
// Function Push_Out_Test
//#####################################################################
void Push_Out_Test()
{
    last_frame=120;
    // set up two blocks to be interpenetrating
    int num_blocks=4;
    for(int i=0;i<num_blocks;i++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
        rigid_body.Frame()=FRAME<TV>(TV((T)1.8*(i+1),(T)1.6*(i+1),(T)0),ROTATION<TV>((T)pi/10*(i+1),TV::Axis_Vector(0)));
        rigid_body.Set_Mass(1000);
        if(i==2) rigid_body.is_static=true;
        /*        rigid_body.Twist().linear=TV(-1,0,0);*/}

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)1.5,(T)0))),true,true,1000);
    tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume,solids_parameters.triangle_collision_parameters);

    TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
        RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)-1.7,(T)2,(T)0))),true,true,1000);
    tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume2,solids_parameters.triangle_collision_parameters);

//    int n=10;
//    for(int i=0;i<n;i++) for(int j=0;j<n;j++) for(int k=0;k<n;k++) tests.Add_Rigid_Body("subdivided_box",1,(T)0).Frame().t=TV(1.8*(i+1),1.*j,1.75*k);
}
//#####################################################################
// Function Coupled_Stack
//#####################################################################
void Coupled_Stack()
{
    last_frame=240;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    dynamic_subsampling=true;
    for(int i=0;i<6;i++){
        if(i%2==0){
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",1,(T).4);
            rigid_body.Frame().t=TV(0,(T)(2*i+1),0);
            rigid_body.Set_Mass(6000);}
        else{
            tests.Create_Mattress(GRID<TV>(VECTOR<int,3>(7,7,7),RANGE<TV>(TV(-(T)1,(T)(2*i),-(T)1),TV((T)1,(T)(2*i+2),(T)1))),1000);}}
 
    maximum_fall_speed=105;
//    T v=maximum_fall_speed,g=solids_parameters.gravity,t=4;
    RIGID_BODY<TV>& sphere=tests.Add_Rigid_Body("sphere",(T)1,(T).2);
    sphere.Frame().t=TV(0,1,82);
    sphere.Twist().linear=TV(0,0,-25);
    sphere.Update_Bounding_Box();
    sphere.Set_Mass(6000);

    tests.Add_Ground((T).4);
}
//#####################################################################
// Function Find_Placement
//#####################################################################
FRAME<TV> Find_Placement(RANDOM_NUMBERS<T>& random,const RANGE<TV>& bounding_box,ARRAY<ORIENTED_BOX<TV> >& bounding_boxes,const RANGE<TV>& world,bool want_rotate)
{
    for(int i=0;i<10000;i++){
        FRAME<TV> frame;
        if(want_rotate) frame.r=random.template Get_Rotation<TV>();
        ORIENTED_BOX<TV> oriented_box(bounding_box,frame.r);
        RANGE<TV> new_box(oriented_box.Bounding_Box());
        frame.t=random.Get_Uniform_Vector(world.min_corner-new_box.min_corner,world.max_corner-new_box.max_corner);
        oriented_box.corner+=frame.t;
        bool okay=true;
        for(int j=0;j<bounding_boxes.m;j++) if(oriented_box.Intersection(bounding_boxes(j))){okay=false;break;}
        if(okay){
            bounding_boxes.Append(oriented_box);
            return frame;}}
    PHYSBAM_FATAL_ERROR("Could not find suitable placement");
}
//#####################################################################
// TRAMPOLINE_CLIPPING
//#####################################################################
class TRAMPOLINE_CLIPPING:public TRIANGULATED_SURFACE_CLIPPING_HELPER<T>
{
public:
    TV axis;
    T radius_squared;
    TRAMPOLINE_CLIPPING(const TV& axis_input,const T radius)
        :axis(axis_input.Normalized()),radius_squared(sqr(radius))
    {}

    bool Prune(const TV& pt) const
    {return pt.Projected_Orthogonal_To_Unit_Direction(axis).Magnitude_Squared()>radius_squared;}

    void operator()(TRIANGULATED_SURFACE<T>& surface) const
    {
        for(int i=surface.mesh.elements.m-1;i>=0;i--){const VECTOR<int,3>& elements(surface.mesh.elements(i));
            for(int j=0;j<3;j++) if(Prune(surface.particles.X(elements(j)))){surface.mesh.elements.Remove_Index_Lazy(i);break;}}
        surface.Discard_Valence_Zero_Particles_And_Renumber();
        static_cast<DEFORMABLE_PARTICLES<TV>&>(surface.particles).mass*=20;
    }
};
//#####################################################################
// Function Trampoline
//#####################################################################
void Trampoline()
{
    last_frame=120;
    RANDOM_NUMBERS<T> random;
    random.Set_Seed(1234);
    solids_parameters.cfl=(T)1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)7;
//    solids_parameters.enforce_repulsions_in_cg=true;
//    solids_parameters.implicit_solve_parameters.cg_projection_iterations=25;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    aspect_ratio=1;
    number_side_panels=75;
    stiffness_multiplier=1000;
    damping_multiplier=(T)1;
    bending_stiffness_multiplier=10;
    bending_damping_multiplier=(T)1;
//    int n=2,levels=4;
    T separation=(T)2.4;
    torus_rotation=ROTATION<TV>((T)half_pi,TV(1,0,0));
    torus_min=TV(0,3,0);
    torus_max=TV(0,5,0);
    period=2;
    start_time=(T)2;

    RIGID_BODY<TV>& torus=tests.Add_Analytic_Torus((T).3,(T)7,12,48);
    torus.Is_Kinematic()=true;
    kinematic_id=torus.particle_index;
    torus.Frame().t=torus_min;
    torus.Frame().r=torus_rotation;

    TORUS<T> torus_geometry(static_cast<ANALYTIC_IMPLICIT_OBJECT<TORUS<T> >&>(*torus.implicit_object->object_space_implicit_object).analytic);
    TRAMPOLINE_CLIPPING clipping(torus.World_Space_Vector(torus_geometry.axis),torus_geometry.inner_radius+torus_geometry.outer_radius);
    TRIANGULATED_SURFACE<T>& cloth=tests.Create_Cloth_Panel(number_side_panels,(T)15,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,3,0))),&clipping);
    (void)cloth;

    tests.Bind_Particles_In_Rigid_Body(torus);

    ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
    RANGE<TV> world(TV(-2*separation,4,-2*separation),TV(2*separation,15,2*separation));

     for(int i=0;i<3;i++){
         RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",(T).6,(T).2);
         rigid_body.Update_Bounding_Box();
         rigid_body.Frame()=Find_Placement(random,rigid_body.axis_aligned_bounding_box,bounding_boxes,world,false);
         rigid_body.Set_Mass(1);}
     for(int i=0;i<3;i++){
         RANGE<TV> box(-TV((T)3,(T)3,(T)1),TV((T)3,(T)3,(T)1));
         RIGID_BODY_STATE<TV> state(Find_Placement(random,box*(T).6,bounding_boxes,world,true));
         TETRAHEDRALIZED_VOLUME<T>& object=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_1K.tet",state,false,false,1000,(T).6);
         SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(object,(T)10);}
     tests.Add_Ground((T)1,-(T)2);
}
//#####################################################################
// Function Sliding_Test
//#####################################################################
void Sliding_Test()
{
    last_frame=96;
    RIGID_BODY<TV>* rigid_body=0;

    const char* boxfile="box";
    T box_size=(T)0.05;
    T ground_angle_deg=20;
    ground_angle_rad=ground_angle_deg*((T)pi/180);
    mu=(T).2;

    T initial_speed_down_ramp1=(T)0;

    T initial_speed_down_ramp2=(T)0;

    LOG::cout<<"Using angle "<<ground_angle_deg<<" deg ("<<ground_angle_rad<<") (critical coeff="<<tan(ground_angle_rad)<<std::endl;

    rigid_body=&tests.Add_Rigid_Body(boxfile,box_size,mu);
    rigid_body->Frame().t=TV(-box_size*sin(ground_angle_rad),box_size*cos(ground_angle_rad)+(T).2,0);
    rigid_body->Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));
    rigid_body->Twist().linear=rigid_body->Frame().r.Rotate(TV(-initial_speed_down_ramp1,0,0));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name=STRING_UTILITIES::string_sprintf("box 1 mu=%.2f",mu);

    rigid_body=&tests.Add_Rigid_Body(boxfile,box_size,mu);
    rigid_body->Frame().t=TV(-box_size*sin(ground_angle_rad),box_size*cos(ground_angle_rad),8*box_size);
    rigid_body->Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Twist().linear=rigid_body->Frame().r.Rotate(TV(-initial_speed_down_ramp2,0,0));
    rigid_body->name=STRING_UTILITIES::string_sprintf("box 2 mu=%.2f",mu);

    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();solid_body_collection.deformable_body_collection.Add_Structure(&free_particles);
    int particle_rest=solid_body_collection.deformable_body_collection.particles.Add_Element(),particle_fall=solid_body_collection.deformable_body_collection.particles.Add_Element();
    solid_body_collection.deformable_body_collection.particles.mass(particle_rest)=1;solid_body_collection.deformable_body_collection.particles.mass(particle_fall)=1;
    solid_body_collection.deformable_body_collection.particles.X(particle_rest)=TV(0,0,(T).8);solid_body_collection.deformable_body_collection.particles.X(particle_fall)=TV(0,(T).2,(T)1.4);
    free_particles.nodes.Append(particle_rest);free_particles.nodes.Append(particle_fall);

    rigid_body=&tests.Add_Rigid_Body("sphere",(T).1,0);
    rigid_body->Frame().t=TV(-box_size*sin(ground_angle_rad),box_size*cos(ground_angle_rad)+(T).2,0);
    rigid_body->name="analytic solution";
    rigid_body->Is_Kinematic()=true;
    kinematic_id=rigid_body->particle_index;

    RIGID_BODY<TV>& ground=tests.Add_Ground(mu,0,0);
    ground.Set_Coefficient_Of_Rolling_Friction(1);
    ground.Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));
}
//#####################################################################
// Function Point_Joint_Coupled_Drop
//#####################################################################
void Point_Joint_Coupled_Drop(bool make_body_static)
{
    last_frame=120;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body1.Frame().t=TV(0,2,0);
    rigid_body1.Set_Mass(1000);
    rigid_body1.name="parent";
    if(make_body_static) rigid_body1.is_static=true;

    RIGID_BODY<TV>& rigid_body2=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body2.Frame().t=TV(0,4,2);
    rigid_body2.Set_Mass(1000);
    rigid_body2.name="child";

    POINT_JOINT<TV>* joint=new POINT_JOINT<TV>();
    arb.joint_mesh.Add_Articulation(rigid_body1.particle_index,rigid_body2.particle_index,joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(1,1,1)));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(1,-1,-1)));

    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)5.5,(T)3.5))),true,true,1000);
    tests.Bind_Particles_In_Rigid_Body(rigid_body2);
    tests.Add_Ground((T).5,-2,1);
}
//#####################################################################
// Function Torsion_Spring
//#####################################################################
void Torsion_Spring()
{
    last_frame=120;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.implicit_solve_parameters.cg_projection_iterations=3;
    RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body1.Frame().t=TV(0,2,0);
    rigid_body1.Set_Mass(30);
    rigid_body1.name="parent";

    RIGID_BODY<TV>& rigid_body2=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body2.Frame().t=TV((T)2.1,2,0);
    rigid_body2.name="child";
    rigid_body2.Set_Mass(30);
    rigid_body2.is_static=true;

    JOINT<TV>* joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(rigid_body1.particle_index,rigid_body2.particle_index,joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV((T)1.05,0,0),ROTATION<TV>()));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(-(T)1.05,0,0),ROTATION<TV>()));

    RIGID_BODY<TV>& rigid_body3=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body3.Frame().t=TV(-8,2,0);
    rigid_body3.Set_Mass(30);
    rigid_body3.name="kinematic";
    rigid_body3.Is_Kinematic()=true;
    kinematic_id=rigid_body3.particle_index;
    curve.Add_Control_Point((T)0,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>()));
    curve.Add_Control_Point((T)1,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>()));
    curve.Add_Control_Point((T)2,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>((T)(.5*pi),TV(1,0,0))));
    curve.Add_Control_Point((T)3,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>((T)(pi),TV(1,0,0))));
    curve.Add_Control_Point((T)4,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>((T)(1.5*pi),TV(1,0,0))));
    curve.Add_Control_Point((T)5,FRAME<TV>(rigid_body3.Frame().t,ROTATION<TV>((T)(2*pi),TV(1,0,0))));

    RIGID_BODY<TV>& rigid_body4=tests.Add_Rigid_Body("subdivided_box",1,(T).5);
    rigid_body4.Frame().t=TV((T)1.05,(T)4.05,0);
    rigid_body4.name="drop";
    rigid_body4.Set_Mass(30);

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    tests.Create_Mattress(GRID<TV>(VECTOR<int,3>(33,5,5),RANGE<TV>(TV(-8,(T)1.5,-.5),TV(0,(T)2.5,.5))),1000);

    tests.Bind_Particles_In_Rigid_Body(rigid_body1);
    tests.Bind_Particles_In_Rigid_Body(rigid_body3);
    tests.Add_Ground((T).5,-2,1);
}
//#####################################################################
// Function Floppy_Fish
//#####################################################################
void Initialize_Joint_Between(JOINT<TV>* joint,const RIGID_BODY<TV>& parent,const RIGID_BODY<TV>& child,TV up)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    arb.joint_mesh.Add_Articulation(parent.particle_index,child.particle_index,joint);
    const FRAME<TV>& pf=parent.Frame();
    const FRAME<TV>& cf=child.Frame();
    TV x=(cf.t-pf.t).Normalized();
    TV y=up.Projected_Orthogonal_To_Unit_Direction(x).Normalized();
    TV z=TV::Cross_Product(x,y);
    FRAME<TV> J((T).5*(cf.t+pf.t),ROTATION<TV>(MATRIX<T,3>(x,y,z)));
    joint->Set_Parent_To_Joint_Frame(J.Inverse_Times(pf));
    joint->Set_Child_To_Joint_Frame(J.Inverse_Times(cf));
}
//#####################################################################
// Function Floppy_Fish
//#####################################################################
void Floppy_Fish()
{
    last_frame=240;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.enforce_poststabilization_in_cg=false;
    solids_parameters.cfl=(T)2;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    T friction=(T).2;
    TV center(0,3,0);
//    T scale=(T).55;
//    T spacing=(T).25;
//    T plank_length=(T)1;
//    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);

    ARRAY<RIGID_BODY<TV>*> bones;

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Frame()=FRAME<TV>(TV(-(T)2,(T)3,0));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Frame()=FRAME<TV>(TV(-(T).75,(T)3,0));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T)1,friction));
    bones.Last()->Frame()=FRAME<TV>(TV((T).5,(T)3,0));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T).7,friction));
    bones.Last()->Frame()=FRAME<TV>(TV((T)1.7,(T)3,0));

    bones.Append(&tests.Add_Rigid_Body("miniplank25wide2",(T).5,friction));
    bones.Last()->Frame()=FRAME<TV>(TV((T)2.5,(T)3,0));

    T density=1000;
    for(int i=0;i<bones.m;i++){
        bones(i)->Set_Mass(density*bones(i)->Volume());}

    T joint_strengths[4]={(T)1000,(T)1000,(T)400,(T)200};
    for(int i=2;i<=bones.m;i++){
        JOINT<TV>* joint=new POINT_JOINT<TV>;
        Initialize_Joint_Between(joint,*bones(i-1),*bones(i),TV(0,0,1));
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(joint_strengths[i-2]);joint_function->Set_Target_Angle(ROTATION<TV>());}

//    // make the joints
//    tests.PD_Curl(scale,center-(T).5*(extent-scale*plank_length)*TV::Axis_Vector(0),ROTATION<TV>(),25,number_of_joints,false);

    // add the fish
    RIGID_BODY_STATE<TV> fish_state(FRAME<TV>(TV(0,3,0),ROTATION<TV>((T)half_pi,TV(1,0,0))));
    if(fish_mattress){
        RANGE<TV> box(-TV((T)3.93278,(T)1.07277,(T)0.384066),TV((T)2.68344,(T)1.1747,(T)0.384353));VECTOR<int,3> counts(20,15,5);
        GRID<TV> mattress_grid=GRID<TV>(counts,box);
        tests.Create_Mattress(mattress_grid,true,&fish_state,1000);}
    else tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",fish_state,true,false,1000,(T)1);

    // binding the deformable particles to the rigid bodies
    for(int p=0;p<rigid_body_collection.rigid_body_particles.Size();p++) tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

    ground=&tests.Add_Ground(friction,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Twisting_Chain
//#####################################################################
void Twisting_Chain(bool apply_twist)
{
    last_frame=120;
    if(apply_twist) last_frame=1200;
//    dynamic_subsampling=true;

    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.cfl=(T)1;
    int n=3;
    RIGID_BODY<TV>* rigid_body=0;
    for(int i=0;i<3*n+1;i++){
        FRAME<TV> frame(TV((T)1.3*(i+1),6,0),ROTATION<TV>((i+1)*(T)pi/2,TV(1,0,0)));
        if(i%3==1){
            rigid_body=&tests.Add_Rigid_Body("torus",1,(T).5);
            rigid_body->Frame()=frame*FRAME<TV>(TV(),ROTATION<TV>((T)pi/2,TV(1,0,0)));
            rigid_body->Set_Mass(30);
            if(i==1 || (i==3*n+1 && !apply_twist)) rigid_body->is_static=true;
            else if(i==3*n+1 && apply_twist){rigid_body->Is_Kinematic()=true;kinematic_id=rigid_body->particle_index;}}
        else if(i%3==2) tests.Make_Lathe_Chain(frame,(T).3);
        else tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",RIGID_BODY_STATE<TV>(frame),true,true,1000,(T).5);}
    if(apply_twist){
        curve.Add_Control_Point(0,rigid_body->Frame());
        for(int i=0;i<=23;i++) curve.Add_Control_Point((T)i+2,rigid_body->Frame()*FRAME<TV>(TV(),ROTATION<TV>((T)(.5*pi*(i+1)),TV(1,0,0))));
        for(int i=0;i<=23;i++) curve.Add_Control_Point(-(T)i+49,rigid_body->Frame()*FRAME<TV>(TV(),ROTATION<TV>((T)(.5*pi*(i+1)),TV(1,0,0))));}

    tests.Add_Ground((T).5,-2,1);
}
//#####################################################################
// Function Brush_And_Wheel
//#####################################################################
void Brush_And_Wheel()
{
    last_frame=120;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
//  solids_parameters.cfl=(T)5;

    RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1,2,48,4);
    cylinder.Frame().t=TV(0,1,0);
    cylinder.Twist().angular=TV(0,0,1);
    cylinder.Update_Angular_Momentum();
    JOINT<TV>* joint=new ANGLE_JOINT<TV>();arb.joint_mesh.Add_Articulation(cylinder.particle_index,int(3),joint);
    joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,0,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));
    joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,3,0),ROTATION<TV>(-(T)pi/2,TV(0,1,0))));

    RIGID_BODY<TV>& box=tests.Add_Analytic_Box(TV(2,2,2));
    box.Frame().t=TV((T)4.5,(T)3.6,0);
    box.Is_Kinematic()=true;
    curve.Add_Control_Point(0,FRAME<TV>(TV((T)4.5,(T)3.6,0)));
    curve.Add_Control_Point((T).5,FRAME<TV>(TV((T)4.5,(T)3.6,0)));
    kinematic_id=box.particle_index;

    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    tests.Create_Mattress(GRID<TV>(VECTOR<int,3>(33,5,5),RANGE<TV>(TV(-1,(T)3.1,-(T).5),TV((T)4.5,(T)4.1,(T).5))),1000);
    tests.Bind_Particles_In_Rigid_Body(box);
    tests.Add_Ground((T).5,-2,1);
}
//#####################################################################
// Function Random_Weights
//#####################################################################
TV Random_Weights(RANDOM_NUMBERS<T>& random)
{
    T x=random.Get_Uniform_Number((T)0,(T)1),y=random.Get_Uniform_Number((T)0,(T)1),s=x+y;
    if(s<=1) return TV(x,y,1-s);
    return TV(1-x,1-y,s-1);
}
//#####################################################################
// Function Ring_Drop
//#####################################################################
void Ring_Drop()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    last_frame=30+(int)(450*num_objects_multiplier);
    solids_parameters.cfl=5;
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
    solids_parameters.triangle_collision_parameters.check_mesh_for_self_intersection=false;
    solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg=false;
    solids_parameters.use_post_cg_constraints=false;

    dynamic_subsampling=true;

    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.enforce_poststabilization_in_cg=true;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=30;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)5;
    solid_body_collection.collision_body_list.Set_Collision_Body_Thickness((T).1);
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-3;
//    solids_parameters.collisions_small_number=(T)1e-5;
//    solids_parameters.self_collision_friction_coefficient=(T)0;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;

    T mu=(T)0.6;
    int poles=5;

    RANDOM_NUMBERS<T> random;
    random.Set_Seed(1236);

    ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
    RANGE<TV> world(TV(-(T)14,(T)30,-(T)14),TV((T)14,(T)(30+(int)(52*num_objects_multiplier)),(T)14));
//    RANGE<TV> world(TV(-(T)5,(T)30,-(T)5),TV((T)5,(T)(30+(int)(300*num_objects_multiplier)),(T)5));

    T g=(T)9.8,h=world.min_corner.y,base_t=sqrt(2*h/g);
    maximum_fall_speed=sqrt(2*g*h);

    int num_objects=(int)(20*num_objects_multiplier);

    // Poles
    for(int i=0;i<poles;i++)
        for(int j=0;j<poles;j++){
            TV position((i-(poles+1)/(T)2)*7,10,(j-(poles+1)/(T)2)*7);
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/medium_cylinder",1,mu);
            rigid_body.name=STRING_UTILITIES::string_sprintf("pole %d %d",i,j);
            rigid_body.is_static=true;
            rigid_body.Frame().t=position;}

    // Rigid spheres
    for(int i=0;i<num_objects;i++){
        T scale=(T)1.75;
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("sphere",scale,mu);
        rigid_body.Update_Bounding_Box();
        rigid_body.Frame()=Find_Placement(random,rigid_body.axis_aligned_bounding_box,bounding_boxes,world,true);
        rigid_body.Set_Mass(ring_mass);
        T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
        if(clamp_time>0) rigid_body_clamp_time.Append(PAIR<RIGID_BODY<TV>*,T>(&rigid_body,clamp_time));}

    // Rigid tori
    for(int i=0;i<num_objects;i++){
        T scale=(T).75;
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("Rings_Test/ring_revolve",scale,mu);

        rigid_body.Update_Bounding_Box();
        rigid_body.Frame()=Find_Placement(random,rigid_body.axis_aligned_bounding_box,bounding_boxes,world,true);
        rigid_body.Set_Mass(ring_mass);
        T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
        if(clamp_time>0) rigid_body_clamp_time.Append(PAIR<RIGID_BODY<TV>*,T>(&rigid_body,clamp_time));}

    // Articulated lathe chains
    for(int i=0;i<num_objects;i++){
        T scale=(T).65;
        RANGE<TV> box(-TV((T)4.5,(T)4.5,(T)1),TV((T)4.5,(T)4.5,(T)1));
        FRAME<TV> frame=Find_Placement(random,box*scale,bounding_boxes,world,true);
        tests.Make_Lathe_Chain(frame,scale,mu);
        int last=solid_body_collection.rigid_body_collection.rigid_body_particles.Size();
        T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
        for(int id=last-5;id<=last;id++){
            solid_body_collection.rigid_body_collection.Rigid_Body(id).Set_Mass(ring_mass/6);
            if(clamp_time>0) rigid_body_clamp_time.Append(PAIR<RIGID_BODY<TV>*,T>(&solid_body_collection.rigid_body_collection.Rigid_Body(id),clamp_time));}}

    // Deformable tori
    TETRAHEDRALIZED_VOLUME<T>* object=0;
    for(int i=0;i<num_objects;i++){
        if(parameter==1 || parameter==2){
            T scale=(T)1;
            RANGE<TV> box(-TV((T)2.5,(T)2.5,(T).5),TV((T)2.5,(T)2.5,(T).5));
            RIGID_BODY_STATE<TV> state(Find_Placement(random,box*scale,bounding_boxes,world,true));
            if(parameter==1) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_1K.tet",state,false,false,1000,scale);
            else if(parameter==2) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_3K.tet",state,false,false,1000,scale);}
        else if(parameter>=3 && parameter<=5){
            T scale=(T)1;
            RANGE<TV> box(-TV((T)3,(T)3,(T)1),TV((T)3,(T)3,(T)1));
            RIGID_BODY_STATE<TV> state(Find_Placement(random,box*scale,bounding_boxes,world,true));
            if(parameter==3) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_3788.tet",state,false,false,1000,scale);
            else if(parameter==4) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_4474.tet",state,false,false,1000,scale);
            else if(parameter==5) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_5428.tet",state,false,false,1000,scale);}
        else if(parameter==6){
            T scale=(T)8;
            RANGE<TV> box(-TV((T)0.375547,(T)0.375547,(T)0.12588),TV((T)0.375547,(T)0.375547,(T)0.12588));
            RIGID_BODY_STATE<TV> state(Find_Placement(random,box*scale,bounding_boxes,world,true));
            object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus.tet",state,false,false,1000,scale);}

        T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
        if(clamp_time>0) structure_clamp_time.Append(PAIR<TETRAHEDRALIZED_VOLUME<T>*,T>(object,clamp_time));
        object->mesh.Initialize_Segment_Mesh();
        SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
        segmented_curve->mesh.elements=object->mesh.segment_mesh->elements;
        deformable_body_collection.Add_Structure(segmented_curve);
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*object,ring_mass/object->Total_Volume());}

    // Deformable spheres
    for(int i=0;i<num_objects;i++){
        T scale=(T)1.75;
        RANGE<TV> box(-TV((T)1,(T)1,(T)1),TV((T)1,(T)1,(T)1));
        RIGID_BODY_STATE<TV> state(Find_Placement(random,box*scale,bounding_boxes,world,true));
        if(parameter==1) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",state,false,false,1000,scale);
        else if(parameter==2) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_2544.tet",state,false,false,1000,scale);
        else if(parameter==3) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_3672.tet",state,false,false,1000,scale);
        else if(parameter==4 || parameter==5 || parameter==6) object=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_4296.tet",state,false,false,1000,scale);
        T d=max((T)0,bounding_boxes.Last().Center().y-world.min_corner.y),clamp_time=d<h?0:base_t+(d-h)/maximum_fall_speed;
        if(clamp_time>0) structure_clamp_time.Append(PAIR<TETRAHEDRALIZED_VOLUME<T>*,T>(object,clamp_time));
        object->mesh.Initialize_Segment_Mesh();
        SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
        segmented_curve->mesh.elements=object->mesh.segment_mesh->elements;
        deformable_body_collection.Add_Structure(segmented_curve);
        SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*object,ring_mass/object->Total_Volume());}

    tests.Add_Ground(mu,(T)0,(T).5,(T)3);
}
//#####################################################################
// Function Cubes_Friction
//#####################################################################
void Cubes_Friction()
{
    last_frame=96;
    RIGID_BODY<TV>* rigid_body=0;

    const char* boxfile="box";
    T box_size=(T)0.05;
    T ground_angle_deg=10;
    ground_angle_rad=ground_angle_deg*((T)pi/180);
    mu=(T).2;
    solids_parameters.cfl=(T)1;

    T initial_speed_down_ramp=(T)1;

    LOG::cout<<"Using angle "<<ground_angle_deg<<" deg ("<<ground_angle_rad<<") (critical coeff="<<tan(ground_angle_rad)<<std::endl;

    FRAME<TV> frame(TV(-box_size*sin(ground_angle_rad),box_size*cos(ground_angle_rad),0),ROTATION<TV>(ground_angle_rad,TV(0,0,1)));
    TWIST<TV> twist(frame.r.Rotate(TV(-initial_speed_down_ramp,0,0)),TV());

    rigid_body=&tests.Add_Rigid_Body(boxfile,box_size,mu);
    rigid_body->Frame()=frame;
    rigid_body->Twist()=twist;
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->name="box";

    frame.t.z+=(T).2;
    RIGID_BODY_STATE<TV> state(frame,twist);
    tests.Create_Mattress(GRID<TV>(VECTOR<int,3>(5,5,5),RANGE<TV>(TV(-(T)1,-(T)1,-(T)1),TV((T)1,(T)1,(T)1))*box_size),false,&state);

    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();solid_body_collection.deformable_body_collection.Add_Structure(&free_particles);
    int particle_rest=solid_body_collection.deformable_body_collection.particles.Add_Element();
    solid_body_collection.deformable_body_collection.particles.mass(particle_rest)=1;
    solid_body_collection.deformable_body_collection.particles.X(particle_rest)=TV(0,0,(T).4);
    free_particles.nodes.Append(particle_rest);

    frame.t.z+=(T).2;
    rigid_body=&tests.Add_Rigid_Body("sphere",(T).1,0);
    rigid_body->Frame().t=frame.t;
    rigid_body->name="analytic solution";
    rigid_body->Is_Kinematic()=true;
    kinematic_id=rigid_body->particle_index;

    RIGID_BODY<TV>& ground=tests.Add_Ground(mu,0,0);
    ground.Set_Coefficient_Of_Rolling_Friction(1);
    ground.Frame().r=ROTATION<TV>(ground_angle_rad,TV(0,0,1));
}
//#####################################################################
// Function PD_Curl
//#####################################################################
void PD_Snake(const T scale,const TV shift,const ROTATION<TV> orient,const T k_p,const int number_of_joints,const T space_adjustment,const T friction=.5)
{
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    RIGID_BODY<TV> *parent_body=0,*child_body=0;
    T cheight=(T)0;
    initial_orientation=ROTATION<TV>((T)half_pi,TV(0,0,1));

    // Create first body
    parent_body=&tests.Add_Rigid_Body("cyllink",scale*(T).2,friction);
    parent_body->Frame().t=shift;
    parent_body->Frame().r=orient*initial_orientation;
    parent_body->name="parent";
    parent_body->Set_Mass(50);

    // Add children and joints
    T desired_x=(T)two_pi/(T)(number_of_joints+1);
    for(int i=0;i<number_of_joints;i++){
        cheight+=scale*((T)1.25+space_adjustment);
        child_body=&tests.Add_Rigid_Body("cyllink",scale*(T).2,friction);
        child_body->Frame().t=orient.Rotate(TV(cheight,0,0))+shift;
        child_body->Frame().r=orient*initial_orientation;
        child_body->Set_Coefficient_Of_Restitution((T)0.5);
        child_body->name=STRING_UTILITIES::string_sprintf("child_%d",i);
        child_body->Set_Mass(50);

        ROTATION<TV> desired_rotation=ROTATION<TV>(desired_x,TV());

        JOINT<TV>* joint=new POINT_JOINT<TV>();arb.joint_mesh.Add_Articulation(child_body->particle_index-1,child_body->particle_index,joint);
        JOINT_FUNCTION<TV>* joint_function=arb.Create_Joint_Function(joint->id_number);
        joint_function->Set_k_p(k_p);joint_function->Set_Target_Angle(desired_rotation);
        joint->Set_Joint_To_Parent_Frame(FRAME<TV>(TV(0,-scale*((T).625+space_adjustment/2),0),ROTATION<TV>(-0*(T)pi/2,TV(0,1,0))));
        joint->Set_Joint_To_Child_Frame(FRAME<TV>(TV(0,scale*((T).625+space_adjustment/2),0),ROTATION<TV>(-0*(T)pi/2,TV(0,1,0))));

        parent_body=child_body;}
}
//#####################################################################
// Function Sidewinding
//#####################################################################
void Sidewinding()
{
    last_frame=1200;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;
    number_of_joints=12;
    TV center(0,(T)1.5,0);
    T scale=(T).55;
    T spacing=(T).25;
    T plank_length=(T)1;
    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);
    ROTATION<TV> snake_rotation(-(T)pi/8,TV(0,1,0));

    joint_separation=(T)1.25*scale;
    wavelength=(T)8*joint_separation;
    period=(T)2;
    height_amplitude=(T).4;
    bend_amplitude=(T)1.4;
    start_time=(T)1.5;
    T friction=(T).8;

    // make the joints
    PD_Snake(scale,center-(T).5*(extent-scale*plank_length)*snake_rotation.Rotate(TV::Axis_Vector(0)),snake_rotation,2000,number_of_joints,0);

    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/snake_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(center,snake_rotation)),false,false,1000,(T).5);

    // binding the deformable particles to the rigid bodies
    for(int p=0;p<rigid_body_collection.rigid_body_particles.Size();p++) tests.Bind_Particles_In_Rigid_Body(rigid_body_collection.Rigid_Body(p));

//     {RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T).3,friction);
//     rigid_body.Frame().t=TV((T)0,(T).3,-10);
//     rigid_body.name="obstruction");
//     rigid_body.Set_Mass(30);}

//     {RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("subdivided_box",(T).3,friction);
//     rigid_body.Frame().t=TV((T)-2,(T).3,-10);
//     rigid_body.name="obstruction");
//     rigid_body.Set_Mass(30);}

    int steps[]={1,2,3,4,5,4,3,3,4,5,6,2,1};
    int count=sizeof steps/sizeof*steps;

    for(int i=0;i<count;i++) for(int j=0;j<6;j++){
        RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("plank",(T)1,friction);
        rigid_body.Frame().t=TV((T)(-(T)35+10*j),(T).2*steps[i-1]-(T).1,(T)(-4-2*i));
        rigid_body.Frame().r=ROTATION<TV>((T)half_pi,TV(0,1,0));
        rigid_body.name="obstruction";
        rigid_body.is_static=true;
        rigid_body.Set_Mass(30);}

    ground=&tests.Add_Ground(friction,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Row_Of_Spheres
//#####################################################################
void Row_Of_Spheres()
{
    last_frame=120;
    RIGID_BODY<TV>* rigid_body=0;
    TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=0;
    INCOMPRESSIBLE_FINITE_VOLUME<TV,3>* fvm=0;
    solids_parameters.cfl=(T)1;

    rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).2);
    rigid_body->Frame().t=TV(-30,0,0);
    rigid_body->Twist().linear=TV(50,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Mass(1000);
    rigid_body->name="in";

    rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).2);
    rigid_body->Frame().t=TV(0,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Mass(1000);
    rigid_body->name="r 1";

    tetrahedralized_volume=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(2,0,0))),true,true,1000);
    fvm=Create_Incompressible_Finite_Volume(*tetrahedralized_volume);
    fvm->mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    solid_body_collection.solid_force_collection.Add_Force(fvm);

    rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).2);
    rigid_body->Frame().t=TV(4,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Mass(1000);
    rigid_body->name="r 3";

    tetrahedralized_volume=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(6,0,0))),true,true,1000);
    fvm=Create_Incompressible_Finite_Volume(*tetrahedralized_volume);
    fvm->mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    solid_body_collection.solid_force_collection.Add_Force(fvm);

    rigid_body=&tests.Add_Rigid_Body("sphere",(T)1,(T).2);
    rigid_body->Frame().t=TV(8,0,0);
    rigid_body->Set_Coefficient_Of_Restitution(0);
    rigid_body->Set_Mass(1000);
    rigid_body->name="out";

    RIGID_BODY<TV>& ground=tests.Add_Ground((T).2,-(T)1,0);
    ground.Set_Coefficient_Of_Rolling_Friction(1);
}
//#####################################################################
// Function Add_Maggot
//#####################################################################
void Add_Maggot(const T scale,const RIGID_BODY_STATE<TV>& state,const std::string& resolution,T friction=(T).2)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    number_of_joints=2;
    T spacing=(T).05;
    T plank_length=(T)1;
    T extent=scale*(number_of_joints*(plank_length+spacing)+plank_length);
    ROTATION<TV> snake_rotation(-(T)pi/8,TV(0,1,0));

    joint_separation=(T)1.05*scale;

    // make the joints
    PD_Snake(scale,state.frame.t-(T).5*(extent-scale*plank_length)*state.frame.r.Rotate(TV::Axis_Vector(0)),state.frame.r,2000*scale,number_of_joints,-(T).2,friction);

    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_"+resolution+"K.tet",state,false,false,1000,(T)1.2*scale);

    // binding the deformable particles to the rigid bodies
    for(int p=rigid_body_collection.rigid_body_particles.Size()-number_of_joints;p<=rigid_body_collection.rigid_body_particles.Size();p++){
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        rigid_body.Twist().angular=state.twist.angular;
        rigid_body.Twist().linear=state.Pointwise_Object_Velocity(rigid_body.Frame().t);
        tests.Bind_Particles_In_Rigid_Body(rigid_body);}
}
//#####################################################################
// Function Maggot
//#####################################################################
void Maggot()
{
    last_frame=360;
    solids_parameters.cfl=(T)1;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    wavelength=(T)8*joint_separation;
    period=(T)1;
    start_time=(T)1.5;

    Add_Maggot((T).65,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1.5,0))),"8");

    ground=&tests.Add_Ground((T).8,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
}
//#####################################################################
// Function Bowl_Of_Maggots
//#####################################################################
void Bowl_Of_Maggots(bool use_blocks)
{
    last_frame=360;
    solids_parameters.cfl=use_blocks?(T).5:(T)1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=15;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)3;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    if(use_blocks) dynamic_subsampling=false;
    refinement_distance=(T).07;
    solids_parameters.enforce_poststabilization_in_cg=false;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    wavelength=(T)8*joint_separation;
    period=(T)1;
    height_amplitude=(T).4;
    bend_amplitude=(T)1.4;
    start_time=(T)0;

    RIGID_BODY<TV>& bowl=tests.Add_Rigid_Body("bowl",(T)5,(T).2);
    bowl.Frame().t=TV(0,(T)0.1,0);
    bowl.Frame().r=ROTATION<TV>(-(T)half_pi,TV(1,0,0));
    bowl.name="bowl";
    bowl.is_static=true;

    random_numbers.Set_Seed(1225);
    ARRAY<ORIENTED_BOX<TV> > bounding_boxes;
    RANGE<TV> maggot_box(-TV((T)1.5,(T).5,(T).5),TV((T)1.5,(T).5,(T).5));
    RANGE<TV> block_box(-TV(1,1,1),TV(1,1,1));
    RANGE<TV> torus_box(-TV((T)1.25,(T).25,(T)1.25),TV((T)1.25,(T).25,(T)1.25));
    RANGE<TV> world(TV(-(T)3.5,(T)2,-(T)3.5),TV((T)3.5,(use_blocks?(T)8:(T)6),(T)3.5));
    T scale=(T).65;
    T block_scale=(T).2;
    T torus_scale=(T).7;
    for(int i=0;i<20;i++){
        FRAME<TV> frame=Find_Placement(random_numbers,maggot_box*scale*(T)1.2,bounding_boxes,world,true);
        periods.Append(random_numbers.Get_Uniform_Number((T).7,(T)1.3));
        phases.Append(random_numbers.Get_Uniform_Number((T)0,(T)two_pi));
        Add_Maggot(scale,RIGID_BODY_STATE<TV>(frame),use_blocks?"3":"3");
        if(use_blocks){
            RIGID_BODY<TV>& block=tests.Add_Rigid_Body("torus",torus_scale,(T).2);
            FRAME<TV> frame=Find_Placement(random_numbers,torus_box*torus_scale,bounding_boxes,world,true);
            block.Frame()=frame;
            block.name="block";
            block.Set_Mass(ring_mass);}}
    (void)block_scale;
    ground=&tests.Add_Ground((T).8,0,0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=false;
}
//#####################################################################
// Function Trapped_Maggot
//#####################################################################
void Trapped_Maggot()
{
    last_frame=360;
    solids_parameters.cfl=(T).5;

    ARTICULATED_RIGID_BODY<TV>& arb=solid_body_collection.rigid_body_collection.articulated_rigid_body;

    arb.Set_Iterative_Tolerance((T)1e-4);
    // prestabilization settings
    arb.Set_Contact_Level_Iterations(5);
    arb.Set_Shock_Propagation_Level_Iterations(5);
    arb.Set_Use_Shock_Propagation(false);
    arb.Set_Do_Final_Pass(false);
    // poststabilization settings
    arb.Set_Poststabilization_Iterations(5);
    arb.poststabilization_projection_iterations=2;

    wavelength=(T)8*joint_separation;
    period=(T)1;
    start_time=(T)5;

    Add_Maggot((T).65,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T).4,0),ROTATION<TV>((T)half_pi,TV(1,0,0)))),"8");

    RIGID_BODY<TV>& block=tests.Add_Rigid_Body("torus",5,(T).6);
    block.Set_Coefficient_Of_Restitution(0);
    block.Frame().t=TV((T)5,(T)1.75,0);
    block.Frame().r=ROTATION<TV>(-(T)pi*(T).025,TV(0,0,1));
    block.name="heavy object";
    block.Set_Mass(1000);

    ground=&tests.Add_Ground((T).8,0,0);
    ground->Set_Coefficient_Of_Restitution(0);

    arb.Use_PD_Actuators();
    arb.global_post_stabilization=true;
    arb.poststabilization_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
}
//#####################################################################
// Function Initialize_Dynamic_Subsampling
//#####################################################################
void Initialize_Dynamic_Subsampling()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    for(int i=0;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++){
        tetrahedralized_volume->Update_Number_Nodes();
        if(!tetrahedralized_volume->triangulated_surface) tetrahedralized_volume->Initialize_Triangulated_Surface();
        surface_elements.Append_Elements(tetrahedralized_volume->triangulated_surface->mesh.elements);}

    for(int i=0;TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
        triangulated_surface->Update_Number_Nodes();
        surface_elements.Append_Elements(triangulated_surface->mesh.elements);}

    surface_elements.Flattened().Get_Unique(surface_particles);

    FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create(deformable_body_collection.particles);
    deformable_body_collection.Add_Structure(&free_particles);

    old_number_particles=deformable_body_collection.particles.Size();
    triangle_free_particles.Resize(surface_elements.m);

//    solids_parameters.write_static_variables_every_frame=true;
}
//#####################################################################
// Function Add_Subsamples
//#####################################################################
void Add_Subsamples(const int surface_triangle_index,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    const VECTOR<int,3>& triangle=surface_elements(surface_triangle_index);
    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    for(int i=0;i<subsamples;i++){
        VECTOR<T,3> weights;do{weights=random_numbers.Get_Uniform_Vector(TV(),TV(1,1,1));}while(weights.x+weights.y>=1);weights.z=(T)1-weights.x-weights.y;
        int hard_bound_particle=particles.Add_Element_From_Deletion_List();
        new_binding_list.Append(new LINEAR_BINDING<TV,3>(particles,hard_bound_particle,triangle,weights));
        int soft_bound_particle=particles.Add_Element_From_Deletion_List();
        new_soft_bindings.Append(VECTOR<int,2>(soft_bound_particle,hard_bound_particle));
        particle_subsamples.Append(soft_bound_particle);}
}
//#####################################################################
// Function Delete_Subsamples
//#####################################################################
void Delete_Subsamples(const int surface_triangle_index)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    // only need to delete the particles
    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    for(int i=0;i<particle_subsamples.m;i++){
        int soft_bound_particle=particle_subsamples(i);
        int hard_bound_particle=soft_bindings.Parent(soft_bound_particle);
        particles.Add_To_Deletion_List(soft_bound_particle);
        particles.Add_To_Deletion_List(hard_bound_particle);}
    particle_subsamples.Remove_All();
}
//#####################################################################
// Function Persist_Subsamples
//#####################################################################
void Persist_Subsamples(const int surface_triangle_index,ARRAY<BINDING<TV>*>& new_binding_list,ARRAY<VECTOR<int,2> >& new_soft_bindings)
{
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    ARRAY<int>& particle_subsamples=triangle_free_particles(surface_triangle_index);
    if(!particle_subsamples.m) return Add_Subsamples(surface_triangle_index,new_binding_list,new_soft_bindings);

    assert(particle_subsamples.m==subsamples);
    for(int i=0;i<particle_subsamples.m;i++){
        int soft_bound_particle=particle_subsamples(i);
        // save soft binding
        const VECTOR<int,2>& soft_binding=soft_bindings.bindings(soft_bindings.Soft_Binding(soft_bound_particle));
        new_soft_bindings.Append(soft_binding);
        // save hard binding
        int hard_binding_index=binding_list.binding_index_from_particle_index(soft_binding.y);
        new_binding_list.Append(binding_list.bindings(hard_binding_index));
        binding_list.bindings(hard_binding_index)=0;} // so that binding won't be deleted when binding_list is rebuilt
}
//#####################################################################
// Function Update_Subsamples
//#####################################################################
void Update_Subsamples()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list=deformable_body_collection.collisions.collision_body_list;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    FREE_PARTICLES<TV>& free_particles=deformable_body_collection.template Find_Structure<FREE_PARTICLES<TV>&>();

    ARRAY<BINDING<TV>*> new_binding_list;
    ARRAY<VECTOR<int,2> > new_soft_bindings;

    LOG::SCOPE("Update_Subsamples");

    if(temporarily_disable_dynamic_subsampling) return;

    binding_list.Clamp_Particles_To_Embedded_Positions();

    // the minimum distance of each particle to a collision object
    particle_distances.Resize(old_number_particles);
    INDIRECT_ARRAY<ARRAY<T>,ARRAY<int>&> subset=particle_distances.Subset(surface_particles);subset.Fill(FLT_MAX);

    for(int i=0;i<surface_particles.m;i++){int p=surface_particles(i);
        for(COLLISION_GEOMETRY_ID body(0);body<collision_body_list.Size();body++)
            particle_distances(p)=PhysBAM::min(particle_distances(p),collision_body_list(body).Implicit_Geometry_Extended_Value(particles.X(p)));}

    // iterate over surface elements
    for(int t=0;t<surface_elements.m;t++){
        const VECTOR<int,3>& triangle=surface_elements(t);
        T triangle_distance=particle_distances.Subset(triangle).Min();
        if(triangle_distance<refinement_distance) Persist_Subsamples(t,new_binding_list,new_soft_bindings);
        else Delete_Subsamples(t);}

    // rebuild free particles
    free_particles.nodes.Remove_All();
    for(int t=0;t<surface_elements.m;t++){
        ARRAY<int>& subsamples=triangle_free_particles(t);
        if(subsamples.m) free_particles.nodes.Append_Elements(subsamples);}

    // rebuild hard bindings
    binding_list.Clean_Memory();
    for(int b=0;b<new_binding_list.m;b++) binding_list.Add_Binding(new_binding_list(b));

    // rebuild soft bindings
    soft_bindings.Clean_Memory();
    for(int b=0;b<new_soft_bindings.m;b++) soft_bindings.Add_Binding(new_soft_bindings(b),true);

    binding_list.Clamp_Particles_To_Embedded_Positions();binding_list.Clamp_Particles_To_Embedded_Velocities();
    soft_bindings.Clamp_Particles_To_Embedded_Positions();soft_bindings.Clamp_Particles_To_Embedded_Velocities();

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();
    
    // correct mass
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // reinitialize deformable object collisions
    deformable_body_collection.collisions.Initialize_Object_Collisions(solids_parameters.deformable_object_collision_parameters.collide_with_interior,
        solids_parameters.deformable_object_collision_parameters.collision_tolerance,
        solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects,
        solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions,
        solids_parameters.deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity);

    solid_body_collection.deformable_body_collection.triangle_repulsions.Set_Repulsion_Thickness(solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness);

    // reinitialize fragments
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Push_Out_Test
//#####################################################################
void Push_Out_Test2()
{
    last_frame=360;
    RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
    RIGID_BODY<TV>& rigid_body2=tests.Add_Rigid_Body("subdivided_box",1,(T)0);
    rigid_body1.Frame().t=TV((T).5,10,0);
    rigid_body2.Frame().t=TV((T)-.5,11,0);
}
//#####################################################################
// Function Sphere_Fall
//#####################################################################
void Sphere_Fall()
{
    last_frame=360;
    solids_parameters.triangle_collision_parameters.total_collision_loops=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=20;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.enforce_repulsions_in_cg=false;
    solids_parameters.enforce_poststabilization_in_cg=true;
    solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=30;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness*=(T)5;
    solid_body_collection.collision_body_list.Set_Collision_Body_Thickness((T).1);
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-3;
//    solids_parameters.collisions_small_number=(T)1e-5;
//    solids_parameters.self_collision_friction_coefficient=(T)0;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;

    bool automatically_add_to_collision_structures=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    TETRAHEDRALIZED_VOLUME<T>& object=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
    object.mesh.Initialize_Segment_Mesh();
    SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
    segmented_curve->mesh.elements=object.mesh.segment_mesh->elements;
    deformable_body_collection.Add_Structure(segmented_curve);
    tests.Add_Ground();

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.particles.mass*=100;
}
//#####################################################################
// Function Sphere_Block
//#####################################################################
void Sphere_Block(bool deformable_on_top)
{
    last_frame=120;
    solids_parameters.cfl=2;
    solids_parameters.triangle_collision_parameters.perform_self_collision=true;
    const char* sphere="/Tetrahedralized_Volumes/sphere_coarse.tet";

    tests.Create_Tetrahedralized_Volume(data_directory+sphere,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)(deformable_on_top?6:3),0))),true,true,1000);
    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("sphere",1,(T).5);
    rigid_body->Set_Mass(200);
    rigid_body->Frame().t=TV(0,(T)(deformable_on_top?3:6),0);
    rigid_body->Set_Coefficient_Of_Restitution((T)0.5);
    rigid_body->name="box";
    tests.Add_Ground();
}
//#####################################################################
// Function Two_Way_Tori
//#####################################################################
void Two_Way_Tori()
{
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    last_frame=240;
    solids_parameters.cfl=2;
    EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=tests.Create_Embedded_Tetrahedralized_Volume(TORUS<T>(TV(),TV(0,0,1),(T).3,(T).6),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true);
    embedding.Update_Binding_List_From_Embedding(solid_body_collection.deformable_body_collection,false);
    tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,soft_bindings);
    embedding.Update_Number_Nodes();
    tests.Initialize_Tetrahedron_Collisions(1,embedding.embedded_object.simplicial_object,solids_parameters.triangle_collision_parameters,&embedding.material_surface);
    RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("torus",1,(T).5);
    rigid_body->Set_Mass(1000);
    rigid_body->Frame().t=TV(0,(T)2,0);
    rigid_body->Frame().r=ROTATION<TV>((T)half_pi,TV(0,0,1));
    tests.Add_Ground();
    if(use_forces_for_drift){
        soft_bindings.use_impulses_for_collisions.Fill(false);
        soft_bindings.Initialize_Binding_Mesh();
        solid_body_collection.solid_force_collection.Add_Force(Create_Edge_Binding_Springs(solid_body_collection.deformable_body_collection.particles,*soft_bindings.binding_mesh,(T)1e6,(T)1));}
}
//#####################################################################
};
}
#endif
