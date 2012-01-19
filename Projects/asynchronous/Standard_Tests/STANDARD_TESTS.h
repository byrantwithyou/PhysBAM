//#####################################################################
// Copyright 2008, Nipun Kwatra, Craig Schroeder, Wen Zheng.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Sphere free falling to the ground
//   2. Constrained Sphere
//   3. Cloth Curtain and Ball
//   4. Constrained Cloth
//   5. A single hair
//   6. A single hair with control
//   7. Cloth and spinning sphere
//   8. Wardrobe cloth test
//   9. Cloth constrained on left side and pulled on other sides
//  10. Soft bound copy of sphere surface
//  11. Soft bound copy of sphere surface using stored data
//  12. Copy of sphere surface attached to one ring of tet boundary
//  13. Copy of sphere surface attached to one ring of tet boundary using stored data
//  15. Same as test 14 but two-way coupled asynchronous 
//  16. Comparison of without, with boundary only, and with all particles to be asynchronous
//  17. Comparison of asynchronous layered boxes
//  18. Armadillo
//  19. Analytic gravity test
//  20. Analytic ether drag test
//  21. Analytic spring test
//  22. Analytic one spring test
//  23. Adaptive Asynchronous
//  24. Asynchronous projection
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_UTILITIES.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Polynomials/QUADRATIC.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/TETRAHEDRON.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Implicit_Objects/ANALYTIC_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/AXIAL_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/DEFORMABLE_GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_TET_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/ANGLE_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/CONSTRAINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/JOINT_FUNCTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/POINT_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/PRISMATIC_TWIST_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Joints/RIGID_JOINT.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/ASYNCHRONOUS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/SOLIDS_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_DRIVER_UNIFORM.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T_input,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef VECTOR<int,2> Vi2;typedef VECTOR<int,3> Vi3;typedef VECTOR<int,4> Vi4;
    typedef typename TV::SPIN T_SPIN;
public:
    SOLIDS_STANDARD_TESTS<TV> tests;

    int parameter;

    RIGID_BODY_COLLECTION<TV>& rigid_body_collection;

    T damping_multiplier;
    T stiffness_multiplier;
    T overdamping_fraction;
    T bending_stiffness_multiplier;
    T bending_damping_multiplier;
    T axial_stiffness_multiplier;
    T axial_damping_multiplier;
    T soft_bindings_multiplier;
    T soft_surface_multiplier;
    T soft_bound_edge_damping;
    T soft_binding_damping;
    RIGID_BODY<TV>* ground;
    int number_side_panels;
    T aspect_ratio,side_length;
    int cloth_triangles;
    bool fully_implicit,coarse_fully_implicit,fine_fully_implicit;
    ARRAY<int> constrained_point;
    int animated_particle;
    INTERPOLATION_CURVE<T,TV> animation_curve;
    T hair_stiffness,hair_damping;

    // test 8
    T test_8_constrained_off;
    T test_8_friction_off;
    T test_8_wind_off;

    // Test 9
    int test_9_skipped_particles_number;

    // Test 10, 11, 12, 13
    ARRAY<INTERPOLATION_CURVE<T,TV> > saved_V,saved_X;
    int sim_length;
    T sim_switch_time;
    T subtract_time;
    HASHTABLE<int,int> particle_map; // Maps particle to copy particle
    int num_controlled_particles;
    int number_of_spheres;
    ARRAY<SEGMENTED_CURVE<TV>*> new_boundary_segmented_curves;
    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> primary_tetrahedralized_volumes;
    ARRAY<SEGMENTED_CURVE<TV>*> primary_segmented_curves;
    ARRAY<TRIANGULATED_SURFACE<T>*> new_boundary_surfaces;

    // Test 12, 13
    ARRAY<SEGMENTED_CURVE<TV>*> boundary_one_ring_segmented_curves;
    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> boundary_tetrahedralized_volumes;

    // Test 14
    SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >* driver;
    int number_of_rings;
    bool invert_asynchronous_implicit;
    int type_of_asynchronous_configuration;
    T substeps_per_frame;

    // Test 15, 16, 17
    ASYNCHRONOUS_EVOLUTION<TV>* asynchronous_evolution;
    T projection_rigidity;
    ARRAY<TETRAHEDRALIZED_VOLUME<T>*> complementary_boundary_tetrahedralized_volumes;
    ARRAY<int> gravity_particles;
    T height_offset;
    T coverage_percent;
    bool split_forces_into_fine_coarse;
    int model_num;
    bool use_async;

    ARRAY<int> int_lists[3];

    // Test 17
    ARRAY<SEGMENTED_CURVE<TV>*> complementary_boundary_segmented_curves;

    // analytic
    T mass_analytic,stiffness_analytic,damping_analytic,restlength_analytic,initial_displacement_analytic,initial_velocity_analytic,orthogonal_velocity;
    bool use_orthogonal_velocity;
    int mixed_particle_index,implicit_particle_index,finescale_particle_index;
    TV initial_X_mixed_particle,initial_X_implicit_particle,initial_X_finescale_particle;
    
    // test 24
    bool treat_bottom_async,treat_left_async,test_implicit_in_explicit_out;
    T make_it_bad;

    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::output_directory;using BASE::restart;using BASE::frame_rate;
    using BASE::stream_type;using BASE::Time_At_Frame;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),parameter(1),
        rigid_body_collection(solid_body_collection.rigid_body_collection),animated_particle(0),
        driver(0),asynchronous_evolution(0)
    {
        fluids_parameters.simulate=false;
        solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;

        hair_stiffness=(T)100;
        hair_damping=(T)10;

        solids_parameters.cfl=4;
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        //solids_parameters.allow_intersections=true;
        //solids_parameters.allow_intersections_tolerance=10000;
        solids_parameters.triangle_collision_parameters.output_interaction_pairs=true;
        solids_parameters.rigid_body_collision_parameters.use_push_out=true;
        solids_parameters.use_rigid_deformable_contact=true;
        solids_parameters.rigid_body_collision_parameters.collision_bounding_box_thickness=(T)1e-3;
        solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
        solids_parameters.verbose_dt=true;
        solids_parameters.triangle_collision_parameters.repulsion_pair_attractions_threshold=-(T).3;
        solids_parameters.enforce_repulsions_in_cg=false;
        solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.implicit_solve_parameters.cg_iterations=10000;
        solids_parameters.implicit_solve_parameters.cg_restart_iterations=100;
        solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > F,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Advance_One_Time_Step_End_Callback(const T dt,const T time) {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Option_Argument("-sp","Use spatial partition for collision body collisions");
    parse_args->Add_Double_Argument("-clothcfl",4.,"Cloth CFL");
    parse_args->Add_Double_Argument("-cgsolids",1e-3,"CG tolerance for backward Euler");
    parse_args->Add_Double_Argument("-substeps_per_frame",(T)0,"set the async sim to take dt=1/(substeps_per_frame*frame_rate).");
    parse_args->Add_Double_Argument("-overdamping_fraction",(T)3,"over damping for edge and altitude springs");
    parse_args->Add_Double_Argument("-bend_damping",(T)1,"multiplier on bending spring damping.");
    parse_args->Add_Double_Argument("-bend_stiffness",(T)1,"multiplier on bending spring stiffness.");
    parse_args->Add_Double_Argument("-axial_damping",(T)1,"multiplier on axial bending spring damping.");
    parse_args->Add_Double_Argument("-axial_stiffness",(T)1,"multiplier on axial bending spring stiffness.");
    parse_args->Add_Double_Argument("-edge_damping",(T)1,"multiplier on edge spring damping.");
    parse_args->Add_Double_Argument("-edge_stiffness",(T)1,"multiplier on edge spring stiffness.");
    parse_args->Add_Double_Argument("-panel_multiplier",(T)1,"scale resolution.");
    parse_args->Add_Option_Argument("-fully_implicit","use fully implicit forces");
    parse_args->Add_Option_Argument("-no_coarse_fully_implicit","use fully implicit coarse forces");
    parse_args->Add_Option_Argument("-fine_fully_implicit","use fully implicit fine forces");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Option_Argument("-use_trapezoid","use trapezoid rule instead of backward Euler for velocity update");
    parse_args->Add_Double_Argument("-soft_bindings_multiplier",(T)1,"scale binding stiffness.");
    parse_args->Add_Double_Argument("-soft_surface_multiplier",(T)1,"scale surface stiffness.");
    parse_args->Add_Double_Argument("-soft_binding_damping",(T)1,"scale binding damping.");
    parse_args->Add_Double_Argument("-soft_bound_edge_damping",(T)1,"scale soft bound surface damping.");
    parse_args->Add_Integer_Argument("-number_of_spheres",1,"used by tests 10,11,12,13,14");
    parse_args->Add_Integer_Argument("-number_of_rings",1,"used by tests 14, negative means all rings");
    parse_args->Add_Option_Argument("-invert","invert asynchronous and implicit elements in test 14.");
    parse_args->Add_Integer_Argument("-type_of_config",1,"used by tests 14. 1: boundary; 2: left half; 3: top half.");
    parse_args->Add_Double_Argument("-height_offset",0,"offset sphere height in tests 15,16.");
    parse_args->Add_Option_Argument("-orthogonal_velocity","give an orthogonal velocity in test 21");
    parse_args->Add_Double_Argument("-coverage_percent",1,"the percentage of areas that covered by asyn boundary in tests 15,16.");    
    parse_args->Add_Option_Argument("-async","to use asynchronous or not");
    parse_args->Add_Integer_Argument("-model_num",1,"input model file");
    parse_args->Add_Double_Argument("-projection_rigidity",(T).5,"parameter for asynchronous evolution. 1 projects to fully rigid.");
    parse_args->Add_Option_Argument("-imp_out_exp_in_test","test implicit boundary and explicit inside");    
    parse_args->Add_Option_Argument("-bottom_async","treat bottom asynchronously");    
    parse_args->Add_Option_Argument("-left_async","treat left side asynchronously");    
    parse_args->Add_Double_Argument("-make_it_bad",1,"create a really short edge, and give how bad it is");    
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    if(parse_args->Is_Value_Set("-parameter")){parameter=parse_args->Get_Integer_Value("-parameter");output_directory+=STRING_UTILITIES::string_sprintf("_param%i",parameter);}
    if(parse_args->Is_Value_Set("-sp")){solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects=true;}
    substeps_per_frame=(T)parse_args->Get_Double_Value("-substeps_per_frame");
    fully_implicit=parse_args->Is_Value_Set("-fully_implicit");
    fine_fully_implicit=parse_args->Is_Value_Set("-fine_fully_implicit");
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
    solids_parameters.use_trapezoidal_rule_for_velocities=parse_args->Get_Option_Value("-use_trapezoid");
    overdamping_fraction=(T)parse_args->Get_Double_Value("-overdamping_fraction");
    bending_stiffness_multiplier=(T)parse_args->Get_Double_Value("-bend_stiffness");
    bending_damping_multiplier=(T)parse_args->Get_Double_Value("-bend_damping");
    axial_stiffness_multiplier=(T)parse_args->Get_Double_Value("-axial_stiffness");
    axial_damping_multiplier=(T)parse_args->Get_Double_Value("-axial_damping");
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-edge_stiffness");
    damping_multiplier=(T)parse_args->Get_Double_Value("-edge_damping");
    T panel_multiplier=(T)parse_args->Get_Double_Value("-panel_multiplier");
    soft_bindings_multiplier=(T)parse_args->Get_Double_Value("-soft_bindings_multiplier");
    soft_surface_multiplier=(T)parse_args->Get_Double_Value("-soft_surface_multiplier");
    soft_binding_damping=(T)parse_args->Get_Double_Value("-soft_binding_damping");
    soft_bound_edge_damping=(T)parse_args->Get_Double_Value("-soft_bound_edge_damping");
    number_of_spheres=parse_args->Get_Integer_Value("-number_of_spheres");
    number_of_rings=parse_args->Get_Integer_Value("-number_of_rings");
    invert_asynchronous_implicit=parse_args->Get_Option_Value("-invert");
    type_of_asynchronous_configuration=parse_args->Get_Integer_Value("-type_of_config");
    height_offset=(T)parse_args->Get_Double_Value("-height_offset");
    use_orthogonal_velocity=parse_args->Get_Option_Value("-orthogonal_velocity");
    coverage_percent=(T)parse_args->Get_Double_Value("-coverage_percent");coverage_percent=clamp(coverage_percent,(T)0,(T)1);
    model_num=parse_args->Get_Integer_Value("-model_num");
    use_async=parse_args->Is_Value_Set("-async");
    projection_rigidity=(T)parse_args->Get_Double_Value("-projection_rigidity");
    output_directory+=STRING_UTILITIES::string_sprintf("_prigid_%f",projection_rigidity);
    split_forces_into_fine_coarse=true;

    test_implicit_in_explicit_out=parse_args->Get_Option_Value("-imp_out_exp_in_test");
    treat_bottom_async=parse_args->Get_Option_Value("-bottom_async");
    treat_left_async=parse_args->Get_Option_Value("-left_async");
    make_it_bad=(T)parse_args->Get_Double_Value("-make_it_bad");make_it_bad=clamp(make_it_bad,(T)0,(T)1);
    if(parse_args->Is_Value_Set("-no_coarse_fully_implicit")) coarse_fully_implicit=false;
    else coarse_fully_implicit=true;
    number_side_panels=(int)(40*panel_multiplier);
    aspect_ratio=(T)1.7;
    side_length=(T)1.0;
    if(parse_args->Is_Value_Set("-cgsolids")) solids_parameters.implicit_solve_parameters.cg_tolerance=(T)parse_args->Get_Double_Value("-cgsolids");
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Initialize_Volume_Forces
//#####################################################################
void Initialize_Volume_Forces(T stiffness,T damping,ARRAY<TETRAHEDRALIZED_VOLUME<T>*>& structures_for_forces,int force_type)
{
    DEFORMABLES_FORCES<TV>* force=0;
    for(int i=0;i<structures_for_forces.m;i++){
        TETRAHEDRALIZED_VOLUME<T>* structure=structures_for_forces(i);
        if(!structure->mesh.elements.m) continue;
        force=Create_Altitude_Springs(*structure,(T)stiffness/(1+sqrt((T)2)),(T)damping);
        if(force_type){
            ARRAY<int> affected_particle_indices,affected_rigid_body_particle_indices;
            structure->mesh.elements.Flattened().Get_Unique(affected_particle_indices);
            if(force_type==1) asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,false);
            else if(force_type==2) asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);
            else PHYSBAM_FATAL_ERROR("Unrecognized force type for force initialization");}
        else solid_body_collection.Add_Force(force);}
}
//#####################################################################
// Function Initialize_Edge_Forces
//#####################################################################
void Initialize_Edge_Forces(T stiffness,T damping,ARRAY<SEGMENTED_CURVE<TV>*>& structures_for_forces,int force_type)
{
    DEFORMABLES_FORCES<TV>* force=0;
    for(int i=0;i<structures_for_forces.m;i++){
        SEGMENTED_CURVE<TV>* structure=structures_for_forces(i);
        if(!structure->mesh.elements.m) continue;
        force=Create_Edge_Springs(*structure,(T)stiffness/(1+sqrt((T)2)),(T)damping);
        if(force_type){
            ARRAY<int> affected_particle_indices,affected_rigid_body_particle_indices;
            structure->mesh.elements.Flattened().Get_Unique(affected_particle_indices);
            if(force_type==1) asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,true);
            else if(force_type==2) asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);
            else PHYSBAM_FATAL_ERROR("Unrecognized force type for force initialization");}
        else solid_body_collection.Add_Force(force);}
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
        case 1: Sphere_Fall();break;
        case 2: Constrained_Sphere();break;
        case 3: Cloth_Ball();break;
        case 4: Constrained_Cloth();break;
        case 5:
        case 6: Single_Hair();break;
        case 7: Cloth_And_Spinning_Sphere();break;
        case 8: Wardrobe_Cloth();break;
        case 9: Cloth_Stretch();break;
        case 10: 
        case 11:
        case 12:
        case 13: Soft_Bound_Sphere_Surface();break;
        case 14:
        case 15: 
        case 16: Asynchronous_Sphere();break;
        case 17: Asynchronous_Layered_Box();break;
        case 18: Asynchronous_Armadillo();break;
        case 19: Gravity_Test();break;
        case 20: Ether_Drag_Test();break;
        case 21: Spring_Test();break;
        case 22: One_Spring_Test();break;
        case 23: Adaptive_Asynchronous();break;
        case 24: Asynchronous_Projected_Sphere();break;
      default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized test number %d",test_number));}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();
    
    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();
        
    // Add forces
    // test_number<->force: swtich table
    // ------------------------- 0 | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 | 11 | 12 | 13 | 14 | 15 | 16 | 17 | 18 | 19 | 20 | 21 | 22 | 23 | 24 |
    bool volume_force_switch[]= {0,  1,  1,  1,  1,  0,  0,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0,   0,   0,   0,   0};
    bool curve_force_switch[]=  {0,  1,  1,  1,  1,  0,  0,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   0,   1,   0,   0,   0,   0,   0,   0};
    bool surface_force_switch[]={0,  0,  0,  1,  1,  0,  0,  1,  1,  1,   0,   0,   0,   0,   0,   0,   0,   0,   1,   0,   0,   0,   0,   0,   0};
    bool gravity_force_switch[]={0,  1,  1,  1,  1,  1,  1,  1,  1,  1,   1,   1,   1,   1,   1,   1,   1,   1,   0,   0,   0,   0,   0,   0,   0};

    // add volume and edge forces
    if(volume_force_switch[test_number] && curve_force_switch[test_number]){
        // setup stiffness
        T stiffness=(T)2e6;
        if(test_number==1) stiffness=(T)1e4;
        if(test_number==10 || test_number==11 || test_number==12 || test_number==13 || test_number==14 || test_number==15  || test_number==16 || test_number==24) stiffness=(T)1e4;
        if(test_number==2) stiffness=(T)1e4;
        if(test_number==18) stiffness=(T)4e5;
        stiffness*=stiffness_multiplier;
        
        // get volumes and curves
        if(!primary_tetrahedralized_volumes.m) for(int i=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(i);i++)
                primary_tetrahedralized_volumes.Append(tetrahedralized_volume);
        ARRAY<TETRAHEDRALIZED_VOLUME<T>*>& tetrahedralized_volumes_for_forces=(test_number==15 || test_number==16)?complementary_boundary_tetrahedralized_volumes:primary_tetrahedralized_volumes;
        if(!primary_segmented_curves.m) for(int i=1;SEGMENTED_CURVE<TV>* segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>*>(i);i++)
                primary_segmented_curves.Append(segmented_curve);
        if(!primary_segmented_curves.m) for(int i=0;i<tetrahedralized_volumes_for_forces.m;i++){
                TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=tetrahedralized_volumes_for_forces(i);
                tetrahedralized_volume->mesh.Initialize_Segment_Mesh();
                SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
                segmented_curve->mesh.elements=tetrahedralized_volume->mesh.segment_mesh->elements;
                segmented_curve->Update_Number_Nodes();
                deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
                primary_segmented_curves.Append(segmented_curve);}
        ARRAY<SEGMENTED_CURVE<TV>*>& segmented_curves_for_forces=(test_number==15 || test_number==16)?complementary_boundary_segmented_curves:primary_segmented_curves;

        int force_type=(asynchronous_evolution && split_forces_into_fine_coarse);
        Initialize_Volume_Forces(stiffness,overdamping_fraction,tetrahedralized_volumes_for_forces,force_type);
        Initialize_Edge_Forces(stiffness,overdamping_fraction,segmented_curves_for_forces,force_type);}

    // add surface forces
    if(surface_force_switch[test_number]){
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            T axial_stiffness=axial_stiffness_multiplier*2/(1+sqrt((T)2)),axial_damping=axial_damping_multiplier*8;
            for(int i=1;TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>(i);i++){
                solid_body_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping));
                solid_body_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));
                solid_body_collection.Add_Force(Create_Axial_Bending_Springs(*triangulated_surface,(T).01,axial_stiffness,axial_damping));
                if(test_number==8){
                    WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(*triangulated_surface,solid_body_collection.rigid_body_collection);
                    solid_body_collection.Add_Force(drag);
                    drag->Use_Linear_Normal_Viscosity((T).001);drag->Use_Constant_Wind(0,TV((T).001,(T).0001,(T).001));}
                PHYSBAM_DEBUG_PRINT("Spring stiffnesses",linear_stiffness,linear_damping,bending_stiffness,bending_damping);}}

    // add gravity 
    if(gravity_force_switch[test_number]){
        if(asynchronous_evolution && split_forces_into_fine_coarse && gravity_particles.m){
            ARRAY<int> affected_rigid_body_particle_indices;
            if(asynchronous_evolution->use_projection){
                // When using projection add gravity as finescale
                asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,
                            &gravity_particles),gravity_particles,affected_rigid_body_particle_indices,coarse_fully_implicit);}
            else{
                asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,
                            &gravity_particles),gravity_particles,affected_rigid_body_particle_indices,coarse_fully_implicit,false);}}
        else tests.Add_Gravity();}

    // Set fully_implicit flag
    if(!(asynchronous_evolution && split_forces_into_fine_coarse) && fully_implicit) 
        for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;

    // add collision structures
    if(test_number==10 || test_number==11 || test_number==12 || test_number==13 || test_number==14 || test_number==15 || test_number==16 
            || test_number==17 || test_number==23 || test_number==24){
        for(int i=0;i<primary_tetrahedralized_volumes.m;i++){
            deformable_body_collection.collisions.collision_structures.Append(primary_tetrahedralized_volumes(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision){
                LOG::cout<<"Adding structure "<<i<<std::endl;
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(primary_tetrahedralized_volumes(i));}}}
    else{
        for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(!dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.deformable_geometry.structures(i))){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision && (!dynamic_cast<FREE_PARTICLES<TV>*>(deformable_body_collection.deformable_geometry.structures(i))))
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}}

    // correct mass
    //binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();

    // initialize
    if(asynchronous_evolution) asynchronous_evolution->Initialize();
    if(test_number==15 || test_number==16) Initialize_Sphere_Analytic_Test();

    // disable strain rate CFL for all forces
    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->limit_time_step_by_strain_rate=false;
}
//#####################################################################
// Function Cloth_Balll
//#####################################################################
void Cloth_Ball()
{
    last_frame=300;
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=1;
    solids_evolution->solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    cloth_triangles=2*number_side_panels*(int)(number_side_panels*aspect_ratio);
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    tests.Add_Ground();
    RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T).25,(T)0);
    body.X().z=(T).5;
    body.Is_Kinematic()=true;
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    constrained_point.Append(1);
    constrained_point.Append(1+m*(n-1));
}
//#####################################################################
// Function Constrained_Cloth
//#####################################################################
void Constrained_Cloth()
{
    last_frame=300;
    cloth_triangles=2*number_side_panels*(int)(number_side_panels*aspect_ratio);
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
    int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
    constrained_point.Append(1);
    constrained_point.Append(1+m*(n-1));
    animated_particle=1311;
    TV initial_point((T)1.7,(T).5,(T).45);
    animation_curve.Add_Control_Point(0,initial_point);
    animation_curve.Add_Control_Point(2,initial_point);
    animation_curve.Add_Control_Point((T)3.3,initial_point-TV(0,10,0));
    animation_curve.Add_Control_Point((T)3.6,initial_point+TV(0,10,0));
    animation_curve.Add_Control_Point((T)3.8,initial_point-TV(0,10,0));
    animation_curve.Add_Control_Point(4,initial_point);
    animation_curve.Add_Control_Point(6,initial_point);
    animation_curve.Add_Control_Point(7,initial_point+TV(0,0,12));
    animation_curve.Add_Control_Point((T)7.5,initial_point-TV(0,0,12));
    animation_curve.Add_Control_Point((T)7.7,initial_point+TV(0,0,12));
    animation_curve.Add_Control_Point((T)7.9,initial_point-TV(0,0,12));
    animation_curve.Add_Control_Point(8,initial_point);
    animation_curve.Add_Control_Point(10,initial_point);
}
//#####################################################################
// Function Cloth_And_Spinning_Sphere
//#####################################################################
void Cloth_And_Spinning_Sphere()
{
    frame_rate=60;
    last_frame=(int)(20*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    aspect_ratio=(T)1.0;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=(T).03;
    solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions=false;

    tests.Create_Cloth_Panel(number_side_panels,(T)1.9*side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).56,0))));
    tests.Add_Ground((T).1);
    RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T).25,(T)1);
    sphere_body.X()=TV(0,(T).30,0);
    sphere_body.Is_Kinematic()=true;
    RIGID_BODY<TV>& body=tests.Add_Rigid_Body("cut_pyramid",(T).1,(T).1);
    body.X()=TV((T)-.65,(T).05,(T).65);body.Rotation()=ROTATION<TV>((T)-pi/2,TV(1,0,0))*body.Rotation();
    body.is_static=true;
}
//#####################################################################
// Function Wardrobe_Cloth
//#####################################################################
void Wardrobe_Cloth()
{
    frame_rate=120;
    last_frame=(int)(20*frame_rate);
    solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
    solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    solids_parameters.implicit_solve_parameters.cg_iterations=200;
    test_8_constrained_off=(T).3;
    test_8_friction_off=(T)3;
    test_8_wind_off=(T)1000;

    tests.Create_Cloth_Panel(number_side_panels,(T)4,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4,(T)-.25))));
    tests.Add_Ground((T)0.2,(T).1);
    RIGID_BODY<TV>& body=tests.Add_Rigid_Body("wardrobe",(T).25,(T)10);
    body.X()=TV(0,(T)1.5,(T).9-(T)0.194);body.Rotation()=ROTATION<TV>((T)-pi/2,TV(1,0,0));
    body.is_static=true;
}
//#####################################################################
// Function Sphere_Fall
//#####################################################################
void Sphere_Fall()
{
    last_frame=300;
    // param 2 = asynchronous with projections

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    solids_parameters.triangle_collision_parameters.perform_self_collision=number_of_spheres>=2;

    for(int i=0;i<number_of_spheres;i++){
        primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)i,(T)(2*i+1),0))),true,true,1000));}

    tests.Add_Ground();

    if(parameter==2){
        asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);
        solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
        split_forces_into_fine_coarse=false;
        
        for(int i=0;i<number_of_spheres;i++){
            TETRAHEDRALIZED_VOLUME<T>& volume=*primary_tetrahedralized_volumes(i);
            volume.mesh.Initialize_Node_On_Boundary();
            ARRAY<int> blob_particles;
            for(int i=0;i<deformable_body_collection.particles.array_collection->Size();i++) if(!(*volume.mesh.node_on_boundary)(i)) blob_particles.Append(i);
            asynchronous_evolution->Add_Blob_From_Particles(blob_particles);}}
}
//#####################################################################
// Function Constrained_Sphere
//#####################################################################
void Constrained_Sphere()
{
    last_frame=300;
    tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
    constrained_point.Append(1407);
    constrained_point.Append(2510);
    constrained_point.Append(1417);
    constrained_point.Append(4271);
    constrained_point.Append(4266);
    animated_particle=1901;
    animation_curve.Add_Control_Point(0,TV(0,2,0));
    animation_curve.Add_Control_Point(2,TV(0,2,0));
    animation_curve.Add_Control_Point((T)3.3,TV(0,-10,0));
    animation_curve.Add_Control_Point((T)3.6,TV(0,10,0));
    animation_curve.Add_Control_Point((T)3.8,TV(0,-10,0));
    animation_curve.Add_Control_Point(4,TV(0,2,0));
    animation_curve.Add_Control_Point(6,TV(0,2,0));
    animation_curve.Add_Control_Point(7,TV(0,2,12));
    animation_curve.Add_Control_Point((T)7.5,TV(0,2,-12));
    animation_curve.Add_Control_Point((T)7.7,TV(0,2,12));
    animation_curve.Add_Control_Point((T)7.9,TV(0,2,-12));
    animation_curve.Add_Control_Point(8,TV(0,2,0));
    animation_curve.Add_Control_Point(10,TV(0,2,0));
}
//#####################################################################
// Function Single_Hair
//#####################################################################
void Single_Hair()
{
    last_frame=300;
    solids_parameters.cfl=4;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
    solids_parameters.implicit_solve_parameters.cg_iterations=1000;
    solids_parameters.implicit_solve_parameters.cg_restart_iterations=100;
    //solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;

    int n_segments=100;
    T perturb_threshold=(T)1e-4;
    ARRAY<TV> points;
    ARRAY<bool> segment_perturbed;
    for(int i=0;i<n_segments+2;i++) points.Append(TV(T(T(i-1)/T(n_segments)),0,0));
    segment_perturbed.Append(false);
    for(int i=1;i<points.m;i++){
        if(i>1 && TRIANGLE_3D<T>(points(i-1),points(i),points(i+1)).Area()<perturb_threshold) segment_perturbed.Append(true);
        else if(i<points.m-1 && TRIANGLE_3D<T>(points(i),points(i+1),points(i+2)).Area()<perturb_threshold) segment_perturbed.Append(true);
        else segment_perturbed.Append(false);}

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>& edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& extra_edges =*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& bending_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    SEGMENTED_CURVE<TV>& torsion_edges=*SEGMENTED_CURVE<TV>::Create(particles);
    TRIANGULATED_SURFACE<T>& triangles=*TRIANGULATED_SURFACE<T>::Create(particles);
    TETRAHEDRALIZED_VOLUME<T>& volume=*TETRAHEDRALIZED_VOLUME<T>::Create(particles);

    int i=2,e0,r0,e1,r1=0,r2=0,r3;
    int last_particle=particles.array_collection->Add_Element();
    particles.X(last_particle)=points(1);
    ARRAY<ARRAY<int> > previous;
    previous.Append(ARRAY<int>());
    TV perturb=(points(1)-points(2)).Orthogonal_Vector();

    int perturb_amount=1;
    while(i<=points.m){
        if(segment_perturbed(i)){
            TV segment_direction=(points(i-1)-points(i)).Normalized();
            perturb=ROTATION<TV>((T)pi/2,segment_direction).Rotate(perturb);
            perturb.Project_Orthogonal_To_Unit_Direction(segment_direction);

            e0=particles.array_collection->Add_Element();
            r0=particles.array_collection->Add_Element();
            particles.X(e0)=(points(i-1)+points(i))*(T).5+perturb*(T)perturb_amount;
            particles.X(r0)=points(i);
            edges.mesh.elements.Append(Vi2(last_particle,r0));
            extra_edges.mesh.elements.Append(Vi2(last_particle,e0));
            extra_edges.mesh.elements.Append(Vi2(e0,r0));
        
            if(previous(i-1).m==3){
                r2=previous(i-1)(1);
                e1=previous(i-1)(2);
                r1=previous(i-1)(3);
                torsion_edges.mesh.elements.Append(Vi2(r2,e0));
                volume.mesh.elements.Append(Vi4(r2,e1,r1,e0));
                torsion_edges.mesh.elements.Append(Vi2(e1,r0));
                volume.mesh.elements.Append(Vi4(e1,r1,r0,e0));
                triangles.mesh.elements.Append(Vi3(r1,e0,r0));
                bending_edges.mesh.elements.Append(Vi2(e1,e0));}
            else if(previous(i-1).m==2){
                r2=previous(i-1)(1);
                r1=previous(i-1)(2);
                torsion_edges.mesh.elements.Append(Vi2(e0,r2));
                volume.mesh.elements.Append(Vi4(r2,r1,e0,r0));
                bending_edges.mesh.elements.Append(Vi2(last_particle,r0));
                if(i-2>=0 && previous(i-2).m==2){
                    r3=previous(i-2)(1);
                    torsion_edges.mesh.elements.Append(Vi2(r3,r0));
                    volume.mesh.elements.Append(Vi4(r3,r2,r1,r0));}
                triangles.mesh.elements.Append(Vi3(r1,e0,r0));
                triangles.mesh.elements.Append(Vi3(r2,r1,r0));}
            
            int previous_index=previous.Append(ARRAY<int>());
            previous(previous_index).Append(last_particle);
            previous(previous_index).Append(e0);
            previous(previous_index).Append(r0);
            last_particle=r0;
            i+=1;}
    else{
        perturb=(points(i-1)-points(i)).Orthogonal_Vector();
        r0=particles.array_collection->Add_Element();
        particles.X(r0)=points(i);
        edges.mesh.elements.Append(Vi2(last_particle,r0));

        triangles.mesh.elements.Append(Vi3(r2,r1,r0));
        if(previous(i-1).m==2){
            r2=previous(i-1)(1);
            r1=previous(i-1)(2);
            bending_edges.mesh.elements.Append(Vi2(r2,r0));
            if(i-2>=0){
                if(previous(i-2).m==2){
                    r3=previous(i-2)(1);
                    torsion_edges.mesh.elements.Append(Vi2(r3,r0));
                    volume.mesh.elements.Append(Vi4(r3,r2,r1,r0));}
                else if(previous(i-2).m==3){
                    r3=previous(i-2)(1);
                    torsion_edges.mesh.elements.Append(Vi2(r3,r0));
                    volume.mesh.elements.Append(Vi4(r3,r2,r1,r0));}}}
        else if(previous(i-1).m==3){
            r2=previous(i-1)(1);
            e1=previous(i-1)(2);
            r1=previous(i-1)(3);
            bending_edges.mesh.elements.Append(Vi2(r2,r0));
            torsion_edges.mesh.elements.Append(Vi2(e1,r0));
            volume.mesh.elements.Append(Vi4(r2,e1,r1,r0));}
        
        
        int previous_index=previous.Append(ARRAY<int>());
        previous(previous_index).Append(last_particle);
        previous(previous_index).Append(r0);
        last_particle=r0;
        i+=1;}}

    // Fix tetrahedra orientation
    i=1;
    for(int t=0;t<volume.mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume.mesh.elements(t);
        if(TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) exchange(nodes[3],nodes[4]);
        i++;}
    for(int t=0;t<volume.mesh.elements.m;t++){
        VECTOR<int,4>& nodes=volume.mesh.elements(t);
        if(TETRAHEDRON<T>(particles.X.Subset(nodes)).Signed_Volume()<0) PHYSBAM_FATAL_ERROR();}

    //correct number nodes
    edges.Update_Number_Nodes();
    extra_edges.Update_Number_Nodes();
    bending_edges.Update_Number_Nodes();
    torsion_edges.Update_Number_Nodes();
    triangles.Update_Number_Nodes();
    volume.Update_Number_Nodes();
    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(volume,density,true);
    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);soft_bindings.Set_Mass_From_Effective_Mass();

    //forces
    T restlength_clamp=(T)1e-4;
    LINEAR_SPRINGS<TV>* edge_springs=Create_Edge_Springs(edges,hair_stiffness,hair_damping,false,(T).1,true,(T)0,true);
    edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(edge_springs);
    LINEAR_SPRINGS<TV>* extra_edge_springs=Create_Edge_Springs(extra_edges,hair_stiffness,hair_damping,false,(T).1,true,(T)0,true);
    extra_edge_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(extra_edge_springs);
    LINEAR_SPRINGS<TV>* bending_springs=Create_Edge_Springs(bending_edges,hair_stiffness,hair_damping,false,(T).1,true,(T)0,true);
    bending_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(bending_springs);
    LINEAR_SPRINGS<TV>* torsion_springs=Create_Edge_Springs(torsion_edges,hair_stiffness,hair_damping,false,(T).1,true,(T)0,true);
    torsion_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(torsion_springs);
    LINEAR_TET_SPRINGS<T> *tet_springs=Create_Tet_Springs(volume,hair_stiffness,hair_damping,false,(T).1,false,(T).1,true,(T)0,true);
    tet_springs->Clamp_Restlength(restlength_clamp);
    solid_body_collection.Add_Force(tet_springs);
    
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&edges);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&extra_edges);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&bending_edges);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&torsion_edges);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&triangles);
    solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(&volume);
        
    constrained_point.Append(1);

    if(test_number==6){
        animated_particle=100;
        animation_curve.Add_Control_Point(0,TV(0,2,0));//particles.X(100)(1),particles.X(100)(2),particles.X(100)(3)));
        animation_curve.Add_Control_Point(2,TV(0,2,0));
        animation_curve.Add_Control_Point((T)3.3,TV(0,-10,0));
        animation_curve.Add_Control_Point((T)3.6,TV(0,10,0));
        animation_curve.Add_Control_Point((T)3.8,TV(0,-10,0));
        animation_curve.Add_Control_Point(4,TV(0,2,0));
        animation_curve.Add_Control_Point(6,TV(0,2,0));
        animation_curve.Add_Control_Point(7,TV(0,2,12));
        animation_curve.Add_Control_Point((T)7.5,TV(0,2,-12));
        animation_curve.Add_Control_Point((T)7.7,TV(0,2,12));
        animation_curve.Add_Control_Point((T)7.9,TV(0,2,-12));
        animation_curve.Add_Control_Point(8,TV(0,2,0));
        animation_curve.Add_Control_Point(10,TV(0,2,0));}
}
//#####################################################################
// Function Cloth_Stretch
//#####################################################################
void Cloth_Stretch()
{
    last_frame=3000;
    test_9_skipped_particles_number=1;
    cloth_triangles=2*number_side_panels*(int)(number_side_panels*aspect_ratio);
    tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
}
//#####################################################################
// Function Write_Saved_Interpolation_Curve_To_File
//#####################################################################
void Write_Saved_Interpolation_Curve_To_File(const std::string& filename,ARRAY<INTERPOLATION_CURVE<T,TV> >& curve)
{
    std::ostream* output_raw=FILE_UTILITIES::Safe_Open_Output(filename);
    TYPED_OSTREAM output(*output_raw,stream_type);
    Write_Binary(output,curve.m);
    for(int i=0;i<curve.m;i++){
        Write_Binary(output,curve(i).control_points.m); 
        for(int j=0;j<curve(i).control_points.m;j++){
            Write_Binary(output,curve(i).control_points(j).t); 
            Write_Binary(output,curve(i).control_points(j).value(1));
            Write_Binary(output,curve(i).control_points(j).value(2));
            Write_Binary(output,curve(i).control_points(j).value(3));}}
    delete output_raw;
}
//#####################################################################
// Function Read_Saved_Interpolation_Curve_To_File
//#####################################################################
void Read_Saved_Interpolation_Curve_From_File(const std::string& filename,ARRAY<INTERPOLATION_CURVE<T,TV> >& curve)
{
    std::istream* input_raw=FILE_UTILITIES::Safe_Open_Input(filename);
    TYPED_ISTREAM input(*input_raw,stream_type);
    int m,n;
    Read_Binary(input,m);
    curve.Resize(m);
    for(int i=0;i<m;i++){
        Read_Binary(input,n);
        for(int j=0;j<n;j++){
            T t;TV value;
            Read_Binary(input,t);
            Read_Binary(input,value(1));
            Read_Binary(input,value(2));
            Read_Binary(input,value(3));
            curve(i).Add_Control_Point(t,value);}}
    delete input_raw;
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    T print_diagnostics=false;
    if((test_number==15 || test_number==16) && print_diagnostics){ // print stats
        DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
        T time_current=Time_At_Frame(frame);
        TV analytic_velocity=TV(0,-(T)9.8,0)*time_current,analytic_displacement=TV(0,-(T)4.9,0)*sqr(time_current);
        TV position_error_mixed,position_error_implicit,position_error_finescale;
        TV velocity_error_mixed,velocity_error_implicit,velocity_error_finescale;
        if(mixed_particle_index!=-1){
            velocity_error_mixed=deformable_body_collection.particles.V(mixed_particle_index)-analytic_velocity;
            position_error_mixed=deformable_body_collection.particles.X(mixed_particle_index)-initial_X_mixed_particle-analytic_displacement;}
        else LOG::cout<<"no mixed particles"<<std::endl;

        if(implicit_particle_index!=-1){
            velocity_error_implicit=deformable_body_collection.particles.V(implicit_particle_index)-analytic_velocity;
            position_error_implicit=deformable_body_collection.particles.X(implicit_particle_index)-initial_X_implicit_particle-analytic_displacement;}
        else LOG::cout<<"no implicit particles"<<std::endl;

        if(finescale_particle_index!=-1){
            velocity_error_finescale=deformable_body_collection.particles.V(finescale_particle_index)-analytic_velocity;
            position_error_finescale=deformable_body_collection.particles.X(finescale_particle_index)-initial_X_finescale_particle-analytic_displacement;}
        else LOG::cout<<"no explicit particles"<<std::endl;

        LOG::cout<<"Position errors: "<<position_error_finescale<<"   /   "<<position_error_mixed<<"   /   "<<position_error_implicit<<std::endl;
        LOG::cout<<"Velocity errors: "<<velocity_error_finescale<<"   /   "<<velocity_error_mixed<<"   /   "<<velocity_error_implicit<<std::endl;}

    if(test_number==19){
        T time=(T)frame/frame_rate,h=5-(T)4.9*sqr(time),v=-(T)9.8*time;
        LOG::cout<<"Position errors: "<<(particles.X(1)-TV(1,h,0))<<"   /   "<<(particles.X(2)-TV(2,h,0))<<"   /   "<<(particles.X(3)-TV(3,h,0))<<std::endl;
        LOG::cout<<"Velocity errors: "<<(particles.V(1)-TV(0,v,0))<<"   /   "<<(particles.V(2)-TV(0,v,0))<<"   /   "<<(particles.V(3)-TV(0,v,0))<<std::endl;}
    if(test_number==20){
        T time=(T)frame/frame_rate,e=1-exp(-2*time),x=2*e,v=4*exp(-2*time);
        if(parameter==2){
            T g=-(T)9.8;
            x+=-g/(T)28*e+g/(T)14*time;
            v+=g/(T)14*e;}
        LOG::cout<<"Position errors: "<<(particles.X(1)-TV(1,x,0))<<"   /   "<<(particles.X(2)-TV(2,x,0))<<"   /   "<<(particles.X(3)-TV(3,x,0))<<std::endl;
        LOG::cout<<"Velocity errors: "<<(particles.V(1)-TV(0,v,0))<<"   /   "<<(particles.V(2)-TV(0,v,0))<<"   /   "<<(particles.V(3)-TV(0,v,0))<<std::endl;}
    if(test_number==21){
        T t=(T)frame/frame_rate;
        T x0=initial_displacement_analytic,k=stiffness_analytic,m=mass_analytic,l0=restlength_analytic,w=sqrt(k/(m*l0));
        T dx=x0*cos(w*t),v=-x0*w*sin(w*t);
        T h=0,vh=0,z=0,vz=0;
        if(parameter==2){
            h=-(T)4.9*sqr(t);
            vh=-(T)9.8*t;}
        if(use_orthogonal_velocity){
            z=orthogonal_velocity*t;
            vz=orthogonal_velocity;}
        if(damping_analytic){
            T n=sqrt(k*m*(T)l0);
            T c=2*damping_analytic*n;
            QUADRATIC<T> quadratic(m*l0,c,k);
            quadratic.Compute_Roots();
            if(quadratic.roots==0){
                T p=-c/(2*m*l0),e=exp(p*t),w=sqrt(4*m*k*l0-sqr(c))/(2*m*l0),ce=e*cos(w*t),se=e*sin(w*t);
                T b=-p/w;
                dx=x0*(ce+b*se);
                v=-x0*(w-p*b)*se;}
            else if(quadratic.roots==1){
                T e=exp(quadratic.root1*t);
                dx=x0*(1-quadratic.root1*t)*e;
                v=-x0*sqr(quadratic.root1)*t*e;}
            else if(damping_analytic==2){
                T a=-x0*quadratic.root2/(quadratic.root1-quadratic.root2),b=x0-a,e1=exp(quadratic.root1*t),e2=exp(quadratic.root2*t);
                dx=a*e1+b*e2;
                v=a*quadratic.root1*e1+b*quadratic.root2*e2;}
        }
        LOG::cout<<"Position errors: "<<(particles.X(1)-TV(-dx-restlength_analytic,h,z))<<"   /   "<<(particles.X(2)-TV(0,h,z))<<"   /   "<<(particles.X(3)-TV(dx+restlength_analytic,h,z))<<std::endl;
        LOG::cout<<"Velocity errors: "<<(particles.V(1)-TV(-v,vh,vz))<<"   /   "<<(particles.V(2)-TV(0,vh,vz))<<"   /   "<<(particles.V(3)-TV(v,vh,vz))<<std::endl;}
    if(test_number==22){
        PHYSBAM_ASSERT(frame==1);
        T fx=-stiffness_analytic/restlength_analytic;
        T fv=-damping_analytic*sqrt(2*mass_analytic*stiffness_analytic*restlength_analytic)/restlength_analytic;
        VECTOR<T,2> analytic_step;
        if(parameter==1) analytic_step=Analytic_Solution_BE(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==2) analytic_step=Analytic_Solution_TR(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==3) analytic_step=Analytic_Solution_Full_BE(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==4) analytic_step=Analytic_Solution_Full_TR(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==5) analytic_step=Analytic_Solution_Asyn_Stability_Implicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==6) analytic_step=Analytic_Solution_Asyn_Stability_Explicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==7) analytic_step=Analytic_Solution_Asyn_Stability_Mixed(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx,0,fv,fx);
        else if(parameter==8) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Implicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==9) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Explicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==10) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Mixed(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx,0,fv,fx);
        else if(parameter==11) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Implicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==12) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Explicit(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx);
        else if(parameter==13) analytic_step=Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Mixed(initial_displacement_analytic,initial_velocity_analytic,mass_analytic,(T)1/frame_rate,0,fv,fx,0,fv,fx);

        LOG::cout<<"actual:  x "<<particles.X(1).x<<"     v "<<particles.V(1).x<<std::endl;
        LOG::cout<<"analytic:  x "<<analytic_step.x<<"     v "<<analytic_step.y<<std::endl;
        LOG::cout<<"Errors:  x "<<(particles.X(1).x-analytic_step.x)<<"     v "<<(particles.V(1).x-analytic_step.y)<<std::endl;
    }
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    LOG::cout<<"Preprocess Frame "<<frame<<std::endl;
    if(test_number==8 && frame==1){
        RANDOM_NUMBERS<T> random;random.Set_Seed(1823);
        T perturbation_size=side_length/number_side_panels*4;
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        for(int p=0;p<particles.array_collection->Size();p++) particles.X(p).y+=random.Get_Uniform_Number((T)0,perturbation_size);}
    if((test_number==10 || test_number==12) && frame<=sim_length+1)
        for(int i=0;i<num_controlled_particles;i++){
            saved_V(i).Add_Control_Point(Time_At_Frame(frame-1),deformable_body_collection.particles.V(i));
            saved_X(i).Add_Control_Point(Time_At_Frame(frame-1),deformable_body_collection.particles.X(i));}
    if(((test_number==10 || test_number==12) && frame==sim_length+1) || ((test_number==11 || test_number==13) && frame==sim_length+1)){
        if(test_number==10 || test_number==12){
            Write_Saved_Interpolation_Curve_To_File("saved_X",saved_X);
            Write_Saved_Interpolation_Curve_To_File("saved_V",saved_V);}
        if(test_number==11 || test_number==13){
            Read_Saved_Interpolation_Curve_From_File("saved_X",saved_X);
            Read_Saved_Interpolation_Curve_From_File("saved_V",saved_V);}

        for(HASHTABLE_ITERATOR<int,int> it(particle_map);it.Valid();it.Next()){
            deformable_body_collection.particles.X(it.Data())=saved_X(it.Key()).Value(0);
            deformable_body_collection.particles.V(it.Data())=saved_V(it.Key()).Value(0);
            if(test_number!=12) solid_body_collection.deformable_body_collection.soft_bindings.Add_Binding(VECTOR<int,2>(it.Data(),it.Key()),false);}
        for(int i=0;i<num_controlled_particles;i++){
            deformable_body_collection.particles.X(i)=saved_X(i).Value(0);
            deformable_body_collection.particles.V(i)=saved_V(i).Value(0);}
        Activate_Secondary_Simulation();}
    if(asynchronous_evolution) asynchronous_evolution->Preprocess_Frame(frame,driver->example.Time_At_Frame(frame),driver->Time());
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time)
{
    if(asynchronous_evolution) asynchronous_evolution->Preprocess_Substep(dt,time);
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(asynchronous_evolution) asynchronous_evolution->Postprocess_Substep(dt,time);
}
//#####################################################################
// Function Activate_Secondary_Simulation
//#####################################################################
void Activate_Secondary_Simulation()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    COLLISION_GEOMETRY_COLLECTION<TV>& collision_body_list=solid_body_collection.collision_body_list;
    // Remove old forces and add new ones
    solid_body_collection.solids_forces.Remove_All();
    T stiffness=(T)1e4*soft_surface_multiplier;
    if(test_number==12 || test_number==13){
        //tests.Add_Gravity();
        for(int i=0;i<number_of_spheres;i++){
            solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,boundary_one_ring_segmented_curves(i)->mesh));
            solid_body_collection.Add_Force(Create_Altitude_Springs(*boundary_tetrahedralized_volumes(i),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping));
            solid_body_collection.Add_Force(Create_Edge_Springs(*boundary_one_ring_segmented_curves(i),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping));}}
    else{
        for(int i=0;i<number_of_spheres;i++){
            solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,primary_segmented_curves(i)->mesh));
            solid_body_collection.Add_Force(Create_Edge_Springs(*new_boundary_segmented_curves(i),(T)stiffness/(1+sqrt((T)2)),soft_bound_edge_damping));}

        // Add bindings
        T binding_stiffness=(T)1e4*soft_bindings_multiplier;
        solid_body_collection.deformable_body_collection.soft_bindings.use_impulses_for_collisions.Fill(false);
        solid_body_collection.deformable_body_collection.soft_bindings.Initialize_Binding_Mesh();
        solid_body_collection.Add_Force(Create_Edge_Binding_Springs(deformable_body_collection.particles,*solid_body_collection.deformable_body_collection.soft_bindings.binding_mesh,binding_stiffness,(T).01));

        deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);
        solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();}

    // Disable collisions
    deformable_body_collection.collisions.collision_structures.Remove_All();
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Remove_All();
    for(COLLISION_GEOMETRY_ID id(1);id<=collision_body_list.Size();id++) collision_body_list.Remove_Body(id);

    solid_body_collection.Update_Simulated_Particles();

    //if(fully_implicit) for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;

    ground->X().y=-999;
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==8){
        if(time>test_8_wind_off){
            INTERPOLATION_CURVE<T,TV> wind_stop;
            wind_stop.Add_Control_Point(test_8_wind_off-(T).3,TV((T).001,1,(T)-.25));
            wind_stop.Add_Control_Point(test_8_wind_off,TV());
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(0);drag.Use_Constant_Wind(0,TV());}
        else if(time>test_8_friction_off){
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(20);drag.Use_Constant_Wind(0,TV((T).001,1,(T)-.25));}}
    if(asynchronous_evolution) asynchronous_evolution->Update_Time_Varying_Material_Properties(time);
}
int Particle_Index_From_Coordinates(const int i,const int j,const int m)
{
    return i+m*(j-1);
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE 
{
    if(asynchronous_evolution) asynchronous_evolution->Set_External_Velocities(V,velocity_time,current_position_time);
    V.Subset(constrained_point).Fill(TV());
    if(animated_particle) V(animated_particle)=animation_curve.Derivative(velocity_time);
    if(test_number==8){
        if(velocity_time<test_8_constrained_off){
            INTERPOLATION_CURVE<T,T> x;
            x.Add_Control_Point((T)0,(T)5);
            x.Add_Control_Point((T)test_8_constrained_off,(T)0);
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;

            int j=(int)((T).8*n);
            int i=(int)((T).33*m); // (.33,.8)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=x.Value(velocity_time);
            i=(int)((T).66*m);// (.66,.8)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=-x.Value(velocity_time);

            j=(int)(.6*n);
            i=(int)(.25*m);// (.25,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=(T)1.5*x.Value(velocity_time);
            i=(int)(.5*m);// (.5,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=-(T)1.5*x.Value(velocity_time);
            i=(int)(.65*m);// (.65,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=(T)1.5*x.Value(velocity_time);

            j=(int)(.5*n);i=(int)(.33*m);// (.33,.5)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=x.Value(velocity_time);}}
    else if(test_number==9){
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            for(int j=1;j<=n;j+=test_9_skipped_particles_number) V(Particle_Index_From_Coordinates(1,j,m))=TV();
            for(int j=1;j<=n;j+=test_9_skipped_particles_number) V(Particle_Index_From_Coordinates(m,j,m))=TV((T).1,0,0);}
    else if((test_number==10 || test_number==11 || test_number==12 || test_number==13) && velocity_time>sim_switch_time)
            for(int i=0;i<num_controlled_particles;i++) V(i)=saved_V(i).Value(velocity_time-subtract_time);
    if(test_number==22) V(2)=TV();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(asynchronous_evolution) asynchronous_evolution->Set_External_Positions(X,time);
    if(animated_particle) X(animated_particle)=animation_curve.Value(time);
    if((test_number==10 || test_number==11 || test_number==12 || test_number==13) && time>sim_switch_time){
        for(int i=0;i<num_controlled_particles;i++) X(i)=saved_X(i).Value(time-subtract_time);}
    if(asynchronous_evolution) asynchronous_evolution->Set_External_Positions(X,time);
    if(test_number==22) X(2)=TV(-restlength_analytic,0,0);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(asynchronous_evolution) asynchronous_evolution->Zero_Out_Enslaved_Velocity_Nodes(V,velocity_time,current_position_time);
    V.Subset(constrained_point).Fill(TV());
    if(animated_particle) V(animated_particle)=TV();
    if(test_number==8){
        if(velocity_time<test_8_constrained_off){
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;

            int j=(int)((T).8*n);
            int i=(int)((T).33*m); // (.33,.8)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=0;
            i=(int)((T).66*m);// (.66,.8)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=-0;

            j=(int)(.6*n);
            i=(int)(.25*m);// (.25,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;
            i=(int)(.5*m);// (.5,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;
            i=(int)(.65*m);// (.65,.6)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;

            j=(int)(.5*n);i=(int)(.33*m);// (.33,.5)
            V(Particle_Index_From_Coordinates(i,j,m)).z=0;V(Particle_Index_From_Coordinates(i,j,m)).x=0;}}
    else if(test_number==9){
        int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        for(int j=1;j<=n;j+=test_9_skipped_particles_number) V(Particle_Index_From_Coordinates(1,j,m))=TV();
        for(int j=1;j<=n;j+=test_9_skipped_particles_number) V(Particle_Index_From_Coordinates(m,j,m))=TV();}
    else if((test_number==10 || test_number==11 || test_number==12 || test_number==13) && velocity_time>sim_switch_time){
        for(int i=0;i<num_controlled_particles;i++) V(i)=TV();}
    if(test_number==22) V(2)=TV();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==int(2)){
        if(time<2) frame.t=TV(0,0,(T).5);
        else if(time<(T)3.5) frame.t=TV((time-2),(T).5*(time-2),(T).5);
        else if(time<4) frame.t=TV((T)1.5,(T)((T).75-(T)1.5*(time-(T)3.5)),(T).5);
        else frame.t=TV((T)(1.5-1.5*(time-4)),0,(T).5);}
    else if((test_number==7) && id==int(2)){
        frame.t=TV(0,(T).30,0);
        if(time>(T).75) frame.r=ROTATION<TV>::From_Rotation_Vector(time*TV(0,(T)-pi/4,0));}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==3 && id==int(2)){
        if(time<2) twist.linear=TV();
        else if(time<(T)3.5) twist.linear=TV(1,(T).5,0);
        else if(time<4) twist.linear=TV(0,(T)-1.5,0);
        else twist.linear=TV(-1.5,0,0);}
    else if((test_number==7) && id==int(2)){if(time>(T).75){twist.angular=TV(0,(T)-pi/4,0);}}
    else return false;
    return true;
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE
{
    if(asynchronous_evolution) asynchronous_evolution->Limit_Solids_Dt(dt,time);
    if(substeps_per_frame) dt=Inverse(substeps_per_frame*frame_rate);
}
//#####################################################################
// Function Soft_Bound_Sphere_Surface
//#####################################################################
void Soft_Bound_Sphere_Surface()
{
    last_frame=600;
    if(test_number==11 || test_number==13) last_frame=300;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    for(int i=0;i<number_of_spheres;i++){
        primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)i,(T)(2*i+1),0))),true,true,1000));}
    ground=&tests.Add_Ground();

    if(test_number==10 || test_number==12){sim_length=last_frame/2;last_frame=sim_length*2;subtract_time=Time_At_Frame(sim_length);sim_switch_time=Time_At_Frame(sim_length);}
    else if(test_number==11 || test_number==13){sim_length=1;subtract_time=Time_At_Frame(sim_length);sim_switch_time=Time_At_Frame(sim_length);}

    num_controlled_particles=deformable_body_collection.particles.array_collection->Size();
    saved_V.Resize(num_controlled_particles);
    saved_X.Resize(num_controlled_particles);

    Add_Secondary_Simulation_Structures();
}
//#####################################################################
// Function Add_Secondary_Simulation_Structures
//#####################################################################
void Add_Secondary_Simulation_Structures()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    // TODO: What about frame zero.
    LOG::cout<<"MAKING COPY OF MESH"<<std::endl;
    for(int i=0;i<number_of_spheres;i++){
        TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(i);
        TRIANGULATED_SURFACE<T>* new_surface=Add_Copy_of_Boundary_Object(volume,particle_map);
        new_boundary_surfaces.Append(new_surface);}

    if(test_number==10 || test_number==11){
        // Add new_surface segment mesh for edge springs
        for(int i=0;i<number_of_spheres;i++){
            TRIANGULATED_SURFACE<T>* new_surface=new_boundary_surfaces(i);
            new_surface->mesh.Initialize_Segment_Mesh();
            SEGMENTED_CURVE<TV>* new_boundary_segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
            new_boundary_segmented_curve->mesh.elements=new_surface->mesh.segment_mesh->elements;
            new_boundary_segmented_curve->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(new_boundary_segmented_curve);
            new_boundary_segmented_curves.Append(new_boundary_segmented_curve);}}
    
    if(test_number==12 || test_number==13){
        // Add boundary one ring segment mesh for edge springs
        for(int i=0;i<number_of_spheres;i++){
            TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(i);
            volume->mesh.Initialize_Segment_Mesh();
            SEGMENTED_CURVE<TV> *boundary_one_ring_segmented_curve=SEGMENTED_CURVE<TV>::Create(solid_body_collection.deformable_body_collection.particles);
            Add_Mapped_Elements<2>(volume->mesh.segment_mesh->elements,boundary_one_ring_segmented_curve->mesh.elements,particle_map);
            boundary_one_ring_segmented_curve->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(boundary_one_ring_segmented_curve);
            boundary_one_ring_segmented_curves.Append(boundary_one_ring_segmented_curve);
        
            // Add boundary tet mesh for altitude springs
            TETRAHEDRALIZED_VOLUME<T>* boundary_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
            Add_Mapped_Elements<4>(volume->mesh.elements,boundary_tetrahedralized_volume->mesh.elements,particle_map);
            boundary_tetrahedralized_volume->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(boundary_tetrahedralized_volume);
            boundary_tetrahedralized_volumes.Append(boundary_tetrahedralized_volume);}}

    // add structures and rigid bodies to collisions
}
//#####################################################################
// Function Add_Copy_of_Boundary_Object
//#####################################################################
TRIANGULATED_SURFACE<T>* Add_Copy_of_Boundary_Object(TETRAHEDRALIZED_VOLUME<T>* volume,HASHTABLE<int,int>& particle_map)
{
    volume->Update_Number_Nodes();
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    TRIANGULATED_SURFACE<T>* compact_surface_copy=volume->Get_Boundary_Object().Create_Compact_Copy();
    TRIANGULATED_SURFACE<T>* new_surface=&tests.Copy_And_Add_Structure(*compact_surface_copy,0);
    for(int i=0;i<volume->triangulated_surface->mesh.elements.m;i++) for(int j=0;j<3;j++){
        int old_index=volume->triangulated_surface->mesh.elements(i)(j);
        int new_index=new_surface->mesh.elements(i)(j);
        bool new_entry=particle_map.Set(old_index,new_index);
        if(new_entry){
            deformable_body_collection.particles.mass(new_index)=deformable_body_collection.particles.mass(old_index);
            deformable_body_collection.particles.X(new_index).x+=10;}}
    return new_surface;
}
//#####################################################################
// Function Add_Mapped_Elements
//#####################################################################
template<int d>
void Add_Mapped_Elements(const ARRAY<VECTOR<int,d> >& elements,ARRAY<VECTOR<int,d> >& new_elements,const HASHTABLE<int,int>& map)
{
    VECTOR<int,d> new_element;
    for(int i=0;i<elements.m;i++){
        bool element_mapped=false;
        new_element=elements(i);
        for(int j=0;j<d;j++){
            const int *mapped_index=map.Get_Pointer(elements(i)(j));
            if(mapped_index){new_element[j]=*mapped_index;element_mapped=true;}}
        if(element_mapped) new_elements.Append(new_element);}
}
//#####################################################################
// Function Expand_One_Ring_Particles
//#####################################################################
void Expand_One_Ring_Particles(const ARRAY<VECTOR<int,2> >& elements,const HASHTABLE<int,int>& all_old_particles,const HASHTABLE<int,int>& frontier_particles,HASHTABLE<int,int>& new_particles)
{
    VECTOR<int,2> edge;
    for(int i=0;i<elements.m;i++){
        edge=elements(i);
        const int *front_index_1=frontier_particles.Get_Pointer(edge(1));
        const int *front_index_2=frontier_particles.Get_Pointer(edge(2));
        if(front_index_1 && !front_index_2){
            const int *mapped_index_2=all_old_particles.Get_Pointer(edge(2));
            if(!mapped_index_2)
                new_particles.Set(edge(2),edge(2));}
        else if(!front_index_1 && front_index_2){
            const int *mapped_index_1=all_old_particles.Get_Pointer(edge(1));
            if(!mapped_index_1)
                new_particles.Set(edge(1),edge(1));}}
}
//#####################################################################
// Function Expand_N_Rings
//#####################################################################
void Expand_N_Rings(const int n,const ARRAY<VECTOR<int,2> >& elements,HASHTABLE<int,int>& all_old_particles,HASHTABLE<int,int>& frontier_particles)
{
    HASHTABLE<int,int> new_particles;
    for(int i=0;i<n;i++){
        new_particles.Clean_Memory();
        Expand_One_Ring_Particles(elements,all_old_particles,frontier_particles,new_particles);
        for(HASHTABLE_ITERATOR<int,int> it(new_particles);it.Valid();it.Next()) all_old_particles.Set(it.Key(),it.Data());
        frontier_particles=new_particles;}
}
//#####################################################################
// Function Asynchronous_Sphere
//#####################################################################
void Asynchronous_Sphere()
{
    last_frame=300;
    // parameter=2: constrain top point to be fixed

    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    fully_implicit=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=(number_of_spheres>=2);
    ground=&tests.Add_Ground();
    //ground->X().y=-999;

    if(test_number==15 || test_number==16) asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);

    if(test_number==16) number_of_spheres=3;
    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        particle_map.Clean_Memory();
        if(test_number==16){
            switch(sphere_index){
                case 1: number_of_rings=0;break;
                case 2: number_of_rings=1;break;
                case 3: number_of_rings=-1;break;}
            primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)(sphere_index-2)*4+height_offset,(T)4,0))),true,true,1000));}
        else{
            primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)sphere_index,(T)(2*sphere_index+1)+height_offset,0))),true,true,1000));}
        TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(sphere_index);
        volume->mesh.Initialize_Segment_Mesh();
        ARRAY_VIEW<int> flattened_volume=volume->mesh.elements.Flattened();
        HASHTABLE<int> volume_particle_map;
        for(int i=0;i<flattened_volume.Size();i++) volume_particle_map.Set(flattened_volume(i));
        ARRAY<int> all_particles;
        volume_particle_map.Get_Keys(all_particles);
        
        TV center=TV();
        switch(type_of_asynchronous_configuration){
            case 1: //boundary
                if(number_of_rings>0){
                    TRIANGULATED_SURFACE<T>& surface=volume->Get_Boundary_Object();
                    ARRAY_VIEW<int> flattened(surface.mesh.elements.Flattened());
                    PARTICLES<TV>& particles=deformable_body_collection.particles;
                    T top=-FLT_MAX,bottom=FLT_MAX;
                    for(int i=0;i<flattened.Size();i++){int p=flattened(i);
                        if(particles.X(p)(2)>top) top=particles.X(p)(2);
                        if(particles.X(p)(2)<bottom) bottom=particles.X(p)(2);}
                    for(int i=0;i<flattened.Size();i++){int p=flattened(i);
                        if(particles.X(p)(2)<=(top-bottom)*coverage_percent+bottom)
                            particle_map.Set(p,p);}
                    HASHTABLE<int,int> frontier_particle_map=particle_map;
                    Expand_N_Rings(number_of_rings-1,volume->mesh.segment_mesh->elements,particle_map,frontier_particle_map);}
                else if(number_of_rings<0){
                    ARRAY_VIEW<int> flattened(volume->mesh.elements.Flattened());
                    for(int i=0;i<flattened.Size();i++) particle_map.Set(flattened(i),flattened(i));}
                else if(number_of_rings==0){
                    particle_map.Clean_Memory();}
                break;
            case 2: //left
                for(int i=0;i<all_particles.Size();i++) center+=deformable_body_collection.particles.X(all_particles(i));
                center/=(T)all_particles.m;
                for(int i=0;i<all_particles.m;i++) if(deformable_body_collection.particles.X(all_particles(i))(1)<center(1))
                    particle_map.Set(i,i);
                break;
            case 3: //top
                for(int i=0;i<all_particles.m;i++) center+=deformable_body_collection.particles.X(all_particles(i));
                center/=(T)all_particles.m;
                for(int i=0;i<all_particles.m;i++) if(deformable_body_collection.particles.X(all_particles(i))(2)<center(2))
                    particle_map.Set(i,i);
                break;
            default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized config type %d",type_of_asynchronous_configuration));}

        if(invert_asynchronous_implicit){
            HASHTABLE<int,int> copy_particle_map=particle_map;
            particle_map.Clean_Memory();
            for(int i=0;i<all_particles.m;i++){int p=all_particles(i);
                if(!copy_particle_map.Contains(p)) particle_map.Set(p,p);}}
 
        SEGMENTED_CURVE<TV> *boundary_one_ring_segmented_curve=SEGMENTED_CURVE<TV>::Create(solid_body_collection.deformable_body_collection.particles);
        Add_Mapped_Elements<2>(volume->mesh.segment_mesh->elements,boundary_one_ring_segmented_curve->mesh.elements,particle_map);
        boundary_one_ring_segmented_curve->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(boundary_one_ring_segmented_curve);
        boundary_one_ring_segmented_curves.Append(boundary_one_ring_segmented_curve);

        // Add boundary tet mesh for altitude springs
        TETRAHEDRALIZED_VOLUME<T>* boundary_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
        Add_Mapped_Elements<4>(volume->mesh.elements,boundary_tetrahedralized_volume->mesh.elements,particle_map);
        boundary_tetrahedralized_volume->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(boundary_tetrahedralized_volume);
        boundary_tetrahedralized_volumes.Append(boundary_tetrahedralized_volume);
    
        if(test_number==15 || test_number==16){
            SEGMENTED_CURVE<TV>* complementary_boundary_segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
            HASHTABLE<VECTOR<int,2> > edge_mesh_element_map;
            for(int i=0;i<boundary_one_ring_segmented_curve->mesh.elements.m;i++) edge_mesh_element_map.Set(boundary_one_ring_segmented_curve->mesh.elements(i));
            for(int i=0;i<volume->mesh.segment_mesh->elements.m;i++) if(!edge_mesh_element_map.Contains(volume->mesh.segment_mesh->elements(i)))
                complementary_boundary_segmented_curve->mesh.elements.Append(volume->mesh.segment_mesh->elements(i));
            complementary_boundary_segmented_curve->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(complementary_boundary_segmented_curve);
            complementary_boundary_segmented_curves.Append(complementary_boundary_segmented_curve);

            // Add complementary boundary tet mesh for altitude, edge springs
            TETRAHEDRALIZED_VOLUME<T>* complementary_boundary_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
            HASHTABLE<VECTOR<int,4> > volume_mesh_element_map;
            for(int i=0;i<boundary_tetrahedralized_volume->mesh.elements.m;i++) volume_mesh_element_map.Set(boundary_tetrahedralized_volume->mesh.elements(i));
            for(int i=0;i<volume->mesh.elements.m;i++) if(!volume_mesh_element_map.Contains(volume->mesh.elements(i)))
                complementary_boundary_tetrahedralized_volume->mesh.elements.Append(volume->mesh.elements(i));
            complementary_boundary_tetrahedralized_volume->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(complementary_boundary_tetrahedralized_volume);
            complementary_boundary_tetrahedralized_volumes.Append(complementary_boundary_tetrahedralized_volume);}}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add forces
    T stiffness=(T)1e4*soft_surface_multiplier;
    DEFORMABLES_FORCES<TV>* force;
    ARRAY<int> affected_particle_indices,affected_rigid_body_particle_indices;
    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        
        if(!boundary_one_ring_segmented_curves(sphere_index)->mesh.elements.m) continue;

        // Get affected_particle_indices
        boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);

        force=new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,boundary_one_ring_segmented_curves(sphere_index)->mesh);
        if(test_number==14) solid_body_collection.Add_Force(force);
        if(asynchronous_evolution) asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Altitude_Springs(*boundary_tetrahedralized_volumes(sphere_index),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        if(test_number==14) solid_body_collection.Add_Force(force);
        if(asynchronous_evolution) asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Edge_Springs(*boundary_one_ring_segmented_curves(sphere_index),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        if(test_number==14) solid_body_collection.Add_Force(force);
        if(asynchronous_evolution) asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);}

    if(test_number==15 || test_number==16){
        HASHTABLE<int> temp_particle_map;
        for(int i=0;i<number_of_spheres;i++){
            ARRAY_VIEW<int> boundary_particle_index=boundary_tetrahedralized_volumes(i)->mesh.elements.Flattened();
            for(int j=0;j<boundary_particle_index.Size();j++) temp_particle_map.Set(boundary_particle_index(j));}
        for(int i=0;i<deformable_body_collection.particles.array_collection->Size();i++) if(!temp_particle_map.Contains(i))
            gravity_particles.Append(i);}

    if(parameter==2){
        animated_particle=4270;
        animation_curve.Add_Control_Point(0,deformable_body_collection.particles.X(animated_particle));
        animation_curve.Add_Control_Point(1,deformable_body_collection.particles.X(animated_particle));}
}
//#####################################################################
// Function Asynchronous_Armadillo
//#####################################################################
void Asynchronous_Armadillo()
{
    last_frame=3000;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.cfl=(T)2;
    ground=&tests.Add_Ground();
    T armadillo_scale=(T).085;

    TETRAHEDRALIZED_VOLUME<T>* volume=&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_20K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)10,0)),TWIST<TV>(TV(),TV(0,0,1))),true,false,1000,armadillo_scale);
    solid_body_collection.Add_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,volume->mesh,98));
}
//#####################################################################
// Function Add_Layered_Elements
//#####################################################################
template<int d>
void Add_Layered_Elements(const ARRAY<VECTOR<int,d> >& elements,ARRAY<VECTOR<int,d> >& new_elements,ARRAY<VECTOR<int,d> >& new_complementary_elements,const int particle_index_offset,const GRID<TV>& grid,const int total_layer_number,const int asynchronous_layer_interval,const int non_layer_interval)
{
    for(int i=0;i<elements.m;i++){
        int min_layer=total_layer_number;
        for(int j=0;j<d;j++){
            int layer=(elements(i)(j)-particle_index_offset)%(grid.counts.x*grid.counts.y)/grid.counts.x;
            if(layer<min_layer) min_layer=layer;}
        min_layer%=asynchronous_layer_interval+non_layer_interval;
        if(min_layer<asynchronous_layer_interval) new_elements.Append(elements(i));
        else new_complementary_elements.Append(elements(i));}
}
//#####################################################################
// Function Asynchronous_Layered_Box
//#####################################################################
void Asynchronous_Layered_Box()
{
    Asynchronous_Layered_Box();
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    fully_implicit=false;
    ground=&tests.Add_Ground();

    asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);

    int number_of_boxes=3,total_layer_number=10,asynchronous_layer_interval,non_layer_interval;
    int particle_index_offset;
    for(int box_index=0;box_index<number_of_boxes;box_index++){
        switch(box_index){
            case 1: asynchronous_layer_interval=0;non_layer_interval=10;break;
            case 2: asynchronous_layer_interval=1;non_layer_interval=1;break;
            case 3: asynchronous_layer_interval=10;non_layer_interval=0;break;
            default: asynchronous_layer_interval=01;non_layer_interval=-1;break;}
        GRID<TV> grid(VECTOR<int,3>(10,total_layer_number,10),BOX<TV>(TV((T)(-1+3*(box_index-3)),(T).1,-(T)1),TV((T)(1+3*(box_index-3)),(T)1.1,(T)1)));
        particle_index_offset=deformable_body_collection.particles.array_collection->Size()+1;
        primary_tetrahedralized_volumes.Append(&tests.Create_Mattress(grid));
        TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(box_index);
        volume->mesh.Initialize_Segment_Mesh();

        SEGMENTED_CURVE<TV> *asynchronous_segmented_curve=SEGMENTED_CURVE<TV>::Create(solid_body_collection.deformable_body_collection.particles);
        SEGMENTED_CURVE<TV> *complementary_segmented_curve=SEGMENTED_CURVE<TV>::Create(solid_body_collection.deformable_body_collection.particles);
        Add_Layered_Elements<2>(volume->mesh.segment_mesh->elements,asynchronous_segmented_curve->mesh.elements,complementary_segmented_curve->mesh.elements,particle_index_offset,grid,total_layer_number,asynchronous_layer_interval,non_layer_interval);
        asynchronous_segmented_curve->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(asynchronous_segmented_curve);
        boundary_one_ring_segmented_curves.Append(asynchronous_segmented_curve);
        complementary_segmented_curve->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(complementary_segmented_curve);
        complementary_boundary_segmented_curves.Append(complementary_segmented_curve);

        // Add boundary tet mesh for altitude springs
        TETRAHEDRALIZED_VOLUME<T>* asynchronous_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
        TETRAHEDRALIZED_VOLUME<T>* complementary_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
        Add_Layered_Elements<4>(volume->mesh.elements,asynchronous_volume->mesh.elements,complementary_volume->mesh.elements,particle_index_offset,grid,total_layer_number,asynchronous_layer_interval,non_layer_interval);
        asynchronous_volume->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(asynchronous_volume);
        boundary_tetrahedralized_volumes.Append(asynchronous_volume);
        complementary_volume->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(complementary_volume);
        complementary_boundary_tetrahedralized_volumes.Append(complementary_volume);}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();
    
    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add forces
    T stiffness=(T)1e4*soft_surface_multiplier;
    DEFORMABLES_FORCES<TV>* force;
    ARRAY<int> affected_particle_indices,affected_rigid_body_particle_indices;
    for(int box_index=0;box_index<number_of_boxes;box_index++){
        
        if(!boundary_one_ring_segmented_curves(box_index)->mesh.elements.m) continue;

        // Get affected_particle_indices
        boundary_tetrahedralized_volumes(box_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);

        force=new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,boundary_one_ring_segmented_curves(box_index)->mesh);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Altitude_Springs(*boundary_tetrahedralized_volumes(box_index),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Edge_Springs(*boundary_one_ring_segmented_curves(box_index),(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);}

    stiffness=(T)1e4*stiffness_multiplier;
    for(int box_index=0;box_index<number_of_boxes;box_index++){
        
        if(!complementary_boundary_segmented_curves(box_index)->mesh.elements.m) continue;
        // Get affected_particle_indices
        complementary_boundary_tetrahedralized_volumes(box_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);

        force=Create_Altitude_Springs(*complementary_boundary_tetrahedralized_volumes(box_index),(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
        asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,false);

        force=Create_Edge_Springs(*complementary_boundary_segmented_curves(box_index),(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
        asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,true);}

    HASHTABLE<int> temp_particle_map;
    for(int i=0;i<number_of_boxes;i++){
        ARRAY_VIEW<int> boundary_particle_index=boundary_tetrahedralized_volumes(i)->mesh.elements.Flattened();
        for(int j=0;j<boundary_particle_index.Size();j++) temp_particle_map.Set(boundary_particle_index(j));}
    for(int i=0;i<deformable_body_collection.particles.array_collection->Size();i++) if(!temp_particle_map.Contains(i))
        gravity_particles.Append(i);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();
        
    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Set_Driver
//#####################################################################
void Set_Driver(SOLIDS_FLUIDS_DRIVER<TV>* driver_input)
{
    driver=dynamic_cast<SOLIDS_FLUIDS_DRIVER_UNIFORM<GRID<TV> >*>(driver_input);
}
//#####################################################################
// Function Gravity_Test
//#####################################################################
void Gravity_Test()
{
    last_frame=20;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    ARRAY<int> empty_list;
    LOG::cout<<std::setprecision(16);

    asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);

    for(int i=0;i<3;i++){
        int p=particles.array_collection->Add_Element();
        particles.X(p)=TV((T)i,5,0);
        particles.mass(p)=7;
        int_lists[i-1].Append(p);}

    SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(2,3));
    segmented_curve->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);

    asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[0]),int_lists[0],empty_list,fine_fully_implicit);
    asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[1],(T)4.9),int_lists[1],empty_list,fine_fully_implicit);
    asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[1],(T)4.9),int_lists[1],empty_list,coarse_fully_implicit,true);
    asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[2]),int_lists[2],empty_list,coarse_fully_implicit,true);
}
//#####################################################################
// Function Ether_Drag_Test
//#####################################################################
void Ether_Drag_Test()
{
    last_frame=20;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    ARRAY<int> empty_list;
    LOG::cout<<std::setprecision(16);

    asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);

    for(int i=0;i<3;i++){
        int p=particles.array_collection->Add_Element();
        particles.X(p)=TV((T)i,0,0);
        particles.V(p)=TV(0,4,0);
        particles.mass(p)=7;
        int_lists[i-1].Append(p);}

    SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(2,3));
    segmented_curve->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);

    asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_ETHER_DRAG<TV>(particles,&int_lists[0],2),int_lists[0],empty_list,fine_fully_implicit);
    asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_ETHER_DRAG<TV>(particles,&int_lists[1],1),int_lists[1],empty_list,fine_fully_implicit);
    asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_ETHER_DRAG<TV>(particles,&int_lists[1],1),int_lists[1],empty_list,coarse_fully_implicit,true);
    asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_ETHER_DRAG<TV>(particles,&int_lists[2],2),int_lists[2],empty_list,coarse_fully_implicit,true);
    
    if(parameter==2){
        asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[0]),int_lists[0],empty_list,fine_fully_implicit);
        asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[1],(T)4.9),int_lists[1],empty_list,fine_fully_implicit);
        asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[1],(T)4.9),int_lists[1],empty_list,coarse_fully_implicit,false);
        asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[2]),int_lists[2],empty_list,coarse_fully_implicit,false);}
}
//#####################################################################
// Function Spring_Test
//#####################################################################
void Spring_Test()
{
    last_frame=20;
    // parameter=1: two springs along x with initial displacement, no gravity
    // parameter=2: same as 1 with gravity
    // parameter=3: two springs along y with no initial displacement, with gravity
    
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>* segmented_curve=0;
    ARRAY<int> empty_list;
    LOG::cout<<std::setprecision(16);

    asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);
    mass_analytic=(T)7;stiffness_analytic=(T)1e3*stiffness_multiplier;damping_analytic=damping_multiplier;restlength_analytic=(T)1;
    if(parameter==3) initial_displacement_analytic=(T)0;
    else initial_displacement_analytic=(T)1;
    orthogonal_velocity=0;if(use_orthogonal_velocity) orthogonal_velocity=100;

    for(int i=0;i<3;i++){
        int p=particles.array_collection->Add_Element();
        if(parameter==3) particles.X(p)=TV(-1,(i-2)*restlength_analytic,0);
        else particles.X(p)=TV((i-2)*restlength_analytic,0,0);
        particles.mass(p)=mass_analytic;
        particles.V(p)=TV(0,0,orthogonal_velocity);}
    int_lists[0].Append(1);int_lists[0].Append(2);int_lists[1].Append(3);

    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    DEFORMABLES_FORCES<TV>* force;
    segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
    force=Create_Edge_Springs(particles,segmented_curve->mesh,stiffness_analytic,damping_analytic);
    asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,fine_fully_implicit);

    segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(2,3));
    segmented_curve->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
    force=Create_Edge_Springs(particles,segmented_curve->mesh,stiffness_analytic,damping_analytic);
    asynchronous_evolution->Add_Coarsescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,coarse_fully_implicit,true);
    //asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,fine_fully_implicit);

    particles.X(1)-=TV(initial_displacement_analytic,0,0);
    particles.X(3)+=TV(initial_displacement_analytic,0,0);

    if(parameter==2 || parameter==3){
        asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[0]),int_lists[0],empty_list,fine_fully_implicit);
        asynchronous_evolution->Add_Coarsescale_Force(new DEFORMABLE_GRAVITY<TV>(particles,&int_lists[1]),int_lists[1],empty_list,coarse_fully_implicit,false);}
}
//#####################################################################
// Function One_Spring_Test
//#####################################################################
void One_Spring_Test()
{
    last_frame=1;
    if(parameter>=5 && parameter<=10) substeps_per_frame=3;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    SEGMENTED_CURVE<TV>* segmented_curve=0;
    ARRAY<int> empty_list;
    LOG::cout<<std::setprecision(16);

    // 1=pb be, 2=pb tr, 3=pb full be, 4=pb full tr, 5=pb asyn stability implicit, 6=pb asyn stability explicit

    if(parameter>4) asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,false,projection_rigidity);
    mass_analytic=(T)7.3;stiffness_analytic=(T)1.9e3*stiffness_multiplier;damping_analytic=damping_multiplier;restlength_analytic=(T)1.2;initial_displacement_analytic=(T)1.8;
    solids_parameters.use_trapezoidal_rule_for_velocities=(parameter==2 || parameter==4);
    initial_velocity_analytic=(T)2.3;

    // dynamic particle
    int p=particles.array_collection->Add_Element();
    particles.X(p)=TV(0,0,0);
    particles.mass(p)=mass_analytic;
    particles.V(p)=TV(0,0,0);
    int_lists[0].Append(p);

    // static particle
    p=particles.array_collection->Add_Element();
    particles.X(p)=TV(-restlength_analytic,0,0);
    particles.mass(p)=mass_analytic;
    particles.V(p)=TV(0,0,0);
    int_lists[0].Append(p);

    particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);

    DEFORMABLES_FORCES<TV>* force;
    segmented_curve=SEGMENTED_CURVE<TV>::Create(particles);
    segmented_curve->mesh.elements.Append(VECTOR<int,2>(1,2));
    segmented_curve->Update_Number_Nodes();
    deformable_body_collection.deformable_geometry.Add_Structure(segmented_curve);
    force=Create_Edge_Springs(particles,segmented_curve->mesh,stiffness_analytic,damping_analytic);
    if(parameter==3 || parameter==4) force->use_implicit_velocity_independent_forces=true;
    if(parameter<=4) solid_body_collection.Add_Force(force);
    else if(asynchronous_evolution){
        ARRAY<int>* affected_particles=new ARRAY<int>(segmented_curve->mesh.elements.Flattened());
        DEFORMABLES_FORCES<TV>* force_dummy_zero=new DEFORMABLE_GRAVITY<TV>(particles,affected_particles,0);
        substeps_per_frame=3;
        asynchronous_evolution->use_velocity_averaging=false;
        asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=true;
        if(parameter==5) asynchronous_evolution->Add_Coarsescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,true,false);
        if(parameter==6) asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,false);
        if(parameter==7){
            asynchronous_evolution->Add_Coarsescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,true,false);
            force=Create_Edge_Springs(particles,segmented_curve->mesh,stiffness_analytic,damping_analytic);
            asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,false);}
        if(parameter==8||parameter==11){
            asynchronous_evolution->use_velocity_averaging=true;
            if(parameter==8) asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=true;
            else asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=false;
            asynchronous_evolution->Add_Finescale_Force(force_dummy_zero,*affected_particles,empty_list,false);
            asynchronous_evolution->Add_Coarsescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,true,false);}
        if(parameter==9||parameter==12){
            asynchronous_evolution->use_velocity_averaging=true;
            if(parameter==9) asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=true;
            else asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=false;
            asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,false);
            asynchronous_evolution->Add_Coarsescale_Force(force_dummy_zero,*affected_particles,empty_list,true,false);}
        if(parameter==10||parameter==13){
            asynchronous_evolution->use_velocity_averaging=true;
            if(parameter==10) asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=true;
            else asynchronous_evolution->use_squared_ratio_for_implicit_velocity_independent_scale=false;
            asynchronous_evolution->Add_Coarsescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,true,false);
            force=Create_Edge_Springs(particles,segmented_curve->mesh,stiffness_analytic,damping_analytic);
            asynchronous_evolution->Add_Finescale_Force(force,ARRAY<int>(segmented_curve->mesh.elements.Flattened()),empty_list,false);}}

    particles.X(1)=TV(initial_displacement_analytic,0,0);
    particles.V(1)=TV(initial_velocity_analytic,0,0);
}
//#####################################################################
// Function Analytic_Solution_BE
//#####################################################################
VECTOR<T,2> Analytic_Solution_BE(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
  vs=v+(dt/2)*(a+ax*x+av*vs)
  xx=x+dt*vs
  vv=v+dt*(a+ax*xx+av*vv)
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T xx=-(2*x-x*av+2*v+a+ax*x)/(-2+av);
    T vv=-(-2*v+v*av-2*a+a*av-2*ax*x+ax*x*av-2*ax*v-ax*a-ax*ax*x)/(2-3*av+av*av);
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_TR
//#####################################################################
VECTOR<T,2> Analytic_Solution_TR(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
  vs=v+(dt/2)*(a+ax*x+av*vs)
  xx=x+dt*vs
  vss=v+(dt/2)*(a+ax*(x+xx)/2+av*vss)
  vv=vss*2-v
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T xx=-(2*x-x*av+2*v+a+ax*x)/(-2+av);
    T vv=-(-4*v-4*a-4*ax*x-2*v*ax-a*ax-ax*ax*x+2*av*a+2*av*ax*x+v*av*av)/(4-4*av+av*av);
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Full_BE
//#####################################################################
VECTOR<T,2> Analytic_Solution_Full_BE(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
  vs=v+(dt/2)*(a+ax*xh+av*vs)
  xh=x+(dt/2)*vs
  xx=x+dt*vs
  vv=v+dt*(a+ax*xs+av*vv)
  xs=xx+dt*vv
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T xx=-(4*x+x*ax-2*av*x+4*v+2*a)/(-4+ax+2*av);
    T vv=-(-4*v-4*a-3*v*ax-a*ax+2*av*v+2*av*a-4*ax*x-x*ax*ax+2*av*x*ax)/(-6*av+2*av*av+3*ax*av-5*ax+4+ax*ax);
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Full_TR
//#####################################################################
VECTOR<T,2> Analytic_Solution_Full_TR(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
  vs=v+(dt/2)*(a+ax*xh+av*vs)
  xh=x+(dt/2)*vs
  xx=x+dt*vs
  vss=v+(dt/2)*(a+ax*xs+av*vss)
  xs=(x+xx)/2+(dt/2)*vss
  vv=vss*2-v
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T den=-4+ax+2*av;
    T xx=4*x+x*ax-2*av*x+4*v+2*a;
    T vv=-16*v-16*a+8*av*a-16*ax*x+8*av*x*ax-8*v*ax+4*v*ax*av+4*v*av*av+v*ax*ax;
    return VECTOR<T,2>(xx,vv/dt/den)/-den;
}
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Implicit*dt
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Implicit(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
    f=0
    fv=0
    fx=0
    g=a
    gv=av
    gx=ax
    dts=dt/3
    r=3
    t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
    t1txx=x+dts*t1tvs
    t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
    t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
    t2txx=t1txx+dts*t2tvs
    t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
    t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
    t3txbepu=x+dts*r*t3tvbepu
    t3tvpu=(t2tvv+t3tvbepu)/2
    xx=t2txx+t3tvpu*dts
    vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
    t3txbevu=xx+dts*r*vv
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T den=6*ax*ax+12*ax*av+6*av*av-12*ax-12*av+6;
    T xx=5*ax*ax*v+5*ax*ax*x+10*ax*av*v-a*ax+11*ax*av*x+5*av*av*v-a*av+6*av*av*x-11*ax*v-11*ax*x-11*av*v+a-12*av*x+6*v+6*x;
    T vv=-5*ax*ax*v-5*ax*ax*x-5*ax*av*v-5*a*ax-6*ax*av*x-6*a*av+6*ax*x-6*av*v+6*a+6*v;

    return VECTOR<T,2>(xx,vv/dt)/den;
}
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Explicit
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Explicit(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
   f=a
   fv=av
   fx=ax
   g=0
   gv=0
   gx=0
   dts=dt/3
   r=3
   t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
   t1txx=x+dts*t1tvs
   t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
   t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
   t2txx=t1txx+dts*t2tvs
   t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
   t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
   t3txbepu=x+dts*r*t3tvbepu
   t3tvpu=(t2tvv+t3tvbepu)/2
   xx=t2txx+t3tvpu*dts
   vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
   t3txbevu=xx+dts*r*vv
*/

    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T den=sqr((av-6)*sqr(av-3));
    T xx=-(av-3)*(-69984*av*x+36450*av*av*x+81*ax*ax*a+26244*ax*x+52488*x+486*ax*ax*v-9234*av*av*av*x+14094*v*av*av-43740*av*v+12*av*av*ax*ax*v+2*av*av*ax*ax*a+2*av*av*ax*ax*ax*x-18*av*av*av*ax*a-18*av*av*av*ax*ax*x-72*av*av*av*ax*v+54*av*av*av*av*ax*x+54*av*av*av*av*a+1134*av*av*av*av*x-1026*a*av*av*av+81*ax*ax*ax*x+26244*a+7209*av*av*a-54*av*av*av*av*av*x+52488*v+10692*ax*v+2754*ax*a+2754*ax*ax*x-22356*ax*x*av-22356*a*av-162*ax*ax*av*v-27*ax*ax*av*a-27*ax*ax*ax*av*x-1566*ax*av*a-1566*ax*ax*av*x-6156*ax*av*v+297*ax*av*av*a+297*ax*ax*av*av*x+1188*ax*v*av*av+7209*ax*x*av*av-1026*ax*x*av*av*av-2025*av*av*av*v+108*v*av*av*av*av);
    T vv=-(-1296*av*av*av*av*a+12636*av*av*av*a-104976*ax*v-21870*av*av*v-157464*v-62694*av*av*a-40824*ax*a-40824*ax*ax*x-3240*ax*ax*ax*x+104976*av*v-157464*a-157464*ax*x-13608*ax*ax*v+54*av*av*av*av*av*a+1458*av*av*av*v+157464*av*a+157464*av*ax*x-3240*ax*ax*a-54*ax*ax*av*av*av*av*x-54*ax*av*av*av*av*a+18*ax*ax*ax*av*av*av*x+18*ax*ax*av*av*av*a+72*ax*ax*av*av*av*v-2*ax*ax*ax*ax*av*av*x-12*ax*ax*ax*av*av*v-2*ax*ax*ax*av*av*a+75816*ax*av*v+31590*ax*av*a+31590*ax*ax*av*x-62694*av*av*ax*x+162*ax*ax*ax*av*v+27*ax*ax*ax*av*a+27*ax*ax*ax*ax*av*x+7452*ax*ax*av*v+1782*ax*ax*av*a+1782*ax*ax*ax*av*x-8991*ax*av*av*a-8991*ax*ax*av*av*x-19926*av*av*ax*v+12636*av*av*av*ax*x+1134*ax*ax*av*av*av*x-1296*ax*ax*av*av*v-315*ax*ax*av*av*a-315*ax*ax*ax*av*av*x+1134*ax*av*av*av*a+2349*av*av*av*ax*v-1296*av*av*av*av*ax*x-486*ax*ax*ax*v-81*ax*ax*ax*a-81*ax*ax*ax*ax*x-108*av*av*av*av*ax*v+54*av*av*av*av*av*ax*x);
    return VECTOR<T,2>(xx,vv/dt)/(den*54);
}
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Mixed
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Mixed(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx,const T g,const T gv,const T gx)
{
/*
    f=a/2
    fv=av/2
    fx=ax/2
    g=a/2
    gv=av/2
    gx=ax/2
    dts=dt/3
    r=3
    t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
    t1txx=x+dts*t1tvs
    t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
    t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
    t2txx=t1txx+dts*t2tvs
    t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
    t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
    t3txbepu=x+dts*r*t3tvbepu
    t3tvpu=(t2tvv+t3tvbepu)/2
    xx=t2txx+t3tvpu*dts
    vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
    t3txbevu=xx+dts*r*vv
*/
    T a=(f+g)/m*dt*dt,ax=(fx+gx)/m*dt*dt,av=(fv+gv)/m*dt;v*=dt;

    T xx=-((-1959552*av*x+711504*av*av*x-4860*ax*ax*a-279936*ax*x+1679616*x-42768*ax*ax*v-112752*av*av*av*x+322704*v*av*av+3*ax*ax*ax*ax*av*x+36*ax*ax*ax*av*v+3*ax*ax*ax*av*a-648*ax*ax*ax*v-1399680*av*v-54*ax*ax*ax*a-54*ax*ax*ax*ax*x-264*av*av*ax*ax*v-31*av*av*ax*ax*a-31*av*av*ax*ax*ax*x+81*av*av*av*ax*a+81*av*av*av*ax*ax*x+180*av*av*av*ax*v+54*av*av*av*av*ax*x+216*av*av*av*av*a+8100*av*av*av*av*x-7236*a*av*av*av-4860*ax*ax*ax*x+559872*a+84240*av*av*a-216*av*av*av*av*av*x+1679616*v-528768*ax*v-117936*ax*a-117936*ax*ax*x+31104*ax*x*av-388800*a*av+7128*ax*ax*av*v+810*ax*ax*av*a+810*ax*ax*ax*av*x+31320*ax*av*a+31320*ax*ax*av*x+114048*ax*av*v-2700*ax*av*av*a-2700*ax*ax*av*av*x-7236*ax*v*av*av+8424*ax*x*av*av-1404*ax*x*av*av*av-28512*av*av*av*v+864*v*av*av*av*av)/(2592+864*av*av-72*av*av*ax-2808*av-102*av*av*av+3*ax*av*av*av-1296*ax+540*ax*av+4*av*av*av*av)/(-12+av))/54;

    T vv=-(T)2/27*(-8424*av*av*av*av*a+124902*av*av*av*a-839808*ax*v-297432*av*av*v-2519424*v-868968*av*av*a+641520*ax*a+641520*ax*ax*x+146124*ax*ax*ax*x+2099520*av*v-2519424*a-2519424*ax*x+727056*ax*ax*v+216*av*av*av*av*av*a+11664*av*av*av*v+2729376*av*a+2729376*av*ax*x+146124*ax*ax*a-54*ax*ax*av*av*av*av*x-54*ax*av*av*av*av*a-81*ax*ax*ax*av*av*av*x-81*ax*ax*av*av*av*a-180*ax*ax*av*av*av*v+31*ax*ax*ax*ax*av*av*x+264*ax*ax*ax*av*av*v+31*ax*ax*ax*av*av*a+1527984*ax*av*v-88452*ax*av*a-88452*ax*ax*av*x-868968*av*av*ax*x-36*ax*ax*ax*ax*av*v-7452*ax*ax*ax*av*v-837*ax*ax*ax*av*a-837*ax*ax*ax*ax*av*x-134784*ax*ax*av*v-34992*ax*ax*av*a-34992*ax*ax*ax*av*x-5832*ax*av*av*a-5832*ax*ax*av*av*x-354780*av*av*ax*v+124902*av*av*av*ax*x+1377*ax*ax*av*av*av*x+7776*ax*ax*av*av*v+2826*ax*ax*av*av*a+2826*ax*ax*ax*av*av*x+1377*ax*av*av*av*a+29808*av*av*av*ax*v-8424*av*av*av*av*ax*x+48600*ax*ax*ax*v+5346*ax*ax*ax*a+5346*ax*ax*ax*ax*x-3*ax*ax*ax*ax*av*a-864*av*av*av*av*ax*v+216*av*av*av*av*av*ax*x-3*ax*ax*ax*ax*ax*av*x+54*ax*ax*ax*ax*ax*x+648*ax*ax*ax*ax*v+54*ax*ax*ax*ax*a)/(186624-186624*ax-65232*av*av*av+9252*av*av*av*av+217728*ax*av-342144*av-324*ax*ax*av*av*av+9*av*av*av*av*ax*ax-79056*av*av*ax+224208*av*av+4212*ax*ax*av*av+12528*ax*av*av*av-900*ax*av*av*av*av+24*av*av*av*av*av*ax-23328*ax*ax*av-624*av*av*av*av*av+16*av*av*av*av*av*av+46656*ax*ax);

    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Implicit
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Implicit(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
    f=0
    fv=0
    fx=0
    g=a
    gv=av
    gx=ax
    dts=dt/3
    r=3
    t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
    t1txx=x+dts*t1tvs
    t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
    t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
    t2txx=t1txx+dts*t2tvs
    t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
    t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
    t3txbepu=x+dts*t3tvbepu
    t3tvpu=(t2tvv+t3tvbepu)/2
    xx=t2txx+t3tvpu*dts
    vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
    t3txbevu=xx+dts*vv
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;

    T xx=-((18*v-5*ax*v-15*av*v+3*a-3*ax*x+18*x-18*x*av)/(-3+ax+3*av))/6;

    T vv=-((-12*ax*v+3*ax*a+3*ax*ax*x-18*v+18*av*v-18*a+18*av*a+5*ax*ax*v+15*ax*av*v-18*ax*x+18*av*ax*x)/(9-6*ax-18*av+ax*ax+6*av*ax+9*av*av))/2;

    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Explicit
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Explicit(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
    f=a
    fv=av
    fx=ax
    g=0
    gv=0
    gx=0
    dts=dt/3
    r=3
    t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
    t1txx=x+dts*t1tvs
    t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
    t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
    t2txx=t1txx+dts*t2tvs
    t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
    t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
    t3txbepu=x+dts*t3tvbepu
    t3tvpu=(t2tvv+t3tvbepu)/2
    xx=t2txx+t3tvpu*dts
    vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
    t3txbevu=xx+dts*vv
*/

    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;

    T xx=-((52488*x)+486*ax*ax*v+81*ax*ax*a+81*ax*ax*ax*x-1566*ax*ax*av*x+2*ax*ax*av*av*a+2*ax*ax*ax*av*av*x+12*ax*ax*v*av*av-162*ax*ax*av*v-27*ax*ax*av*a-27*ax*ax*ax*av*x-2025*av*av*av*v+54*av*av*av*av*a-22356*av*a+7209*av*av*a-43740*av*v+108*v*av*av*av*av-1026*av*av*av*a+14094*av*av*v-22356*av*ax*x+54*av*av*av*av*ax*x+7209*av*av*ax*x-6156*ax*av*v-1566*ax*av*a-18*ax*av*av*av*a-18*ax*ax*av*av*av*x+297*ax*av*av*a+297*ax*ax*av*av*x-72*av*av*av*ax*v+1134*x*av*av*av*av-54*x*av*av*av*av*av+36450*x*av*av-69984*x*av-9234*x*av*av*av+26244*ax*x+26244*a+1188*ax*av*av*v+52488*v-1026*av*av*av*ax*x+10692*ax*v+2754*ax*a+2754*ax*ax*x)/(162-189*av-15*av*av*av+81*av*av+av*av*av*av)/(-6+av)/54;

    T vv=-(1458*av*av*av*v+75816*ax*av*v+31590*ax*av*a-157464*v-1296*av*av*av*av*a-104976*ax*v+157464*av*a+54*av*av*av*av*av*a-157464*a-12*ax*ax*ax*av*av*v-2*ax*ax*ax*av*av*a-157464*ax*x-62694*av*av*a+104976*av*v+162*ax*ax*ax*av*v+27*ax*ax*ax*av*a+27*ax*ax*ax*ax*av*x-315*ax*ax*av*av*a-315*ax*ax*ax*av*av*x-1296*ax*ax*v*av*av+1134*ax*av*av*av*a+1134*ax*ax*av*av*av*x+72*av*av*av*ax*ax*v+18*av*av*av*ax*ax*a+18*av*av*av*ax*ax*ax*x-54*av*av*av*av*ax*a-54*av*av*av*av*ax*ax*x-108*av*av*av*av*ax*v+54*av*av*av*av*av*ax*x+7452*ax*ax*av*v+1782*ax*ax*av*a+1782*ax*ax*ax*av*x-8991*ax*av*av*a-8991*ax*ax*av*av*x-19926*ax*av*av*v+12636*av*av*av*ax*x+157464*av*ax*x-1296*av*av*av*av*ax*x+2349*av*av*av*ax*v-2*ax*ax*ax*ax*av*av*x-40824*ax*a-40824*ax*ax*x-486*ax*ax*ax*v-81*ax*ax*ax*a-81*ax*ax*ax*ax*x+12636*av*av*av*a+31590*ax*ax*av*x-62694*av*av*ax*x-3240*ax*ax*a-3240*ax*ax*ax*x-21870*av*av*v-13608*ax*ax*v)/(2916+av*av*av*av*av*av-4860*av+3321*av*av-1188*av*av*av+234*av*av*av*av-24*av*av*av*av*av)/54;

    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Mixed
//#####################################################################
VECTOR<T,2> Analytic_Solution_Async_Velocity_Averaging_Nosquared_Implicit_Mixed(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx,const T g,const T gv,const T gx)
{
/*
    f=a/2
    fv=av/2
    fx=ax/2
    g=a/2
    gv=av/2
    gx=ax/2
    dts=dt/3
    r=3
    t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
    t1txx=x+dts*t1tvs
    t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
    t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
    t2txx=t1txx+dts*t2tvs
    t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
    t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
    t3txbepu=x+dts*t3tvbepu
    t3tvpu=(t2tvv+t3tvbepu)/2
    xx=t2txx+t3tvpu*dts
    vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
    t3txbevu=xx+dts*vv
*/
    T a=(f+g)/m*dt*dt,ax=(fx+gx)/m*dt*dt,av=(fv+gv)/m*dt;v*=dt;

    T xx=-((1679616*x)-11664*ax*ax*v-1404*ax*ax*a-1404*ax*ax*ax*x+2808*ax*ax*av*x-7*ax*ax*av*av*a-7*ax*ax*ax*av*av*x-48*ax*ax*v*av*av+1512*ax*ax*av*v+198*ax*ax*av*a+198*ax*ax*ax*av*x-28512*av*av*av*v+216*av*av*av*av*a-388800*av*a+84240*av*av*a-1399680*av*v+864*v*av*av*av*av-7236*av*av*av*a+322704*av*av*v-248832*av*ax*x+162*av*av*av*av*ax*x+58968*av*av*ax*x-25920*ax*av*v+2808*ax*av*a-9*ax*av*av*av*a-9*ax*ax*av*av*av*x+108*ax*av*av*a+108*ax*ax*av*av*x-252*av*av*av*ax*v+8100*x*av*av*av*av-216*x*av*av*av*av*av+711504*x*av*av-1959552*x*av-112752*x*av*av*av+279936*ax*x+559872*a+6372*ax*av*av*v-18*ax*ax*ax*a-18*ax*ax*ax*ax*x-216*ax*ax*ax*v+ax*ax*ax*av*a+ax*ax*ax*ax*av*x+12*ax*ax*ax*av*v+1679616*v-5292*av*av*av*ax*x-62208*ax*v-24624*ax*a-24624*ax*ax*x)/(2592-2808*av-432*ax-102*av*av*av+180*av*ax+864*av*av+4*av*av*av*av+ax*av*av*av-24*ax*av*av)/(-12+av)/54;

    T vv=-2/27*(11664*av*av*av*v+1667952*ax*av*v+261468*ax*av*a-2519424*v-8424*av*av*av*av*a-1679616*ax*v+2729376*av*a+216*av*av*av*av*av*a-2519424*a+48*ax*ax*ax*av*av*v+7*ax*ax*ax*av*av*a-2519424*ax*x-868968*av*av*a+2099520*av*v-1620*ax*ax*ax*av*v-207*ax*ax*ax*av*a-207*ax*ax*ax*ax*av*x-90*ax*ax*av*av*a-90*ax*ax*ax*av*av*x-6480*ax*ax*v*av*av+5427*ax*av*av*av*a+5427*ax*ax*av*av*av*x+252*av*av*av*ax*ax*v+9*av*av*av*ax*ax*a+9*av*av*av*ax*ax*ax*x-162*av*av*av*av*ax*a-162*av*av*av*av*ax*ax*x-864*av*av*av*av*ax*v+216*av*av*av*av*av*ax*x+24624*ax*ax*av*v-3564*ax*ax*av*a-3564*ax*ax*ax*av*x-62208*ax*av*av*a-62208*ax*ax*av*av*x-360612*ax*av*av*v+124902*av*av*av*ax*x+2729376*av*ax*x-8424*av*av*av*av*ax*x+29808*av*av*av*ax*v+7*ax*ax*ax*ax*av*av*x-198288*ax*a-198288*ax*ax*x+13608*ax*ax*ax*v+1566*ax*ax*ax*a+1566*ax*ax*ax*ax*x+124902*av*av*av*a+261468*ax*ax*av*x-868968*av*av*ax*x+33372*ax*ax*a+33372*ax*ax*ax*x-297432*av*av*v+120528*ax*ax*v+216*ax*ax*ax*ax*v+18*ax*ax*ax*ax*ax*x+18*ax*ax*ax*ax*a-12*ax*ax*ax*ax*av*v-ax*ax*ax*ax*a*av-ax*ax*ax*ax*ax*x*av)/(186624-342144*av+4176*ax*av*av*av-26352*ax*av*av+224208*av*av-300*ax*av*av*av*av-36*ax*ax*av*av*av+8*ax*av*av*av*av*av+ax*ax*av*av*av*av-2592*ax*ax*av+468*ax*ax*av*av-65232*av*av*av+9252*av*av*av*av+5184*ax*ax+16*av*av*av*av*av*av-624*av*av*av*av*av-62208*ax+72576*av*ax);

    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Asyn_Stability_Implicit
//#####################################################################
VECTOR<T,2> Analytic_Solution_Asyn_Stability_Implicit(const T x,T v,const T m,const T dt,const T g,const T gv,const T gx)
{
/*
  f=0
  fv=0
  fx=0
  g=a
  gv=av
  gx=ax
  dts=dt/3
  r=3
  t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
  t1txx=x+dts*t1tvs
  t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
  t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
  t2txx=t1txx+dts*t2tvs
  t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
  t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
  t3txbepu=x+dts*r*t3tvbepu
  t3tvpu=t2tvv+dts/2*(f+fx*t2txx+fv*t3tvbepu)+dts/2*r*r*(g+gx*t3txbepu+gv*t3tvbepu)
  xx=t2txx+t3tvpu*dts
  vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
  t3txbevu=xx+dts*r*vv
*/
    T a=g/m*dt*dt,ax=gx/m*dt*dt,av=gv/m*dt;v*=dt;
    T xx=(sqr(ax)*v+sqr(ax)*x+2*ax*av*v-a*ax+3*ax*av*x+sqr(av)*v-a*av+2*sqr(av)*x-3*ax*v-3*ax*x-3*av*v+a-4*av*x+2*v+2*x)/(2*sqr(ax)+4*ax*av+2*sqr(av)-4*ax-4*av+2);
    T vv=(-sqr(ax)*v-sqr(ax)*x-ax*av*v-a*ax-2*ax*av*x-2*a*av+2*ax*x-2*av*v+2*a+2*v)/(2*sqr(ax)+4*ax*av+2*sqr(av)-4*ax-4*av+2);
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Asyn_Stability_Explicit
//#####################################################################
VECTOR<T,2> Analytic_Solution_Asyn_Stability_Explicit(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx)
{
/*
  f=a
  fv=av
  fx=ax
  g=0
  gv=0
  gx=0
  dts=dt/3
  r=3
  t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
  t1txx=x+dts*t1tvs
  t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
  t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
  t2txx=t1txx+dts*t2tvs
  t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
  t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
  t3txbepu=x+dts*r*t3tvbepu
  t3tvpu=t2tvv+dts/2*(f+fx*t2txx+fv*t3tvbepu)+dts/2*r*r*(g+gx*t3txbepu+gv*t3tvbepu)
  xx=t2txx+t3tvpu*dts
  vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
  t3txbevu=xx+dts*r*vv
*/
    T a=f/m*dt*dt,ax=fx/m*dt*dt,av=fv/m*dt;v*=dt;
    T xx=-(81*ax*ax*ax*x+10692*ax*v+52488*x+81*ax*ax*a+486*ax*ax*v+7209*av*av*a-1026*ax*x*av*av*av-162*ax*ax*av*v-27*ax*ax*av*a-27*ax*ax*ax*av*x+297*ax*av*av*a+297*ax*ax*av*av*x-6156*ax*av*v-1566*ax*av*a-1566*ax*ax*av*x+7209*av*av*ax*x+1188*ax*v*av*av-43740*av*v-1026*a*av*av*av-9234*x*av*av*av+14094*v*av*av+108*v*av*av*av*av+26244*a+26244*ax*x+54*av*av*av*av*a+54*av*av*av*av*ax*x+12*ax*ax*av*av*v+2*ax*ax*av*av*a-2025*av*av*av*v-18*ax*ax*av*av*av*x-72*v*ax*av*av*av+2*ax*ax*ax*av*av*x-18*ax*av*av*av*a-54*x*av*av*av*av*av+1134*x*av*av*av*av+52488*v+2754*ax*ax*x+2754*ax*a-69984*x*av-22356*a*av-22356*ax*x*av+36450*x*av*av)/(av*av*av*av+162-15*av*av*av+81*av*av-189*av)/(-6+av)/54;
    T vv=-((104976*av*v-1296*av*av*av*av*ax*x-1296*ax*ax*av*av*v-315*ax*ax*av*av*a-40824*ax*ax*x-157464*v+12636*av*av*av*a-13608*ax*ax*v+157464*av*a-21870*av*av*v-40824*a*ax-62694*av*av*a-1296*av*av*av*av*a-3240*ax*ax*ax*x-157464*a-104976*v*ax+1458*av*av*av*v+54*av*av*av*av*av*a-157464*ax*x-2*av*av*ax*ax*ax*ax*x-2*av*av*ax*ax*ax*a-12*av*av*ax*ax*ax*v+72*av*av*av*ax*ax*v+18*av*av*av*ax*ax*a+18*av*av*av*ax*ax*ax*x-54*av*av*av*av*ax*ax*x-108*av*av*av*av*v*ax+54*av*av*av*av*av*ax*x+7452*ax*ax*av*v+1782*ax*ax*av*a+1782*ax*ax*ax*av*x-8991*ax*av*av*a-8991*ax*ax*av*av*x-19926*v*ax*av*av-54*ax*av*av*av*av*a-486*ax*ax*ax*v-81*ax*ax*ax*ax*x-81*ax*ax*ax*a+162*ax*ax*ax*av*v+27*ax*ax*ax*av*a+27*ax*ax*ax*ax*av*x+157464*av*ax*x-3240*ax*ax*a+12636*av*av*av*ax*x+75816*ax*av*v+31590*ax*av*a+31590*ax*ax*av*x-62694*av*av*ax*x+1134*ax*ax*av*av*av*x+2349*v*ax*av*av*av-315*ax*ax*ax*av*av*x+1134*ax*av*av*av*a)/(2916+3321*av*av-1188*av*av*av+av*av*av*av*av*av-24*av*av*av*av*av+234*av*av*av*av-4860*av))/54;
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Analytic_Solution_Asyn_Stability_Mixed
//#####################################################################
VECTOR<T,2> Analytic_Solution_Asyn_Stability_Mixed(const T x,T v,const T m,const T dt,const T f,const T fv,const T fx,const T g,const T gv,const T gx)
{
/*
  f=a/2
  fv=av/2
  fx=ax/2
  g=a/2
  gv=av/2
  gx=ax/2
  dts=dt/3
  r=3
  t1tvs=v+(dts/2)*(f+fx*x+fv*t1tvs)
  t1txx=x+dts*t1tvs
  t1tvv=v+dts*(f+fx*t1txx+fv*t1tvv)
  t2tvs=t1tvv+(dts/2)*(f+fx*t1txx+fv*t2tvs)
  t2txx=t1txx+dts*t2tvs
  t2tvv=t1tvv+dts*(f+fx*t2txx+fv*t2tvv)
  t3tvbepu=t2tvv+dts*(f+fx*t2txx+fv*t3tvbepu)+dts*r*(g+gx*t3txbepu+gv*t3tvbepu)
  t3txbepu=x+dts*r*t3tvbepu
  t3tvpu=t2tvv+dts/2*(f+fx*t2txx+fv*t3tvbepu)+dts/2*r*r*(g+gx*t3txbepu+gv*t3tvbepu)
  xx=t2txx+t3tvpu*dts
  vv=t2tvv+dts*(f+fx*xx+fv*vv)+dts*r*(g+gx*t3txbevu+gv*vv)
  t3txbevu=xx+dts*r*vv
*/
    T a=(f+g)/m*dt*dt,ax=(fx+gx)/m*dt*dt,av=(fv+gv)/m*dt;v*=dt;
    T xx=-(6804*ax*ax*ax*x-248832*ax*v+1679616*x+6804*ax*ax*a-24*ax*ax*av*av*av*v+34992*ax*ax*v+86184*av*av*a-5*ax*ax*ax*av*av*av*x-5*ax*ax*av*av*av*a-1404*ax*x*av*av*av+3*ax*ax*x*av*av*av*av-7128*ax*ax*av*v-1674*ax*ax*av*a-1674*ax*ax*ax*av*x+756*ax*av*av*a+756*ax*ax*av*av*x+145152*ax*av*v-3672*ax*av*a-3672*ax*ax*av*x+10368*av*av*ax*x-22140*ax*v*av*av+av*av*ax*ax*ax*ax*x+av*av*ax*ax*ax*a+12*av*av*ax*ax*ax*v-1119744*av*v-7236*a*av*av*av-112752*x*av*av*av+276048*v*av*av+1944*ax*ax*ax*v+162*ax*ax*ax*ax*x+162*ax*ax*ax*a+864*v*av*av*av*av-324*ax*ax*ax*av*v-27*ax*ax*ax*av*a-27*ax*ax*ax*ax*av*x+839808*a+216*av*av*av*av*a+54*av*av*av*av*ax*x+672*ax*ax*av*av*v+155*ax*ax*av*av*a-26568*av*av*av*v+3*ax*av*av*av*av*a-81*ax*ax*av*av*av*x+1476*v*ax*av*av*av+155*ax*ax*ax*av*av*x-81*ax*av*av*av*a-216*x*av*av*av*av*av+8100*x*av*av*av*av-36*ax*v*av*av*av*av+1679616*v+22032*ax*ax*x+22032*ax*a-1959552*x*av-435456*a*av-15552*ax*x*av+711504*x*av*av)/(2592-2808*av-102*av*av*av+4*av*av*av*av-72*ax*av*av+864*av*av-1296*ax+3*ax*av*av*av+540*ax*av)/(-12+av)/54;
    T vv=-(T)2/27*(2099520*av*v-8424*av*av*av*av*ax*x+22680*ax*ax*av*av*v-630*ax*ax*av*av*a+361584*ax*ax*x-2519424*v+124902*av*av*av*a+447120*ax*ax*v+2729376*av*a-297432*av*av*v+361584*a*ax-868968*av*av*a-8424*av*av*av*av*a+6156*ax*ax*ax*x-2519424*a-12*ax*ax*ax*ax*av*av*v-ax*ax*ax*ax*av*av*a+27*ax*ax*ax*ax*ax*av*x+5*ax*ax*ax*ax*av*av*av*x-ax*ax*ax*ax*ax*av*av*x-839808*v*ax+11664*av*av*av*v+216*av*av*av*av*av*a-2519424*ax*x-155*av*av*ax*ax*ax*ax*x-155*av*av*ax*ax*ax*a-672*av*av*ax*ax*ax*v-1476*av*av*av*ax*ax*v+81*av*av*av*ax*ax*a+81*av*av*av*ax*ax*ax*x-54*av*av*av*av*ax*ax*x-864*av*av*av*av*v*ax+216*av*av*av*av*av*ax*x-3*ax*ax*ax*av*av*av*av*x-165888*ax*ax*av*v-7776*ax*av*av*a-7776*ax*ax*av*av*x-308124*v*ax*av*av-54*ax*av*av*av*av*a-29160*ax*ax*ax*v-6318*ax*ax*ax*ax*x-6318*ax*ax*ax*a+6804*ax*ax*ax*av*v+1647*ax*ax*ax*av*a+1647*ax*ax*ax*ax*av*x+2729376*av*ax*x-1944*ax*ax*ax*ax*v-162*ax*ax*ax*ax*ax*x+36*ax*ax*av*av*av*av*v+24*ax*ax*ax*av*av*av*v-3*ax*ax*av*av*av*av*a+5*ax*ax*ax*av*av*av*a+6156*ax*ax*a-162*ax*ax*ax*ax*a+124902*av*av*av*ax*x+1248048*ax*av*v-41796*ax*av*a-41796*ax*ax*av*x-868968*av*av*ax*x+1377*ax*ax*av*av*av*x+27864*v*ax*av*av*av-630*ax*ax*ax*av*av*x+1377*ax*av*av*av*a+27*ax*ax*ax*ax*a*av+324*ax*ax*ax*ax*av*v)/(186624-342144*av+4212*ax*ax*av*av-324*ax*ax*av*av*av+9*ax*ax*av*av*av*av-23328*ax*ax*av-186624*ax+46656*ax*ax+224208*av*av-65232*av*av*av+9252*av*av*av*av+16*av*av*av*av*av*av-624*av*av*av*av*av-79056*ax*av*av+12528*ax*av*av*av+217728*ax*av-900*ax*av*av*av*av+24*av*av*av*av*av*ax);
    return VECTOR<T,2>(xx,vv/dt);
}
//#####################################################################
// Function Initialize_Sphere_Analytic_Test
//#####################################################################
void Initialize_Sphere_Analytic_Test()
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    if(asynchronous_evolution->both_forces_particles_indices.m){
        mixed_particle_index=asynchronous_evolution->both_forces_particles_indices(1);
        initial_X_mixed_particle=particles.X(mixed_particle_index);}
    else mixed_particle_index=-1;
    ARRAY<int> implicit_particle_indices;
    asynchronous_evolution->finescale_forces_particles_map.Get_Complementary_Keys(IDENTITY_ARRAY<int>(particles.array_collection->Size()),implicit_particle_indices);
    if(implicit_particle_indices.m){
        implicit_particle_index=implicit_particle_indices(1);
        initial_X_implicit_particle=particles.X(implicit_particle_index);}
    else implicit_particle_index=-1;
    ARRAY<int> finescale_particle_indices;asynchronous_evolution->finescale_forces_particles_map.Get_Keys(finescale_particle_indices);
    if(finescale_particle_indices.m){
        finescale_particle_index=finescale_particle_indices(1);
        initial_X_finescale_particle=particles.X(finescale_particle_index);}
    else finescale_particle_index=-1;
    LOG::cout<<"finescale_particle_index="<<finescale_particle_index<<",mixed_particle_index="<<mixed_particle_index<<",implicit_particle_index="<<implicit_particle_index<<std::endl;
}
//#####################################################################
// Function Adaptive_Asynchronous
//#####################################################################
void Adaptive_Asynchronous()
{
    last_frame=300;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    fully_implicit=false;
    solids_parameters.triangle_collision_parameters.perform_self_collision=(number_of_spheres>=2);
    ground=&tests.Add_Ground();

    asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);
    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)sphere_index,(T)(2*sphere_index+1)+height_offset,0))),true,true,1000));}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(sphere_index);
        volume->mesh.Initialize_Segment_Mesh();
        ARRAY_VIEW<int> flattened_volume=volume->mesh.elements.Flattened();
        HASHTABLE<int> volume_particle_map;
        for(int i=0;i<flattened_volume.Size();i++) volume_particle_map.Set(flattened_volume(i));
        ARRAY<int> all_particles;
        volume_particle_map.Get_Keys(all_particles);
 
        T stiffness=(T)1e4*soft_surface_multiplier;
        DEFORMABLES_FORCES<TV>* force;
        ARRAY<int> affected_particle_indices=all_particles,affected_rigid_body_particle_indices;
        gravity_particles=all_particles;
        SEGMENTED_CURVE<TV>* segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
        segmented_curve->mesh.elements=volume->mesh.segment_mesh->elements;
        segmented_curve->Update_Number_Nodes();
        // fine scale forces
        force=new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,&gravity_particles);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);

        // coarse scale forces
        force=new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,&gravity_particles);
        asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,false);

        force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,false);

        force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)soft_bound_edge_damping);
        asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,true);

        // boundary particles
        TRIANGULATED_SURFACE<T>& surface=volume->Get_Boundary_Object();
        ARRAY_VIEW<int> flattened(surface.mesh.elements.Flattened());
        for(int i=0;i<flattened.Size();i++){int p=flattened(i); particle_map.Set(p,p);}
        HASHTABLE<int,int> frontier_particle_map=particle_map;
        Expand_N_Rings(1,volume->mesh.segment_mesh->elements,particle_map,frontier_particle_map);
        particle_map.Get_Keys(asynchronous_evolution->coarsescale_particles_pool);}

    asynchronous_evolution->use_adaptive=true;
}
//#####################################################################
// Function Asynchronous_Projected_Sphere
//#####################################################################
void Asynchronous_Projected_Sphere()
{
    last_frame=300;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;

    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=(number_of_spheres>=2 || model_num==2);

    if(use_async && !test_implicit_in_explicit_out) asynchronous_evolution=new ASYNCHRONOUS_EVOLUTION<TV>(solid_body_collection,solids_evolution,solids_parameters.cfl,true,projection_rigidity);
    // Add Bodies
    ground=&tests.Add_Ground();
    // Add tet volumes
    std::string model_file;
    bool use_constant_mass;
    T scale,model_scale,rotation;
    int start_point;
    switch(model_num){
        case 1://sphere
            model_file=data_directory+"/Tetrahedralized_Volumes/sphere.tet";
            use_constant_mass=true;
            scale=(T)1;model_scale=(T)1;
            rotation=(T)0;
            start_point=1;
            break;
        case 2://armadillo
            model_file=data_directory+"/Tetrahedralized_Volumes/armadillo_20K.tet";
            use_constant_mass=false;
            scale=(T)10;model_scale=(T).085;
            rotation=(T)1;
            start_point=-1;
            break;
        case 3://sphere uniform
            model_file=data_directory+"/Tetrahedralized_Volumes/sphere_uniform_52k.tet";
            use_constant_mass=true;
            scale=(T)1;model_scale=(T)1;
            rotation=(T)0;
            start_point=1;
            break;

        default: PHYSBAM_FATAL_ERROR(STRING_UTILITIES::string_sprintf("Unrecognized model number %d",test_number));}

    int middle_sphere_number=number_of_spheres/2+1;
    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        particle_map.Clean_Memory();
        primary_tetrahedralized_volumes.Append(&tests.Create_Tetrahedralized_Volume(model_file,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)(sphere_index-middle_sphere_number)*scale,(T)(2*sphere_index+start_point)*scale+height_offset,0)),TWIST<TV>(TV(),TV(0,0,rotation))),true,use_constant_mass,1000,model_scale));
        TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(sphere_index);

        if(make_it_bad<1){
            volume->mesh.Initialize_Segment_Mesh();
            SEGMENTED_CURVE<TV>* curve_temp=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
            curve_temp->mesh.elements=volume->mesh.segment_mesh->elements;
            curve_temp->Update_Number_Nodes();
            T shortest=FLT_MAX;
            for(int i=0;i<curve_temp->mesh.elements.m;i++){VECTOR<int,2> e=curve_temp->mesh.elements(i);
                T length=(deformable_body_collection.particles.X(e(2))-deformable_body_collection.particles.X(e(1))).Magnitude();
                if(length<shortest) shortest=length;}
            T largest=0;int largest_index=-1;
            for(int i=0;i<volume->mesh.elements.m;i++){
                T v=volume->Volume(i);
                if(v>largest){largest=v;largest_index=i;}}

            LOG::cout<<"ZHW before,tet number: "<<volume->mesh.elements.m<<std::endl;
            T desired_length=shortest*make_it_bad;
            int i,j,k,l;volume->mesh.elements(largest_index).Get(i,j,k,l);
            TV top=deformable_body_collection.particles.X(i);
            TV centroid=volume->Centroid(largest_index),direction=centroid-top;
            T length_scope=direction.Normalize();
            if(length_scope<desired_length) desired_length=length_scope;
            TV new_position=top+desired_length*direction;
            int new_particle=deformable_body_collection.particles.array_collection->Add_Element();
            deformable_body_collection.particles.X(new_particle)=new_position;
            deformable_body_collection.particles.mass(new_particle)=deformable_body_collection.particles.mass(i);
            VECTOR<int,4> new_tet(new_particle,j,k,l);
            volume->mesh.elements(largest_index)=new_tet;
            new_tet=VECTOR<int,4>(i,new_particle,k,l);
            volume->mesh.elements.Append(new_tet);
            new_tet=VECTOR<int,4>(i,j,new_particle,l);
            volume->mesh.elements.Append(new_tet);
            new_tet=VECTOR<int,4> (i,j,k,new_particle);
            volume->mesh.elements.Append(new_tet);
            volume->Update_Number_Nodes();
            LOG::cout<<"ZHW after,tet number: "<<volume->mesh.elements.m<<",length after made_it_bad: "<<(deformable_body_collection.particles.X(i)-deformable_body_collection.particles.X(new_particle)).Magnitude()<<",new particle index: "<<new_particle<<std::endl;}

        volume->mesh.Initialize_Segment_Mesh();
        SEGMENTED_CURVE<TV>* curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
        curve->mesh.elements=volume->mesh.segment_mesh->elements;
        curve->Update_Number_Nodes();
        deformable_body_collection.deformable_geometry.Add_Structure(curve);
        primary_segmented_curves.Append(curve);
        if(use_async || test_implicit_in_explicit_out){
            // Fill particle map
            if(treat_bottom_async){
                PARTICLES<TV>& particles=deformable_body_collection.particles;
                T top=-FLT_MAX,bottom=FLT_MAX;
                for(int i=0;i<particles.array_collection->Size();i++){int p=i;
                    if(particles.X(p)(2)>top) top=particles.X(p)(2);
                    if(particles.X(p)(2)<bottom) bottom=particles.X(p)(2);}
                for(int i=0;i<particles.array_collection->Size();i++){int p=i;
                    if(particles.X(p)(2)<=(top-bottom)*coverage_percent+bottom)
                        particle_map.Set(p,p);}}
            else if(treat_left_async){
                PARTICLES<TV>& particles=deformable_body_collection.particles;
                T right=-FLT_MAX,left=FLT_MAX;
                for(int i=0;i<particles.array_collection->Size();i++){int p=i;
                    if(particles.X(p)(1)>right) right=particles.X(p)(1);
                    if(particles.X(p)(1)<left) left=particles.X(p)(1);}
                for(int i=0;i<particles.array_collection->Size();i++){int p=i;
                    if(particles.X(p)(1)<=(right-left)*coverage_percent+left)
                        particle_map.Set(p,p);}}
            else{
                TRIANGULATED_SURFACE<T>& surface=volume->Get_Boundary_Object();
                ARRAY_VIEW<int> flattened(surface.mesh.elements.Flattened());
                for(int i=0;i<flattened.Size();i++){int p=flattened(i); particle_map.Set(p,p);}}
            
            HASHTABLE<int,int> frontier_particle_map=particle_map;
            Expand_N_Rings(number_of_rings-1,volume->mesh.segment_mesh->elements,particle_map,frontier_particle_map);
        
            // Add boundary curves
            SEGMENTED_CURVE<TV> *boundary_one_ring_segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
            Add_Mapped_Elements<2>(volume->mesh.segment_mesh->elements,boundary_one_ring_segmented_curve->mesh.elements,particle_map);
            boundary_one_ring_segmented_curve->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(boundary_one_ring_segmented_curve);
            boundary_one_ring_segmented_curves.Append(boundary_one_ring_segmented_curve);

            // Add boundary tets
            TETRAHEDRALIZED_VOLUME<T>* boundary_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
            Add_Mapped_Elements<4>(volume->mesh.elements,boundary_tetrahedralized_volume->mesh.elements,particle_map);
            boundary_tetrahedralized_volume->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(boundary_tetrahedralized_volume);
            boundary_tetrahedralized_volumes.Append(boundary_tetrahedralized_volume);

            // Add complementary curves
            SEGMENTED_CURVE<TV>* complementary_boundary_segmented_curve=SEGMENTED_CURVE<TV>::Create(deformable_body_collection.particles);
            HASHTABLE<VECTOR<int,2> > edge_mesh_element_map;
            for(int i=0;i<boundary_one_ring_segmented_curve->mesh.elements.m;i++) edge_mesh_element_map.Set(boundary_one_ring_segmented_curve->mesh.elements(i));
            for(int i=0;i<volume->mesh.segment_mesh->elements.m;i++) if(!edge_mesh_element_map.Contains(volume->mesh.segment_mesh->elements(i)))
                complementary_boundary_segmented_curve->mesh.elements.Append(volume->mesh.segment_mesh->elements(i));
            complementary_boundary_segmented_curve->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(complementary_boundary_segmented_curve);
            complementary_boundary_segmented_curves.Append(complementary_boundary_segmented_curve);

            // Add complementary boundary tet mesh for altitude, edge springs
            TETRAHEDRALIZED_VOLUME<T>* complementary_boundary_tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(deformable_body_collection.particles);
            HASHTABLE<VECTOR<int,4> > volume_mesh_element_map;
            for(int i=0;i<boundary_tetrahedralized_volume->mesh.elements.m;i++) volume_mesh_element_map.Set(boundary_tetrahedralized_volume->mesh.elements(i));
            for(int i=0;i<volume->mesh.elements.m;i++) if(!volume_mesh_element_map.Contains(volume->mesh.elements(i)))
                complementary_boundary_tetrahedralized_volume->mesh.elements.Append(volume->mesh.elements(i));
            complementary_boundary_tetrahedralized_volume->Update_Number_Nodes();
            deformable_body_collection.deformable_geometry.Add_Structure(complementary_boundary_tetrahedralized_volume);
            complementary_boundary_tetrahedralized_volumes.Append(complementary_boundary_tetrahedralized_volume);}}

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Mass_To_Parents();
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(deformable_body_collection.particles.mass);
    deformable_body_collection.particles.Compute_Auxiliary_Attributes(solid_body_collection.deformable_body_collection.soft_bindings);solid_body_collection.deformable_body_collection.soft_bindings.Set_Mass_From_Effective_Mass();

    // Add forces
    T stiffness=(T)1e4*soft_surface_multiplier*scale;
    DEFORMABLES_FORCES<TV>* force;
    ARRAY<int> affected_particle_indices,affected_rigid_body_particle_indices;
    for(int sphere_index=0;sphere_index<number_of_spheres;sphere_index++){
        if(!test_implicit_in_explicit_out){
            if(!use_async){
                TETRAHEDRALIZED_VOLUME<T>* volume=primary_tetrahedralized_volumes(sphere_index);
                force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                solid_body_collection.Add_Force(force);

                SEGMENTED_CURVE<TV>* segmented_curve=primary_segmented_curves(sphere_index);
                force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                solid_body_collection.Add_Force(force);}
            if(use_async){
                // fine scale forces
                if(boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.m){
                    boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);

                    TETRAHEDRALIZED_VOLUME<T>* volume=boundary_tetrahedralized_volumes(sphere_index);
                    force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);
                    
                    SEGMENTED_CURVE<TV>* segmented_curve=boundary_one_ring_segmented_curves(sphere_index);
                    force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    asynchronous_evolution->Add_Finescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);}

                // coarse scale forces
                if(complementary_boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.m){
                    complementary_boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);
                    
                    TETRAHEDRALIZED_VOLUME<T>* volume=complementary_boundary_tetrahedralized_volumes(sphere_index);
                    force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,false);
                    
                    SEGMENTED_CURVE<TV>* segmented_curve=complementary_boundary_segmented_curves(sphere_index);
                    force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    asynchronous_evolution->Add_Coarsescale_Force(force,affected_particle_indices,affected_rigid_body_particle_indices,coarse_fully_implicit,true);}}}
        else{
            fully_implicit=0;
            // boundary forces
            if(boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.m){
                boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);

                TETRAHEDRALIZED_VOLUME<T>* volume=boundary_tetrahedralized_volumes(sphere_index);
                force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                force->use_implicit_velocity_independent_forces=true;
                solid_body_collection.Add_Force(force);
                
                SEGMENTED_CURVE<TV>* segmented_curve=boundary_one_ring_segmented_curves(sphere_index);
                force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                force->use_implicit_velocity_independent_forces=true;
                solid_body_collection.Add_Force(force);
                
                // inside forces
                if(complementary_boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.m){
                    complementary_boundary_tetrahedralized_volumes(sphere_index)->mesh.elements.Flattened().Get_Unique(affected_particle_indices);
                    
                    TETRAHEDRALIZED_VOLUME<T>* volume=complementary_boundary_tetrahedralized_volumes(sphere_index);
                    force=Create_Altitude_Springs(*volume,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    force->use_implicit_velocity_independent_forces=false;
                    solid_body_collection.Add_Force(force);
                    
                    SEGMENTED_CURVE<TV>* segmented_curve=complementary_boundary_segmented_curves(sphere_index);
                    force=Create_Edge_Springs(*segmented_curve,(T)stiffness/(1+sqrt((T)2)),(T)overdamping_fraction);
                    force->use_implicit_velocity_independent_forces=false;
                    solid_body_collection.Add_Force(force);}}}}
        
    // add gravity
    if(use_async && !test_implicit_in_explicit_out){
        affected_particle_indices.Remove_All();
        asynchronous_evolution->Add_Finescale_Force(new DEFORMABLE_GRAVITY<TV>(deformable_body_collection.particles,
                true,(T)9.8),affected_particle_indices,affected_rigid_body_particle_indices,fine_fully_implicit);}
    else tests.Add_Gravity();
}
//#####################################################################
};
}
#endif
