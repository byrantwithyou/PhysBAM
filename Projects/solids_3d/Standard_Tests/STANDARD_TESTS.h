//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Jonathan Su.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//   1. Sphere free falling to the ground
//   2. Torus free falling to the ground
//   3. Embedded sphere free falling to the ground
//   4. Embedded torus free falling to the ground
//   5. Cloth Curtain and Ball
//   6. Quasistatic mattress deforming under gravity
//   7. Quasistatic mattress stretched
//   8. Quasistatic mattress compressed to inversion
//   9. Quasistatic mattress twisted
//  10. Cloth Curtain and Two Spheres
//  11. Two hexlinks to test if interpenetrating collision objects work
//  12. Cloth Curtain and hex link multiple level set test case
//  13. Cloth and spinning sphere
//  14. Collapsing curtain
//  15. Torus, red-green mesh with T junctions (distribution of forces, not masses)
//  16. Torus, suspended from a T junction
//  17. Torus, suspended from a T junction, with drift
//  18. Torus, suspended from a non-mesh embedded node
//  19. Two Cloth Curtains and a Sphere
//  20. Torus, with a bound copy of its surface
//  21. Cloth Curtain and segments
//  22. Free Particle Sphere and Cloth Curtain
//  23. Two Rigid Boats
//  24. Incompressible sphere free falling to the ground.
//  25. Simple cloth ramp with single point and segment
//  26. Embedded Collisions
//  27. cylinder and sphere cloth for pinching (tests protectors)
//  28. point moving down level set incline plane.
//  29. Two Cloth Curtains on a plane
//  30. Wardrobe cloth test
//  31. Falling cloths tests
//  32. Cloth twister torture test
//  33. Single point triangle
//  34. Single segment segment
//  35. Falling cloth "leaves"
//  36. strip of cloth suspendded with large gravity and mass
//  37. Two-sphere drop
//  38. Hair from a square test
//  39. <none>
//  40. <none>
//  41. Simple hanging curtain cloth test
//  42. Deformable object with springs repelling it from a cube rigid body
//  43. <none>
//  44. Testing embedded curves in a tet volume.
//  45. Testing cloth thickness
//  46. Cloth on sphere
//  47. Cloth curtain with simple bending elements
//  48. Rigid sphere falling on curtain suspended by ends
//  49. Falling cube
//  50. Two spheres fall on ground
//  51. Two bound spheres fall on ground
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_CELL.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_FACE.h>
#include <PhysBAM_Tools/Grids_Uniform/UNIFORM_GRID_ITERATOR_NODE.h>
#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Parsing/PARSE_ARGS.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Basic_Geometry/SPHERE.h>
#include <PhysBAM_Geometry/Basic_Geometry/TORUS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_ID.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Spatial_Acceleration/TETRAHEDRON_HIERARCHY.h>
#include <PhysBAM_Geometry/Tessellation/SPHERE_TESSELLATION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/FREE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/LINEAR_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/SPLINE_MODEL.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/AXIAL_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/COLLISION_AREA_PENALTY_FORCE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_ALTITUDE_SPRINGS_S3D.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SEGMENT_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/SIMPLE_TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_ELEMENTS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/TRIANGLE_BENDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_MATERIAL_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Fracture/EMBEDDED_TETRAHEDRALIZED_VOLUME_BOUNDARY_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/WIND_DRAG_3D.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Meshing/RED_GREEN_TRIANGLES.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <boost/math/special_functions/asinh.hpp>
namespace PhysBAM{

using boost::math::asinh;

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
    typedef UNIFORM_GRID_ITERATOR_CELL<TV> CELL_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_FACE<TV> FACE_ITERATOR;
    typedef UNIFORM_GRID_ITERATOR_NODE<TV> NODE_ITERATOR;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::solids_parameters;using BASE::fluids_parameters;using BASE::data_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::output_directory;
    using BASE::stream_type;using BASE::solid_body_collection;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;

    SOLIDS_STANDARD_TESTS<TV> tests;
    int parameter;

    ARRAY<TV> deformable_body_rest_positions;
    ARRAY<TV> deformable_body_rest_positions_array;

    // test 3,4
    bool use_forces_for_drift;

    // test 5,10
    int number_side_panels;
    T aspect_ratio,side_length;
    int cloth_triangles;

    // test 6
    GRID<TV> mattress_grid;

    // tests 7-9
    TWIST<TV> attachment_twist;

    // tests 10
    TV test_10_frame_velocity;

    // test 16-17
    int constrained_particle,suspended_particle,drifting_particle;

    // test 20
    SEGMENT_MESH binding_segment_mesh;

    // test 21
    ARRAY<VECTOR<int,2> > constrained_particles;

    // test 24
    T test_24_poissons_ratio;

    // test 25
    ARRAY<int> cloth_ramp_particles;

    bool no_altitude_springs;
    T stiffness_multiplier,damping_multiplier,bending_stiffness_multiplier,bending_damping_multiplier,planar_damping_multiplier,axial_bending_stiffness_multiplier,axial_bending_damping_multiplier;

    // test 26
    int number_of_boundary_refinements;

    // test 30
    T test_30_constrained_off;
    T test_30_friction_off;
    T test_30_wind_off;

    // test 32
    T test_32_wind_drag_on;
    T test_32_wind_drag_ramp_off;

    // test 38
    GRID<TV> hair_layout_grid;

    // test 44
    ARRAY<int> fixed_particles;

    bool fully_implicit;
    bool print_matrix;
    bool use_axial;
    bool substitute_springs;
    bool test_forces;
    std::string model_mesh;
    int total_loops;
    T cloth_cfl;
    bool use_be;

    bool project_nullspace;
    bool opt_noattractions;
    bool opt_nopertimesteprepulsions;
    bool opt_noglobalrepulsions;
    bool opt_noself;
    bool opt_backward;
    bool opt_velocity_prune;
    bool opt_no_friction;
    bool opt_cloth_triangles;
    bool opt_side_panels;
    bool opt_topological_hierarchy_build_frequency;
    bool opt_repulsion_pair_update_frequency;
    bool opt_residuals;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),parameter(0),use_forces_for_drift(false),
        number_side_panels(40),aspect_ratio((T)1.7),side_length((T)1.0),cloth_triangles(INT_MAX),constrained_particle(0),
        suspended_particle(0),drifting_particle(0),test_24_poissons_ratio((T).5),no_altitude_springs(false),stiffness_multiplier(1),
        damping_multiplier(1),bending_stiffness_multiplier(1),bending_damping_multiplier(1),planar_damping_multiplier(1),
        axial_bending_stiffness_multiplier(1),axial_bending_damping_multiplier(1),fully_implicit(false),print_matrix(false),
        use_axial(false),substitute_springs(false),test_forces(false),total_loops(0),cloth_cfl(4),use_be(false),
        project_nullspace(false),opt_noattractions(false),opt_nopertimesteprepulsions(false),opt_noglobalrepulsions(false),
        opt_noself(false),opt_backward(false),opt_velocity_prune(false),opt_no_friction(false),opt_cloth_triangles(false),
        opt_side_panels(false),opt_topological_hierarchy_build_frequency(false),opt_repulsion_pair_update_frequency(false),
        opt_residuals(false)
    {
    }

    ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
    {
        if(test_forces){
            solid_body_collection.deformable_body_collection.Test_Energy(time);
            solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}
    }
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Test_32_Arc_Length
//#####################################################################
T Test_32_Arc_Length(const T t)
{
    return -(1-2*t)*sqrt((T)2.5-9*t+9*sqr(t))-(T)one_sixth*asinh((T)(3-6*t))+sqrt((T)2.5)+asinh((T)3)/6; // 3/2 x^2
}
//#####################################################################
// Function Test_32_Find_Parameter
//#####################################################################
T Test_32_Find_Parameter(const T desired_s,const T t_guess)
{
    T t=t_guess,t_old=t;
    for(int i=0;i<100;i++){
        T derivative=(T)2*sqrt(10+36*(t-1)*t);
        t=t-(Test_32_Arc_Length(t)-desired_s)/derivative;
        if(abs(t_old-t)<1e-7) return t;}
    return t;
}
//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options()
{
    BASE::Register_Options();
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=4;
    solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T)1e-2;
    solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-6;
    solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
    solids_parameters.triangle_collision_parameters.repulsion_pair_attractions_threshold=-(T).3;

    parse_args->Add("-parameter",&parameter,"parameter","parameter used by multiple tests to change the parameters of the test");
    parse_args->Add("-poisson",&test_24_poissons_ratio,"value","poisson's ratio for test 24");
    parse_args->Add("-stiffen",&stiffness_multiplier,"value","stiffness multiplier for various tests");
    parse_args->Add("-stiffen_bending",&bending_stiffness_multiplier,"value","stiffness multiplier for bending springs in various cloth tests");
    parse_args->Add("-dampen_bending",&bending_damping_multiplier,"value","damping multiplier for bending springs in various cloth tests");
    parse_args->Add("-dampen_planar_bending",&planar_damping_multiplier,"value","planar damping multiplier for bending springs in various cloth tests");
    parse_args->Add("-dampen_axial_bending",&axial_bending_damping_multiplier,"value","axial damping multiplier for bending springs in various cloth tests");
    parse_args->Add("-stiffen_axial_bending",&axial_bending_stiffness_multiplier,"value","axial stiffness multiplier for bending springs in various cloth tests");
    parse_args->Add("-dampen",&damping_multiplier,"value","damping multiplier for various tests");
    parse_args->Add("-noself",&opt_noself,"disable self-collisions");
    parse_args->Add("-noglobalrepulsions",&opt_noglobalrepulsions,"disable global repulsions");
    parse_args->Add("-nopertimesteprepulsions",&opt_nopertimesteprepulsions,"disable per time step repulsions");
    parse_args->Add("-repulsion_pair_update_frequency",&solids_parameters.triangle_collision_parameters.repulsion_pair_update_frequency,
        &opt_repulsion_pair_update_frequency,"steps","How many time steps before repulsion pairs are recomputed");
    parse_args->Add("-topological_hierarchy_build_frequency",&solids_parameters.triangle_collision_parameters.topological_hierarchy_build_frequency,
        &opt_topological_hierarchy_build_frequency,"steps","How many collision steps before a hierarchy topology rebuild");
    parse_args->Add("-side_panels",&number_side_panels,&opt_side_panels,"number","Cloth side panels");
    parse_args->Add("-cloth_triangles",&cloth_triangles,&opt_cloth_triangles,"number ","Cloth number of triangles");
    parse_args->Add("-noalt",&no_altitude_springs,"don't use altitude springs");
    parse_args->Add("-backward",&opt_backward,"use backward Euler evolution");
    parse_args->Add("-printpairs",&solids_parameters.triangle_collision_parameters.output_interaction_pairs,"output interaction pairs");
    parse_args->Add("-spectrum",&solids_parameters.implicit_solve_parameters.spectral_analysis,"run spectral analysis during timestepping");
    parse_args->Add("-residuals",&opt_residuals,"print residuals during timestepping");
    parse_args->Add("-binding_springs",&use_forces_for_drift,"use binding springs for drift particles");
    parse_args->Add("-velocity_prune",&opt_velocity_prune,"turn on velocity culling for collision pairs");
    parse_args->Add("-totalloops",&total_loops,"loops","set total loops");
    parse_args->Add("-noattractions",&opt_noattractions,"disable attractions on repulsion pair inversions");
    parse_args->Add("-attractionthreshold",&solids_parameters.triangle_collision_parameters.repulsion_pair_attractions_threshold,"value","threshold for attractions of inverted repulsion pairs");
    parse_args->Add("-clothcfl",&cloth_cfl,"value","Cloth CFL");
    parse_args->Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
    parse_args->Add("-cgsolids",&solids_parameters.implicit_solve_parameters.cg_tolerance,"tol","CG tolerance for backward Euler");
    parse_args->Add("-fully_implicit",&fully_implicit,"use fully implicit forces");
    parse_args->Add("-half_fully_implicit",&solids_parameters.implicit_solve_parameters.use_half_fully_implicit,"use fully implicit forces for position update");
    parse_args->Add("-use_be",&use_be,"use backward euler");
    parse_args->Add("-print_matrix",&print_matrix,"Print Krylov matrix");
    parse_args->Add("-project_nullspace",&project_nullspace,"project out nullspace");
    parse_args->Add("-use_axial",&use_axial,"use axial bending springs");
    parse_args->Add("-extra_cg",&solids_parameters.use_projections_in_position_update,"use extra projected cg for position update");
    parse_args->Add("-no_friction",&opt_no_friction,"no friction");
    parse_args->Add("-projection_iterations",&solids_parameters.implicit_solve_parameters.cg_projection_iterations,"number","number of iterations used for projection in cg");
    parse_args->Add("-substitute_springs",&substitute_springs,"instead of finite volume, use springs");
    parse_args->Add("-solver_iterations",&solids_parameters.implicit_solve_parameters.cg_iterations,"iterations","number of iterations used for solids system");
    parse_args->Add("-test_forces",&test_forces,"use fully implicit forces");
    parse_args->Add("-model",&model_mesh,"file","use this mesh");
    parse_args->Add("-repulsion_threshold",&solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness,"value","threshold for repulsions");
    parse_args->Add("-collision_threshold",&solids_parameters.triangle_collision_parameters.collisions_collision_thickness,"value","threshold for collisions");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options()
{
    BASE::Parse_Options();
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    frame_rate=24;
    
    T cloth_clamp_fraction=(T).03; // from curtain and ball
    
    if(total_loops!=0) solids_parameters.triangle_collision_parameters.total_collision_loops=total_loops;

    if(no_altitude_springs){
        if(test_number!=5 && !RANGE<VECTOR<int,1> >(10,14).Lazy_Inside(VECTOR<int,1>(test_number)) && test_number!=19 && test_number!=22){
            LOG::cerr<<"-noalt not supported for example "<<test_number<<std::endl;exit(1);}
        output_directory+="_noalt";}
    if(stiffness_multiplier!=1){
        if(test_number!=5 && !RANGE<VECTOR<int,1> >(10,14).Lazy_Inside(VECTOR<int,1>(test_number)) && test_number!=22 && test_number !=30 && test_number!=1){
            LOG::cerr<<"-stiffen not supported for example "<<test_number<<std::endl;exit(1);}
        output_directory+=STRING_UTILITIES::string_sprintf("_stiffen%g",stiffness_multiplier);}
    if(damping_multiplier!=1){
        if(test_number!=5 && !RANGE<VECTOR<int,1> >(10,14).Lazy_Inside(VECTOR<int,1>(test_number)) && test_number!=22 && test_number!=1){
            LOG::cerr<<"-dampen not supported for example "<<test_number<<std::endl;exit(1);}
        output_directory+=STRING_UTILITIES::string_sprintf("_dampen%g",damping_multiplier);}

    solids_parameters.use_trapezoidal_rule_for_velocities=!use_be;
    solids_parameters.use_rigid_deformable_contact=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.triangle_collision_parameters.use_gauss_jacobi=false;
    solids_parameters.triangle_collision_parameters.repulsions_limiter_fraction=1;
    solids_parameters.triangle_collision_parameters.collisions_final_repulsion_limiter_fraction=.1;
    solids_parameters.triangle_collision_parameters.repulsion_pair_attractions_threshold=solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness;
    if(project_nullspace) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solid_body_collection.Print_Residuals(opt_residuals);

    //solids_parameters.triangle_collision_parameters.collisions_nonrigid_collision_attempts=20;
    //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
    //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
    

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 49:
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        case 4:
            solids_parameters.cfl=(T)5;
            break;
        case 50:
        case 51:
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            break;
        case 5:
        case 29:
        case 30:
        case 31:
        case 32:
        case 35:
        case 45:
            solids_parameters.cfl=(T)50;
            frame_rate=60;
            last_frame=(int)(7*frame_rate);
            if(test_number==30){frame_rate=120;last_frame=(int)(20*frame_rate);}
            if(test_number==31) aspect_ratio=(T)1;
            if(test_number==32){last_frame=(int)(15*frame_rate);aspect_ratio=(T)4;}
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            solids_parameters.cfl=cloth_cfl; // was 4
            solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            solids_parameters.implicit_solve_parameters.cg_iterations=200;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            if(test_number==30){
                test_30_constrained_off=(T).3;
                test_30_friction_off=(T)3;
                test_30_wind_off=(T)1000;}
            if(test_number==32){
                test_32_wind_drag_on=(T)2.83;
                test_32_wind_drag_ramp_off=(T)3.017;}
            PHYSBAM_DEBUG_PRINT("Basic settings",solids_parameters.cfl,solids_parameters.implicit_solve_parameters.cg_tolerance,solids_parameters.implicit_solve_parameters.cg_iterations);
            break;
        case 6:
        case 7:
        case 8:
        case 9:
            if(test_number==6) last_frame=(int)frame_rate; else if(test_number==9) last_frame=(int)(10*frame_rate); else last_frame=(int)(5*frame_rate);
            delete solids_evolution;
            solids_evolution=new QUASISTATIC_EVOLUTION<TV>(solids_parameters,solid_body_collection);
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            if(test_number==6) solids_parameters.newton_iterations=1; else solids_parameters.newton_iterations=10;
            solids_parameters.newton_tolerance=(T)1e-2;
            solids_parameters.use_partially_converged_result=true;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            if(test_number==7) attachment_twist.linear=TV(0,0,(T).1); else if(test_number==8) attachment_twist.linear=TV(0,0,(T)-.1);
            if(test_number==9) attachment_twist.angular=TV(0,0,(T).1);
            break;
        case 10:
            test_10_frame_velocity=TV();
        case 11:
            last_frame=(int)(14*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.deformable_object_collision_parameters.maximum_levelset_collision_projection_velocity=(T)2;
            solids_parameters.deformable_object_collision_parameters.use_spatial_partition_for_levelset_collision_objects=false;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            break;
        case 12:
        case 13:
            frame_rate=60;
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            last_frame=(int)(20*frame_rate);
            solids_parameters.cfl=cloth_cfl;
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
            LOG::cout<<"CFL="<<solids_parameters.cfl<<" cg tolerance="<<solids_parameters.implicit_solve_parameters.cg_tolerance<<std::endl;
            aspect_ratio=(T)1.0;
            solids_parameters.implicit_solve_parameters.cg_iterations=200;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.deformable_object_collision_parameters.disable_multiple_levelset_collisions=false;
            PHYSBAM_DEBUG_PRINT("Basic settings",solids_parameters.cfl,solids_parameters.implicit_solve_parameters.cg_tolerance,solids_parameters.implicit_solve_parameters.cg_iterations);
            break;
        case 14:
            last_frame=(int)(7*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            solids_evolution=new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters,solid_body_collection);
            break;
        case 15:
            break;
        case 16:
        case 17:
        case 18:
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 19:
        case 46:
            frame_rate=60;
            last_frame=(int)(3*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_iterations=300;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            break;
        case 20:
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 21:
            last_frame=(int)(3*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            break;
        case 22:
            last_frame=(int)(7*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            break;
        case 23:
            last_frame=(int)(10*frame_rate);
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            //solids_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            solids_parameters.triangle_collision_parameters.turn_off_all_collisions=true;
            break;
        case 24:{
            if(test_24_poissons_ratio!=(T).5) output_directory+=STRING_UTILITIES::string_sprintf("_p%g",test_24_poissons_ratio);
            break;}
        case 25:
            last_frame=(int)(6*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            aspect_ratio=3;
            solids_parameters.triangle_collision_parameters.turn_off_all_collisions=true;
            solid_body_collection.deformable_body_collection.triangle_repulsions.hierarchy_repulsion_thickness_multiplier=10;
            
            //solids_parameters.collisions_repulsion_clamp_fraction=(T1;
            
            //solids_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            //solids_parameters.collisions_repulsion_clamp_fraction=(T).9;
            //solids_parameters.collision_repulsion_spring_multiplier=(T)100;
            //solids_parameters.collisions_repulsion_thickness=(T).1;
            break;
        case 26:
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            number_of_boundary_refinements=3;
            break;
        case 27:
            last_frame=(int)(6*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            aspect_ratio=3;
            break;
        case 28:
            solids_parameters.use_trapezoidal_rule_for_velocities=true;
            frame_rate=30;last_frame=(int)(120*frame_rate);
            solids_parameters.cfl=(T)1;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
            //solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-4;
            break;
        case 33:
        case 34:
            frame_rate=1;
            last_frame=(int)(7*frame_rate);
            aspect_ratio=1;
            solids_parameters.cfl=(T)4.0;
            solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)3.2;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_clamp_fraction=cloth_clamp_fraction;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T).1;
            solids_parameters.triangle_collision_parameters.clamp_repulsion_thickness=false;
            break;
        case 36:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).3;
            last_frame=(int)(64*frame_rate);
            aspect_ratio=16;
            side_length=(T)1/8;
            if(!opt_side_panels) number_side_panels=1;
            solids_parameters.cfl=(T)10.0;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 37:
            last_frame=(int)(5*frame_rate);
            solids_parameters.cfl=(T)2.0;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.rigid_body_evolution_parameters.simulate_rigid_bodies=true;
            break;
        case 38:
            frame_rate=24;
            last_frame=(int)(64*frame_rate);
            solids_parameters.cfl=(T)1;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 41:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;//(int)(200*frame_rate);
            aspect_ratio=1;
            side_length=(T)1;
            if(!opt_side_panels) number_side_panels=50;
            solids_parameters.cfl=(T)1;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 47:
        case 48:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_iterations=600;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;//(int)(200*frame_rate);
            aspect_ratio=1;
            side_length=(T)1;
            if(!opt_side_panels) number_side_panels=2;
            solids_parameters.cfl=(T)1;
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 42:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;
            aspect_ratio=1;
            side_length=.5;
            if(!opt_side_panels) number_side_panels=1;
            solids_parameters.cfl=(T)10;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            break;
        case 44:
            frame_rate=24;
            last_frame=200;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            break;
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

    if(opt_backward){
        if(typeid(*solids_evolution)!=typeid(NEWMARK_EVOLUTION<TV>)){
            LOG::cerr<<"refusing to replace non-Newmark evolution with -backward"<<std::endl;exit(1);}
        solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-2;
        solids_parameters.cfl*=10;
        delete solids_evolution;
        solids_evolution=new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters,solid_body_collection);
        output_directory+="_backward";}
    
    if(opt_noself){
        if(!solids_parameters.triangle_collision_parameters.perform_self_collision){LOG::cerr<<"-noself invalid for examples without self-collisions"<<std::endl;exit(1);}
        solids_parameters.triangle_collision_parameters.perform_self_collision=false;
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        output_directory+="_noself";}
    
    if(opt_noglobalrepulsions){
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        output_directory+="_noglobalrepulsions";}
    
    if(opt_nopertimesteprepulsions){
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
        output_directory+="_nopertimesteprepulsions";}
    
    if(opt_noattractions){
        solids_parameters.triangle_collision_parameters.perform_repulsion_pair_attractions=false;
        output_directory+="_noattractions";}
    
    if(opt_repulsion_pair_update_frequency){
        output_directory+=STRING_UTILITIES::string_sprintf("_repulsionpairupdatefrequency=%d",solids_parameters.triangle_collision_parameters.repulsion_pair_update_frequency);}
    
    if(opt_velocity_prune){
        //solids_parameters.collision_pair_velocity_pruning=true;
        output_directory+="_velocityprune";}
    
    if(opt_topological_hierarchy_build_frequency){
        output_directory+=STRING_UTILITIES::string_sprintf("_topologicalhierarchybuildfrequency=%d",solids_parameters.triangle_collision_parameters.topological_hierarchy_build_frequency);}
    
    if(opt_side_panels){
        output_directory+=STRING_UTILITIES::string_sprintf("_sidepanels=%d",number_side_panels);}
    cloth_triangles=2*number_side_panels*(int)(number_side_panels*aspect_ratio);
    
    if(opt_cloth_triangles){
        if(opt_side_panels) throw std::runtime_error("Cannot have both side panels and cloth triangles set");
        if(test_number==35) cloth_triangles/=100;
        number_side_panels=(int)ceil(sqrt((T).5*(T)cloth_triangles/aspect_ratio));
        int actual_cloth_triangles=2*((int)(aspect_ratio*number_side_panels))*number_side_panels;
        PHYSBAM_DEBUG_PRINT("cloth resolution",cloth_triangles,actual_cloth_triangles,number_side_panels);
        cloth_triangles=actual_cloth_triangles;
        if(test_number==35) cloth_triangles*=100;}

    if(opt_no_friction){
        solids_parameters.no_contact_friction=true;
        solids_parameters.triangle_collision_parameters.self_collision_friction_coefficient=(T)0;}
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=true;
    bool automatically_add_to_triangle_collisions=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    switch(test_number){
        case 1:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 2:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume,solids_parameters.triangle_collision_parameters);
            tests.Add_Ground();
            break;}
        case 3:{
            EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=tests.Create_Embedded_Tetrahedralized_Volume(SPHERE<TV>(TV(),(T).9),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true);
            embedding.Update_Binding_List_From_Embedding(solid_body_collection.deformable_body_collection,false);
            automatically_add_to_triangle_collisions=false;
//            tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,soft_bindings);
            embedding.Update_Number_Nodes();
            tests.Add_Ground();
            break;}
        case 4:{
            EMBEDDED_MATERIAL_SURFACE<TV,3>& embedding=tests.Create_Embedded_Tetrahedralized_Volume(TORUS<T>(TV(),TV(0,0,1),(T).3,(T).6),RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true);
            embedding.Update_Binding_List_From_Embedding(solid_body_collection.deformable_body_collection,false);
//            tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,soft_bindings);
            embedding.Update_Number_Nodes();
            tests.Initialize_Tetrahedron_Collisions(1,embedding.embedded_object.simplicial_object,solids_parameters.triangle_collision_parameters,&embedding.material_surface);
            tests.Add_Ground();
            break;}
        case 5:
        case 45:{
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
            tests.Add_Ground();
            RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T).25,(T)0);
            body.Frame().t.z=(T).5;
            if(test_number==5) rigid_body_collection.rigid_body_particle.kinematic(body.particle_index)=true;
            if(test_number==45){
                body.is_static=true;body.Frame().t=TV(.75,.2,.5);
                deformable_body_collection.collisions.thickness_table=new HASHTABLE<int,T>();
                for(int i=0;i<particles.Size();i++) deformable_body_collection.collisions.thickness_table->Insert(i,.1);}
            break;}
        case 6:
        case 7:
        case 8:
        case 9:
            mattress_grid=GRID<TV>(TV_INT(6,10,20),RANGE<TV>(TV((T)-.25,(T)-.5,(T)-1),TV((T).25,(T).5,(T)1)));
            tests.Create_Mattress(mattress_grid,true,0,density);
            deformable_body_rest_positions=particles.X;
            break;
        case 10:
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
            tests.Add_Rigid_Body("sphere",(T).25,(T)0).Frame().t=TV((T)0,-(T).25,(T).5);
            tests.Add_Rigid_Body("sphere",(T).25,(T)0).Frame().t=TV((T)2,-(T).25,(T).5);
            break;
        case 11:{
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,0);
            RIGID_BODY<TV>& rigid_body1=tests.Add_Rigid_Body("hexlink",(T).25,(T).5);rigid_body1.Frame().t=TV((T).5,-(T).25,(T).5);rigid_body1.Frame().r=ROTATION<TV>((T)(pi/2),TV(0,0,1));
            RIGID_BODY<TV>& rigid_body2=tests.Add_Rigid_Body("hexlink",(T).25,(T).5);rigid_body2.Frame().t=TV((T)1.3,-(T).25,(T).5);rigid_body2.Frame().r=ROTATION<TV>((T)(pi/2),TV(0,0,1));
            tests.Add_Rigid_Body("sphere",(T).25,(T).5).Frame().t=TV((T)1.2,(T).25,(T).5);
            break;}
        case 12:
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).7,0))));
            tests.Add_Ground();
            tests.Add_Rigid_Body("skinnyhexlink",(T).3,(T)0).Frame().t=TV((T).45,(T).3,(T).3);
            tests.Add_Rigid_Body("skinnyhexlink",(T).3,(T)0).Frame().t=TV((T).45,(T).3,-(T).3);
            tests.Add_Rigid_Body("skinnyhexlink",(T).3,(T)0).Frame().t=TV(-(T).45,(T).3,(T).3);
            tests.Add_Rigid_Body("skinnyhexlink",(T).3,(T)0).Frame().t=TV(-(T).45,(T).3,-(T).3);
            tests.Add_Rigid_Body("sphere",(T).15,(T).5).Frame().t=TV(0,(T)1.25,0);
            break;
        case 13:
            {tests.Create_Cloth_Panel(number_side_panels,(T)1.9*side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).56,0))));
            tests.Add_Ground((T).1);
            RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T).25,(T)1);
            sphere_body.Frame().t=TV(0,(T).30,0);
            rigid_body_collection.rigid_body_particle.kinematic(sphere_body.particle_index)=true;
            RIGID_BODY<TV>& body=tests.Add_Rigid_Body("cut_pyramid",(T).1,(T).1);
            body.Frame().t=TV((T)-.65,(T).05,(T).65);body.Frame().r=ROTATION<TV>((T)-pi/2,TV(1,0,0))*body.Frame().r;
            body.is_static=true;}
            break;
        case 14:
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)1,0),ROTATION<TV>((T)pi/2,TV(1,0,0)))));
            break;
        case 15:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/red_green_torus_with_t_junctions/tetrahedralized_volume.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
            tetrahedralized_volume.Update_Number_Nodes();
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Tetrahedralized_Volumes/red_green_torus_with_t_junctions/bindings",binding_list);
            tetrahedralized_volume.mesh.boundary_mesh=new TRIANGLE_MESH();
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Tetrahedralized_Volumes/red_green_torus_with_t_junctions/boundary_mesh",*tetrahedralized_volume.mesh.boundary_mesh);
            tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume,solids_parameters.triangle_collision_parameters);
            tests.Add_Ground();
            break;}
        case 16:
        case 17:
        case 18:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(
                data_directory+"/Tetrahedralized_Volumes/red_green_torus_with_t_junctions/tetrahedralized_volume.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
            FILE_UTILITIES::Read_From_File(stream_type,data_directory+"/Tetrahedralized_Volumes/red_green_torus_with_t_junctions/bindings",binding_list);
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            constrained_particle=segmented_curve.particles.Add_Element();
            particles.mass(constrained_particle)=1;particles.X(constrained_particle)=TV(0,5,0);
            suspended_particle=4531;
            if(test_number==16)
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(constrained_particle,suspended_particle));
            else if(test_number==17){
                drifting_particle=particles.Append(particles,suspended_particle);
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(constrained_particle,drifting_particle));
                binding_segment_mesh.elements.Append(VECTOR<int,2>(suspended_particle,drifting_particle));binding_segment_mesh.Set_Number_Nodes(segmented_curve.particles.Size());}
            else if(test_number==18){
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(constrained_particle,suspended_particle));
                suspended_particle=particles.Add_Element();
                binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,suspended_particle,tetrahedralized_volume.mesh.elements(0),
                        VECTOR<T,4>((T).25*VECTOR<T,4>::All_Ones_Vector())));
                particles.X(suspended_particle)=binding_list.bindings.Last()->Embedded_Position();}
            break;}
        case 19:{
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).7,0))));
            if(parameter<2) tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).75,0))));
            tests.Add_Ground();
            if(parameter==0 || parameter==2){
                RIGID_BODY<TV>& temp_sphere=tests.Add_Rigid_Body("sphere",(T).25,(T).5);
                temp_sphere.Frame().t=TV(0,(T).25,0);
                temp_sphere.is_static=true;}
            break;}
        case 20:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/adaptive_torus_float.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,1000);
            tetrahedralized_volume.Update_Number_Nodes();tetrahedralized_volume.Initialize_Triangulated_Surface();
            tetrahedralized_volume.mesh.Initialize_Node_On_Boundary();ARRAY<bool>& node_on_boundary=*tetrahedralized_volume.mesh.node_on_boundary;
            ARRAY<int> child_particles(tetrahedralized_volume.particles.Size());
            for(int p=0;p<node_on_boundary.m;p++) if(node_on_boundary(p)){
                child_particles(p)=particles.Append(particles,p);
                particles.X(child_particles(p)).x+=1;
                binding_segment_mesh.elements.Append(VECTOR<int,2>(child_particles(p),p));}
            binding_segment_mesh.Set_Number_Nodes(particles.Size());
            TRIANGULATED_SURFACE<T>* triangulated_surface=TRIANGULATED_SURFACE<T>::Create(particles);
            for(int t=0;t<tetrahedralized_volume.triangulated_surface->mesh.elements.m;t++)
                triangulated_surface->mesh.elements.Append(VECTOR<int,3>::Map(child_particles,tetrahedralized_volume.triangulated_surface->mesh.elements(t)));
            deformable_body_collection.deformable_geometry.Add_Structure(triangulated_surface);
            automatically_add_to_collision_structures=false;
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures.Last());
            tests.Add_Ground();
            break;}
        case 21:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)1.0,(T)0))));
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            constrained_particles.Resize(2);int last_particle;
            for(int fragment(0);fragment<2;fragment++){last_particle=0;
                for(T z=(T)-.25;z<=(T).25;z+=(T).02){
                    int new_particle=segmented_curve.particles.Add_Element();
                    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_particle)=
                        static_cast<DEFORMABLE_PARTICLES<TV>&>(triangulated_surface.particles).mass(0);
                    segmented_curve.particles.X(new_particle)=TV((T)(Value(fragment)*.2-.1),(T)1.25,z);
                    if(last_particle) segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_particle,last_particle));
                    else constrained_particles(fragment+1)(0)=new_particle;
                    last_particle=new_particle;}
                constrained_particles(fragment+1)(1)=last_particle;}
            int last_particle_x,last_particle_z;last_particle_x=last_particle_z=0;
            for(T pos=(T)-.25;pos<=(T).25;pos+=(T).02){
                int new_particle_x=segmented_curve.particles.Add_Element();int new_particle_z=segmented_curve.particles.Add_Element();
                static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_particle_x)=
                    static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_particle_z)=
                        static_cast<DEFORMABLE_PARTICLES<TV>&>(triangulated_surface.particles).mass(0);
                segmented_curve.particles.X(new_particle_x)=TV(pos,(T)1.5,(T)0);segmented_curve.particles.X(new_particle_z)=TV((T)0,(T)3.0,pos);
                if(last_particle_x) segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_particle_x,last_particle_x));
                if(last_particle_z) segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_particle_z,last_particle_z));
                last_particle_x=new_particle_x;last_particle_z=new_particle_z;}
            tests.Add_Ground();
            break;}
        case 22:{
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)1.5,(T)0))));
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
            T sphere_radius=(T).25;
            for(T y=-sphere_radius;y<=sphere_radius;y+=(T).01){
                T y_magnitude=abs(y);T disc_radius=sqrt((sphere_radius*sphere_radius)-(y_magnitude*y_magnitude));
                int num_particles=(int)((disc_radius/sphere_radius)*50+1);
                for(T theta=0;theta<2*pi;theta+=2*(T)pi/num_particles){
                    T z=disc_radius*sin(theta);T x=disc_radius*cos(theta);
                    int new_particle=particles.Add_Element();
                    particles.mass(new_particle)=1;particles.X(new_particle)=TV(x,y+(T).5,z);
                    free_particles.nodes.Append(new_particle);}}
            tests.Add_Ground();
            break;}
        case 23:{
            automatically_add_to_collision_structures=false;
            FRAME<TV> frame1((TV((T)0,(T)1,0)),ROTATION<TV>::From_Euler_Angles((T)0,(T).5,(T)0)),
                frame2((TV((T)0,(T).3,0)),(ROTATION<TV>::From_Euler_Angles((T)0,(T)0,(T)0)));
            VECTOR<FRAME<TV>,2> frames(frame1,frame2); // gcc 3.3.2 needs frame1/frame2 pulled out
            for(int frame_index=0;frame_index<frames.m;frame_index++){
                RIGID_BODY<TV>& body=tests.Add_Rigid_Body("Thin_Shells/boat_hires",1,(T).3,false);body.Frame().t=frames(frame_index).t;body.Frame().r=frames(frame_index).r;
                ARRAY<int> particle_indices;
                //TRIANGULATED_SURFACE<T>& surface=*static_cast<TRIANGULATED_SURFACE<T>*>(body.simplicial_object->Append_Particles_And_Create_Copy(particles,&particle_indices));
                TRIANGULATED_SURFACE<T>& old_surface=*body.simplicial_object;TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);
                int offset=particles.Size();particles.Add_Elements(old_surface.particles.Size());
                surface.mesh.Initialize_Mesh_With_Particle_Offset(old_surface.mesh,offset);
                for(int p=0;p<old_surface.particles.Size();p++){
                    particles.X(offset+p)=old_surface.particles.X(p);particles.mass(offset+p)=(T)0;
                    binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,offset+p,rigid_body_collection,body.particle_index,particles.X(offset+p)));}
                surface.Update_Number_Nodes();//surface.density=body.Mass()/surface.Total_Area();
                deformable_body_collection.deformable_geometry.Add_Structure(&surface);}
            tests.Add_Ground();
            tests.Add_Gravity();
            break;}
        case 24:{
            solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
            solids_parameters.cfl=(T)2;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_21K.tet",
            //TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,false,1000);
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
            if(0) tests.Initialize_Tetrahedron_Collisions(1,tetrahedralized_volume,solids_parameters.triangle_collision_parameters);
            tests.Add_Ground();
            break;}
        case 25:{
            // from left, most free particles
            tests.Create_Cloth_Panel(1,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)1,(T).75,0),ROTATION<TV>((T)(-pi/8),TV(0,0,1)))),&cloth_ramp_particles);
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
            int new_free_node=particles.Add_Element();particles.X(new_free_node)=TV(0,(T)1.2,0);free_particles.nodes.Append(new_free_node);
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;particles.mass(new_free_node)=aspect_ratio*sqr(side_length)/(m*n);
            // 2nd
            {tests.Create_Cloth_Panel(1,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)1,(T).75,(T)1.5),ROTATION<TV>((T)(-pi/8),TV(0,0,1)))),&cloth_ramp_particles);
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            int new_edge_node1=segmented_curve.particles.Add_Element();int new_edge_node2=segmented_curve.particles.Add_Element();
            static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=particles.mass(new_free_node);
            segmented_curve.particles.X(new_edge_node1)=TV(0,(T)1.2,(T)1.5-(T).45*side_length);segmented_curve.particles.X(new_edge_node2)=TV(0,(T)1.2,(T)1.5+(T).45*side_length);
            segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));}
            // 3rd
            {tests.Create_Cloth_Panel(1,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)1,(T).75,(T)3),ROTATION<TV>((T)(-pi/8),TV(0,0,1)))),&cloth_ramp_particles);
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            int new_edge_node1=segmented_curve.particles.Add_Element();int new_edge_node2=segmented_curve.particles.Add_Element();
            static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=particles.mass(new_free_node);
            segmented_curve.particles.X(new_edge_node1)=TV((T)-.45*side_length,(T)1.35,(T)3);segmented_curve.particles.X(new_edge_node2)=TV((T).45*side_length,(T)1.35,(T)3);
            segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));}
            // 4th
            {tests.Create_Cloth_Panel(1,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)1,(T).75,(T)4.5),ROTATION<TV>((T)(-pi/8),TV(0,0,1)))),&cloth_ramp_particles);
            TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&surface);
            int new_tri_node1=particles.Add_Element(),new_tri_node2=particles.Add_Element(),new_tri_node3=particles.Add_Element();
            particles.mass(new_tri_node1)=particles.mass(new_tri_node2)=particles.mass(new_tri_node3)=particles.mass(new_free_node);
            particles.X(new_tri_node1)=TV(0,(T)1.2,(T)4.5-(T).45*side_length);
            particles.X(new_tri_node2)=TV(0,(T)1.2,(T)4.5+(T).45*side_length);
            particles.X(new_tri_node3)=TV((T).45*side_length,(T)1.2,(T)4.5);
            surface.mesh.elements.Append(VECTOR<int,3>(new_tri_node1,new_tri_node2,new_tri_node3));}
            // set ramps to inifinite mass
            particles.mass.Subset(cloth_ramp_particles).Fill((T)FLT_MAX);
            tests.Add_Ground();
            break;}
        case 26:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2,0))),true,true,1000);
            // initialize and refine boundary
            EMBEDDING<TV>& embedding=*new EMBEDDING<TV>(particles);
            tetrahedralized_volume.Initialize_Triangulated_Surface();
            embedding.material_surface_mesh.Initialize_Mesh(tetrahedralized_volume.triangulated_surface->mesh);
            delete tetrahedralized_volume.triangulated_surface;tetrahedralized_volume.triangulated_surface=0;
            TRIANGULATED_SURFACE<T>& triangulated_surface=embedding.material_surface;
            RED_GREEN_TRIANGLES<TV> redgreen(triangulated_surface);
            for(int level=0;level<number_of_boundary_refinements;level++) redgreen.Refine_Simplex_List(IDENTITY_ARRAY<>(triangulated_surface.mesh.elements.m));
            deformable_body_collection.deformable_geometry.Add_Structure(&embedding);
            Bind_Redgreen_Segment_Midpoints(redgreen);
            // create soft bindings for all embedded particles
//            tests.Substitute_Soft_Bindings_For_Embedded_Nodes(embedding.material_surface,soft_bindings);
            tests.Add_Ground();
            RIGID_BODY<TV>& rigid_sphere=tests.Add_Rigid_Body("sphere",(T).5,(T).5);
            rigid_sphere.is_static=true;
            break;}
        case 27:{
            RIGID_BODY<TV>& cylinder_body=tests.Add_Rigid_Body("cylinder",(T).95,(T)0);cylinder_body.Frame().t=TV(0,(T)1.5,0);
            RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T)4,(T)0);sphere_body.Frame().t=TV(0,(T)1.75,(T)5.5);
            ARRAY<int> sphere_particles,cylinder_particles;
            deformable_body_collection.deformable_geometry.Add_Structure(cylinder_body.simplicial_object->Append_Particles_And_Create_Copy(deformable_body_collection.particles,&cylinder_particles));
            deformable_body_collection.deformable_geometry.Add_Structure(sphere_body.simplicial_object->Append_Particles_And_Create_Copy(deformable_body_collection.particles,&sphere_particles));
            for(int i=0;i<cylinder_particles.m;i++) deformable_body_collection.particles.X(cylinder_particles(i))=cylinder_body.World_Space_Point(deformable_body_collection.particles.X(cylinder_particles(i)));
            for(int i=0;i<sphere_particles.m;i++) deformable_body_collection.particles.X(sphere_particles(i))=sphere_body.World_Space_Point(deformable_body_collection.particles.X(sphere_particles(i)));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*((TRIANGULATED_SURFACE<T>*)deformable_body_collection.deformable_geometry.structures(0)),density);
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*((TRIANGULATED_SURFACE<T>*)deformable_body_collection.deformable_geometry.structures(1)),density);
            automatically_add_to_collision_structures=false;
            deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
            T protection_thickness=(T)atof(getenv("PROTECTION_THICKNESS"));
            LOG::cout<<"Protection thickness "<<protection_thickness<<std::endl;
            deformable_body_collection.collisions.use_protectors=true;deformable_body_collection.collisions.protection_thickness=protection_thickness;
            deformable_body_collection.collisions.protecting_bodies_of_nodes.Resize(deformable_body_collection.particles.Size());
            for(int i=0;i<cylinder_particles.m;i++){deformable_body_collection.collisions.protecting_bodies_of_nodes(cylinder_particles(i)).Append(COLLISION_GEOMETRY_ID(0));}
            for(int i=0;i<sphere_particles.m;i++){deformable_body_collection.collisions.protecting_bodies_of_nodes(sphere_particles(i)).Append(COLLISION_GEOMETRY_ID(1));}
            break;}
        case 28:{
            // The 2nd point stops by frame 14 at position 0.18655 -.135535 0.5
            // The tet's 6th point (front most) stops by frame 40 at position 0.378322 -0.274866 1
            RIGID_BODY<TV>& ground=tests.Add_Rigid_Body("ground",(T)1,(T).4);ground.Frame().r=ROTATION<TV>((T)(-pi/5),TV(0,0,1));ground.is_static=true;
            ground.coefficient_of_friction=(T).6; // used .9 for stopping case
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
            int particle1=particles.Add_Element();free_particles.nodes.Append(particle1);
            particles.mass(particle1)=1;particles.X(particle1)=TV(0,(T).1,0);particles.V(particle1)=ground.Frame().r.Rotate(TV(0,0,0));
            int particle2=particles.Add_Element();free_particles.nodes.Append(particle2);
            T init_speed=(T)0; // used (T)1 for stopping case
            particles.mass(particle2)=1;particles.X(particle2)=TV(0,0,(T).5);particles.V(particle2)=ground.Frame().r.Rotate(TV(init_speed,0,0));
            particles.Add_Elements(4);
            for(int i=0;i<3;i++){int particle=i+2;T theta=-2*(T)pi/3*(i+1);
                //particles.X(particle)=ground.Frame().r.Rotate(TV((T).2*cos(theta),0,1.1+.2*sin(theta)));
                particles.X(particle)=ground.Frame().r.Rotate(TV((T).2*cos(theta),0,(T)1.1+(T).2*sin(theta)));
                particles.mass(particle)=(T)1;particles.V(particle)=ground.Frame().r.Rotate(TV(init_speed,0,0));
            }
            particles.X(5)=ground.Frame().r.Rotate(TV(-0,(T).2,(T)1.1));particles.mass(5)=(T)1;particles.V(5)=ground.Frame().r.Rotate(TV(init_speed,0,0));
            //particles.X(1)=ground.Frame().r.Rotate(TV(0,(T).1,(T)1.1));particles.mass(1)=(T)1;particles.V(1)=ground.Frame().r.Rotate(TV(1,0,0));
            //particles.X(2)=ground.Frame().r.Rotate(TV(0,(T).1,(T)0.9));particles.mass(2)=(T)1;particles.V(2)=ground.Frame().r.Rotate(TV(1,0,0));
            //particles.X(3)=ground.Frame().r.Rotate(TV(-.1,(T).2,1));particles.mass(3)=(T)1;particles.V(3)=ground.Frame().r.Rotate(TV(1,0,0));
            TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            tetrahedralized_volume->mesh.elements.Append(VECTOR<int,4>(2,3,4,5));
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*tetrahedralized_volume,density);
            deformable_body_collection.deformable_geometry.Add_Structure(tetrahedralized_volume); 
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,1,2),ROTATION<TV>((T)pi/2,TV(0,0,1))),
                    TWIST<TV>(ground.Frame().r.Rotate(TV(1,0,0)),typename TV::SPIN())));}
            break;
        case 29:
            {RIGID_BODY_STATE<TV> state1((FRAME<TV>((TV(0,(T).5,0))))),state2((FRAME<TV>((TV(0,1,0)))));
            tests.Create_Cloth_Panel(5,side_length,aspect_ratio,&state1);
            tests.Create_Cloth_Panel(5,side_length,aspect_ratio,&state2);
            tests.Add_Ground();}
            break;
        case 30:{
            tests.Create_Cloth_Panel(number_side_panels,(T)4,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4,(T)-.25))));
            tests.Add_Ground((T)0.2,(T).1);
            RIGID_BODY<TV>& body=tests.Add_Rigid_Body("wardrobe",(T).25,(T)10);body.Frame().t=TV(0,(T)1.5,(T).9-(T)0.194);body.Frame().r=ROTATION<TV>((T)-pi/2,TV(1,0,0));
            //RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T).2,(T)10000);sphere_body.Frame().t=TV(-.5,(T)1.5,(T).7);}
            }
            break;
        case 31:{
            tests.Create_Cloth_Panel(number_side_panels,(T)1,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)-.25,(T)2.25,0),ROTATION<TV>((T)pi/3,TV((T).4,1,(T).2).Normalized()))));
            //tests.Create_Cloth_Panel(number_side_panels,(T)1,(T)1,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).25,(T)2.5,0),ROTATION<TV>(3*(T)pi/2,TV(0,1,0)))));
            //tests.Create_Cloth_Panel(number_side_panels,(T)1,(T)1,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-.25,(T)2.75,0),ROTATION<TV>(6*(T)pi/6,TV(0,1,0)))));
            //tests.Create_Cloth_Panel(number_side_panels,(T)1,(T)1,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T).25,(T)2.0,0),ROTATION<TV>(3*(T)pi/3,TV(0,1,0)))));
            tests.Add_Ground((T)0.2);}
            //RIGID_BODY<TV>& sphere_body_left=tests.Add_Rigid_Body("sphere",(T).2,(T).05);sphere_body_left.Frame().t=TV(-.75,(T)1.5,0);
            //RIGID_BODY<TV>& sphere_body_right=tests.Add_Rigid_Body("sphere",(T).2,(T).05);sphere_body_right.Frame().t=TV((T).75,(T)1.5,0);}
            break;
        case 32:{
            side_length=(T).3;
            tests.Create_Cloth_Panel(number_side_panels,(T)side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-.25,(T)2.25,0),ROTATION<TV>(0,TV(0,1,0)))));
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;T dy=side_length/(n-1);
            T s=0,t=0;T total_s=Test_32_Arc_Length((T)1),ds=total_s/(m-1);
            for(int i=0;i<m;i++){
                VECTOR<T,2> position;
                if(i>0){t=Test_32_Find_Parameter(s+ds,t);s+=ds;}
                position=VECTOR<T,2>(2*t-1,(T).5*3*sqr(2*t-1)); // 3/2 x^2
                for(int j=0;j<n;j++) particles.X(i+m*j)=TV((T).25*position.x,(T).25*position.y,j*dy-side_length/2);}
            RIGID_BODY<TV>& rigid_body=tests.Add_Rigid_Body("skinnycyllink",(T).5,(T).3);rigid_body.Frame().t=TV(0,(T)0.075,0);rigid_body.Frame().r=ROTATION<TV>((T)pi/2,TV(1,0,0));}
            break;
        case 33:{
            FREE_PARTICLES<TV>& free_particles=*FREE_PARTICLES<TV>::Create();deformable_body_collection.deformable_geometry.Add_Structure(&free_particles);
            int new_free_node=particles.Add_Element();particles.X(new_free_node)=TV((T).225*side_length,(T).025,0);particles.V(new_free_node)=TV(0,(T)-1.5,0);
            free_particles.nodes.Append(new_free_node);
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;particles.mass(new_free_node)=aspect_ratio*sqr(side_length)/(m*n);            
            TRIANGULATED_SURFACE<T>& surface=*TRIANGULATED_SURFACE<T>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&surface);
            int new_tri_node1=particles.Add_Element(),new_tri_node2=particles.Add_Element(),new_tri_node3=particles.Add_Element();
            particles.mass(new_tri_node1)=particles.mass(new_tri_node2)=particles.mass(new_tri_node3)=particles.mass(new_free_node);
            particles.X(new_tri_node1)=TV(0,0,(T)-.45*side_length);particles.X(new_tri_node2)=TV(0,0,(T).45*side_length);particles.X(new_tri_node3)=TV((T).45*side_length,0,0);
            surface.mesh.elements.Append(VECTOR<int,3>(new_tri_node1,new_tri_node2,new_tri_node3));}
            break;
        case 34:{
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            int new_edge_node1=segmented_curve.particles.Add_Element();int new_edge_node2=segmented_curve.particles.Add_Element();
            int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node1)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve.particles).mass(new_edge_node2)=aspect_ratio*sqr(side_length)/(m*n);
            segmented_curve.particles.X(new_edge_node1)=TV((T)-.45,0,0);segmented_curve.particles.X(new_edge_node2)=TV((T).45,0,0);
            segmented_curve.mesh.elements.Append(VECTOR<int,2>(new_edge_node1,new_edge_node2));
            SEGMENTED_CURVE<TV>& segmented_curve2=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve2);
            int new_edge_node3=segmented_curve2.particles.Add_Element();int new_edge_node4=segmented_curve2.particles.Add_Element();
            static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve2.particles).mass(new_edge_node3)=static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve2.particles).mass(new_edge_node4)=aspect_ratio*sqr(side_length)/(m*n);
            segmented_curve2.particles.X(new_edge_node3)=TV(0,(T).025,(T)-.45);segmented_curve2.particles.X(new_edge_node4)=TV(0,(T).025,(T).45);
            static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve2.particles).V(new_edge_node3)=TV(0,(T)-1.5,0);static_cast<DEFORMABLE_PARTICLES<TV>&>(segmented_curve2.particles).V(new_edge_node4)=TV(0,(T)-1.5,(T)0);
            segmented_curve2.mesh.elements.Append(VECTOR<int,2>(new_edge_node3,new_edge_node4));}
            break;
        case 35:{
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness=(T).1;
            solids_parameters.triangle_collision_parameters.collisions_collision_thickness=(T)1e-5;
            RANDOM_NUMBERS<T> random;random.Set_Seed(65539);
            T thickness=(T).015,max_angle=(T)pi/3,base_height=12;
            VECTOR<int,3> counts(5,5,4);
            RANGE<TV> panel_box(0,aspect_ratio*side_length,0,0,0,side_length);panel_box=panel_box.Thickened(thickness)-panel_box.Center();
            GRID<TV> grid(counts,RANGE<TV>());
            ARRAY<FRAME<TV> ,VECTOR<int,3> > frames(grid.Domain_Indices());
            for(T cell_size=1;;cell_size+=(T).01){
                LOG::cout<<"trying cell_size = "<<cell_size<<std::endl;
                grid=GRID<TV>(counts,RANGE<TV>(TV(),TV(counts))*cell_size+TV(0,base_height,0),true);
                for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){VECTOR<int,3> I=iterator.Cell_Index();
                    frames(I).r=ROTATION<TV>::From_Rotation_Vector(max_angle*random.template Get_Vector_In_Unit_Sphere<TV>());
                    frames(I).t=iterator.Location();}
                int intersections;
                for(int attempts=0;attempts<200;attempts++){intersections=0;
                    for(FACE_ITERATOR iterator(grid,0,GRID<TV>::INTERIOR_REGION);iterator.Valid();iterator.Next()){
                        VECTOR<int,3> I=iterator.First_Cell_Index(),J=iterator.Second_Cell_Index();
                        FRAME<TV> frame=frames(I).Inverse_Times(frames(J));
                        if(ORIENTED_BOX<TV>(panel_box,frame).Intersection(panel_box)){
                            frames(I).r=ROTATION<TV>::From_Rotation_Vector(max_angle*random.template Get_Vector_In_Unit_Sphere<TV>());
                            frames(J).r=ROTATION<TV>::From_Rotation_Vector(max_angle*random.template Get_Vector_In_Unit_Sphere<TV>());
                            intersections++;}}
                    if(!intersections) break;}
                LOG::cout<<"intersections = "<<intersections<<std::endl;
                if(!intersections){LOG::cout<<"final cell size = "<<cell_size<<std::endl;break;}}
            for(CELL_ITERATOR iterator(grid);iterator.Valid();iterator.Next()){VECTOR<int,3> I=iterator.Cell_Index();
                tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(frames(I)));}
            TRIANGULATED_SURFACE<T>& merged_surface=*TRIANGULATED_SURFACE<T>::Create(particles);
            for(int s=0;s<deformable_body_collection.deformable_geometry.structures.m;s++)
                merged_surface.mesh.elements.Append_Elements(dynamic_cast<TRIANGULATED_SURFACE<T>&>(*deformable_body_collection.deformable_geometry.structures(s)).mesh.elements);
            deformable_body_collection.deformable_geometry.structures.Clean_Memory();deformable_body_collection.deformable_geometry.Add_Structure(&merged_surface);
            LOG::Stat("total elements",merged_surface.mesh.elements.m);
            tests.Add_Ground((T).3,2);
            break;}
        case 36:
            {TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0),ROTATION<TV>(-(T)pi/2,TV(0,0,1)))));
            for(int i=0;i<deformable_body_collection.particles.Size();i++) deformable_body_collection.particles.X(i).y*=(T).5;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,(T)10);}
            break;
        case 37:{
            for(int i=0;i<5;i++){
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.2*(i+1),0))),true,true,1000);}
            tests.Add_Ground();
            break;}
        case 38:
            {SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);const int strands=25,strand_segments=50;
                hair_layout_grid.Initialize(TV_INT(strands,strands,strand_segments),RANGE<TV>(TV((T)-.5,0,0),TV((T).5,1,2)));particles.Add_Elements(hair_layout_grid.Numbers_Of_Nodes().Product());
            int count=1;
            for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
                for(int ij=1;ij<hair_layout_grid.counts.z;ij++){
                    particles.X(count)=hair_layout_grid.X(i,j,ij);
                    segmented_curve.mesh.elements.Append(VECTOR<int,2>(count,count+1));count++;}
                particles.X(count)=hair_layout_grid.X(i,j,hair_layout_grid.counts.z);count++;}
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(segmented_curve,density);
            RIGID_BODY<TV>& sphere_body=tests.Add_Rigid_Body("sphere",(T)1,(T)0.15);sphere_body.Frame().t=TV(0,(T)-1.2,0);
            break;}
        case 41:
            {TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>());
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,density);
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            i=0;j=0;particles.mass(i+m*j)=FLT_MAX;i=0;j=n-1;particles.mass(i+m*j)=FLT_MAX;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,true);}
            break;
        case 47:
            {TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>());
                /*DEFORMABLE_PARTICLES<TV>* temp_particles=new DEFORMABLE_PARTICLES<TV>();
            TRIANGULATED_SURFACE<T>& base_surface=*TRIANGULATED_SURFACE<T>::Create(*temp_particles);
            temp_particles->Store_Mass();
            TRIANGLE_MESH& mesh=base_surface.mesh;
            mesh.elements.Append(VECTOR<int,3>(0,1,2));
            mesh.elements.Append(VECTOR<int,3>(3,2,1));
            //mesh.elements.Append(VECTOR<int,3>(1,3,4));
            //mesh.elements.Append(VECTOR<int,3>(2,3,5));
            mesh.Set_Number_Nodes(4);
            base_surface.particles.Add_Elements(4);
            base_surface.particles.X(0)=TV(-1,0,0);
            base_surface.particles.X(1)=TV(0,0,1);
            base_surface.particles.X(2)=TV(0,0,-1);
            base_surface.particles.X(3)=TV(1,.5,0);
            //base_surface.particles.X(4)=TV(1,0,2);
            //base_surface.particles.X(5)=TV(1,0,-2);
            TRIANGULATED_SURFACE<T>& surface=tests.Copy_And_Add_Structure(base_surface);
            delete &base_surface.particles;*/
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,density);
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            i=0;j=0;particles.mass(i+m*j)=FLT_MAX;i=0;j=n-1;particles.mass(i+m*j)=FLT_MAX;
            //i=1;for(j=0;j<n;j++) particles.mass(i+m*(j-1))=FLT_MAX;
            //particles.mass(0)=FLT_MAX;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,true);}
            break;
        case 48:
            {RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T).125,(T)0);
            body.Frame().t.y=(T).5;
            TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>());
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,density);
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
            i=0;for(j=0;j<n;j++) particles.mass(i+m*j)=FLT_MAX;
            i=m-1;for(j=0;j<n;j++) particles.mass(i+m*j)=FLT_MAX;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,true);}
            break;
        case 42: 
            {TRIANGULATED_SURFACE<T>& surface=tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,2,0))));
            RIGID_BODY<TV>* rigid_body=&tests.Add_Rigid_Body("box",1,0);rigid_body->is_static=true;
            ARRAY<int> referenced_particles;surface.mesh.elements.Flattened().Get_Unique(referenced_particles);
            SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
            for(int i=0;i<referenced_particles.m;i++){int p=referenced_particles(i);
                int rigid_p=particles.Add_Element();
                TV object_space_surface=rigid_body->implicit_object->Closest_Point_On_Boundary(rigid_body->Object_Space_Point(particles.X(p)));
                RIGID_BODY_BINDING<TV>* binding=new RIGID_BODY_BINDING<TV>(particles,rigid_p,rigid_body_collection,rigid_body->particle_index,object_space_surface);
                binding_list.Add_Binding(binding);
                binding->Clamp_To_Embedded_Position();binding->Clamp_To_Embedded_Velocity();
                segmented_curve.mesh.elements.Append(VECTOR<int,2>(p,rigid_p));
            }

            automatically_add_to_collision_structures=false;
            SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(surface,true);}
            break;
        case 44:{
            int num_volumes=3,num_particles=10;T height=5;
            RIGID_BODY<TV>& body=tests.Add_Rigid_Body("sphere",(T)1,(T)0);
            body.Set_Mass(10000);
            body.Frame().t.x=num_volumes*2+5;body.Frame().t.y=height/2;body.Twist().linear.x=-50;
            for(int n=0;n<num_volumes;n++){
                //Create curve
                SEGMENTED_CURVE<TV>& segmented_curve=*SEGMENTED_CURVE<TV>::Create(particles);deformable_body_collection.deformable_geometry.Add_Structure(&segmented_curve);
                ARRAY<int> particle_indices;
                int base_particle=0;
                for(int j=0;j<num_particles;j++){
                    int p=particles.Add_Element();particle_indices.Append(p);PHYSBAM_ASSERT(p);
                    if(j==0){fixed_particles.Append(p);base_particle=p;}
                    if(particle_indices.m>1) segmented_curve.mesh.elements.Append(VECTOR<int,2>(p,particle_indices(particle_indices.m-1)));
                    particles.X(p)=TV((n-1)*2,(float)j/(float)(num_particles-1)*height,0);
                    deformable_body_rest_positions_array.Append(particles.X(p));}
                //Create tet volume
                TETRAHEDRALIZED_VOLUME<T>* volume=TETRAHEDRALIZED_VOLUME<T>::Create();
                RANGE<TV> grid_domain;
                T thickness=(T).3;int x_edge=25;
                for(int i=0;i<particle_indices.m;i++) grid_domain.Enlarge_To_Include_Point(particles.X(particle_indices(i)));
                GRID<TV> new_grid(TV_INT(grid_domain.Thickened(2*thickness).Edge_Lengths()/(T).02),grid_domain.Thickened(2*thickness));
                ARRAY<T,VECTOR<int,3> > new_phi(new_grid.Domain_Indices());new_phi.Fill(1e10);
                for(CELL_ITERATOR iterator(new_grid);iterator.Valid();iterator.Next()){const VECTOR<int,3> &cell_index=iterator.Cell_Index();
                    for(int i=0;i<particle_indices.m;i++) new_phi(cell_index)=min(new_phi(cell_index),(iterator.Location()-deformable_body_collection.deformable_geometry.particles.X(particle_indices(i))).Magnitude()-thickness);}
                LEVELSET_IMPLICIT_OBJECT<TV> implicit(new_grid,new_phi);
                TV edges=new_grid.Domain().Edge_Lengths();
                TV edge_cells_float=TV((T)x_edge,(T)x_edge*edges.y/edges.x,(T)x_edge*edges.z/edges.x);
                VECTOR<int,3> edge_cells(edge_cells_float);
                volume->Initialize_Cube_Mesh_And_Particles(GRID<TV>(edge_cells,new_grid.Domain()));
                volume->Discard_Tetrahedrons_Outside_Implicit_Surface(implicit);
                volume->Discard_Valence_Zero_Particles_And_Renumber();
                volume->Update_Number_Nodes();
                TETRAHEDRALIZED_VOLUME<T>* volume_appended=(TETRAHEDRALIZED_VOLUME<T>*)volume->Append_Particles_And_Create_Copy(particles);
                volume_appended->Update_Number_Nodes();
                volume_appended->Initialize_Hierarchy();
                SOLIDS_STANDARD_TESTS<TV>::Set_Mass_Of_Particles(*volume,1,true);
                deformable_body_collection.deformable_geometry.Add_Structure(volume_appended);
                ARRAY<int> volume_elements;volume_appended->mesh.elements.Flattened().Get_Unique(volume_elements);
                for(int i=0;i<volume_elements.m;i++){int p=volume_elements(i);
                    if(particles.X(p).y<=particles.X(base_particle).y+1e-4){
                        fixed_particles.Append(p);
                        deformable_body_rest_positions_array.Append(particles.X(p));}}
                //Embeddings
                ARRAY<TRIPLE<int,int,TV> > bindings; 
                ARRAY<int> tets;const T tolerance=(T)1e-4;
                for(int i=0;i<particle_indices.m;i++){int p=particle_indices(i);
                    tets.Remove_All();volume_appended->hierarchy->Intersection_List(particles.X(p),tets,tolerance);bool got_bind=false;
                    for(int tt=0;tt<tets.m;tt++){int t=tets(tt);
                        TV bary=TETRAHEDRON<T>::First_Three_Barycentric_Coordinates(particles.X(p),particles.X.Subset(volume_appended->mesh.elements(t)));
                        if(bary.x>-tolerance && bary.y>-tolerance && bary.z>-tolerance && bary.x+bary.y+bary.z<(T)1+tolerance){bindings.Append(TRIPLE<int,int,TV>(p,t,bary));got_bind=true;break;}}
                    if(!got_bind){LOG::cout<<"no binding on particle "<<p<<std::endl;bindings.Append(TRIPLE<int,int,TV>(p,0,TV(0,0,0)));}}
                for(int i=0;i<bindings.m;i++){
                    if(bindings(i).y==0) continue;
                    VECTOR<int,4> nodes=volume_appended->mesh.elements(bindings(i).y);
                    binding_list.Add_Binding(new LINEAR_BINDING<TV,4>(particles,bindings(i).x,nodes,bindings(i).z));
                    int soft_bound_particle=particles.Add_Element_From_Deletion_List();
                    if(fixed_particles.Contains(bindings(i).x)){
                        fixed_particles.Append(soft_bound_particle);
                        deformable_body_rest_positions_array.Append(particles.X(bindings(i).x));}
                    soft_bindings.Add_Binding(VECTOR<int,2>(soft_bound_particle,bindings(i).x),true);}
                volume_appended->Update_Number_Nodes();}
            break;}
        case 46:{
            tests.Create_Cloth_Panel(number_side_panels,side_length,aspect_ratio,RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).7,0))));
            RIGID_BODY<TV>& tmp_sphere=tests.Add_Rigid_Body("sphere",(T).25,(T).5);
            tmp_sphere.Frame().t=TV(0,(T).25,0);
            tmp_sphere.is_static=true;
            break;}
        case 49:{
            tests.Create_Mattress(GRID<TV>(TV_INT()+10,RANGE<TV>::Centered_Box()),true,0,density);
            tests.Add_Ground();
            break;}
        case 50:{
            for(int i=0;i<parameter;i++)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere_coarse.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.1*(i+1),0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 51:{
            solids_parameters.use_rigid_deformable_contact=false;
            automatically_add_to_collision_structures=false;
            automatically_add_to_triangle_collisions=false;
            if(!parameter) parameter=2;
            for(int i=0;i<parameter;i++){
                TRIANGULATED_SURFACE<T>* surface=0,*new_surface=0;
                if(model_mesh!=""){
                    surface=TRIANGULATED_SURFACE<T>::Create();
                    FILE_UTILITIES::Read_From_File<float>(model_mesh,*surface);
                    surface->Update_Bounding_Box();
                    surface->particles.X-=surface->bounding_box->Center();
                    surface->particles.X/=surface->bounding_box->Edge_Lengths().y;
                    surface->particles.X+=TV(0,i*(T)1.2+(T).6,0);}
                else{
                    SPHERE<TV> sphere(TV(0,i*2.1+2,0),1);
                    surface=TESSELLATION::Generate_Triangles(sphere,4);}
                TETRAHEDRALIZED_VOLUME<T>* new_volume=0;
                ARRAY<int> surface_particle_map;
                tests.Create_Regular_Embedded_Surface(binding_list,soft_bindings,*surface,density,512,1e-3,surface_particle_map,&new_surface,&new_volume,false);
                deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures.Last());
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures.Last());}
            tests.Add_Ground();
            break;}
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    if(automatically_add_to_triangle_collisions) solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    //solid_body_collection.deformable_body_collection.collisions.Use_Structure_Skip_Collision_Body();
    //solid_body_collection.deformable_body_collection.collisions.Use_Structure_Skip_Collision_Body();

    // correct number nodes
    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

    // correct mass
    binding_list.Distribute_Mass_To_Parents();
    binding_list.Clear_Hard_Bound_Particles(particles.mass);
    particles.Compute_Auxiliary_Attributes(soft_bindings);
    soft_bindings.Set_Mass_From_Effective_Mass();
}
//#####################################################################
// Function Initialize_Bodies
//#####################################################################
void Initialize_Bodies() PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;

    Get_Initial_Data();

    switch(test_number){
        case 17:{
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(deformable_body_collection.particles,binding_segment_mesh,1e5);
            spring_force->Clamp_Restlength(1);
            solid_body_collection.Add_Force(spring_force);}
        case 16:
        case 18:{
            SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            solid_body_collection.Add_Force(Create_Edge_Springs(deformable_body_collection.particles,segmented_curve.mesh,(T)2e4));}
        case 1:
        case 2:
        case 15:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            break;}
        case 3:
        case 4:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<EMBEDDED_MATERIAL_SURFACE<TV,3>&>().embedded_object.simplicial_object;
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume.mesh,0));
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            if(use_forces_for_drift){
                soft_bindings.use_impulses_for_collisions.Fill(false);
                soft_bindings.Initialize_Binding_Mesh();
                solid_body_collection.Add_Force(Create_Edge_Binding_Springs(deformable_body_collection.particles,*soft_bindings.binding_mesh,(T)1e6,(T)1));}
            break;}
        case 5:
        case 10:
        case 11:
        case 12:
        case 13:
        case 14:
        case 22:
        case 45:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            if(triangulated_surface.mesh.elements.m != cloth_triangles){
                LOG::cerr<<"we got "<<triangulated_surface.mesh.elements.m<<" and expected "<<cloth_triangles<<" with side_length="<<side_length<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
            if(fully_implicit){
                if(use_axial){
                    T axial_bending_stiffness=axial_bending_stiffness_multiplier*2/(1+sqrt((T)2)),axial_bending_damping=axial_bending_damping_multiplier*8;
                    solid_body_collection.Add_Force(Create_Axial_Bending_Springs(triangulated_surface,(T).01,axial_bending_stiffness,axial_bending_damping));}
                T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
                solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
                PHYSBAM_DEBUG_PRINT("Spring stiffnesses",linear_stiffness,linear_damping,bending_stiffness,bending_damping);}
            else{
                if(!no_altitude_springs) solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_surface,stiffness_multiplier*2*4/(1+sqrt((T)2)),damping_multiplier*4));
                TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(triangulated_surface);
                bend->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,(T).1);
                bend->use_force_differential=false;
                solid_body_collection.Add_Force(bend);}
            break;}
        case 30:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i))){
                if(surface->mesh.elements.m != cloth_triangles) PHYSBAM_FATAL_ERROR();
                solid_body_collection.Add_Force(Create_Edge_Springs(*surface,linear_stiffness,linear_damping)); // were *2 and *10
                solid_body_collection.Add_Force(Create_Bending_Springs(*surface,bending_stiffness,bending_damping));
                WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(*surface,solid_body_collection.rigid_body_collection);
                solid_body_collection.Add_Force(drag);
                drag->Use_Linear_Normal_Viscosity((T).001);drag->Use_Constant_Wind(0,TV((T).001,(T).0001,(T).001));}
            break;}
        case 31:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i))){
                if(surface->mesh.elements.m != cloth_triangles) PHYSBAM_FATAL_ERROR();
                solid_body_collection.Add_Force(Create_Edge_Springs(*surface,linear_stiffness,linear_damping)); // were *2 and *10
                solid_body_collection.Add_Force(Create_Bending_Springs(*surface,bending_stiffness,bending_damping));
                WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(*surface,solid_body_collection.rigid_body_collection);
                solid_body_collection.Add_Force(drag);
                drag->Use_Linear_Normal_Viscosity(17);drag->Use_Constant_Wind(0,TV((T).001,1,(T).2));}
            break;}
        case 35:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            GRID<TV>* wind_grid=new GRID<TV>(TV_INT()+20,RANGE<TV>::Bounding_Box(particles.X));
            ARRAY<TV,VECTOR<int,3> >* wind_V=new ARRAY<TV,VECTOR<int,3> >(wind_grid->Domain_Indices());
            TV wind_center=wind_grid->Domain().Center();
            INTERPOLATION_CURVE<T,T> height_falloff;
            height_falloff.Add_Control_Point(3,0);
            height_falloff.Add_Control_Point(4,1);
            for(NODE_ITERATOR iterator(*wind_grid);iterator.Valid();iterator.Next()){VECTOR<int,3> I=iterator.Node_Index();TV X=iterator.Location();
                TV inward=wind_center-X;inward.y=0;
                (*wind_V)(I)=height_falloff.Value(X.y)*((T).2*inward+TV(0,1,0));}
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i))){
                if(surface->mesh.elements.m != 100*cloth_triangles) PHYSBAM_FATAL_ERROR();
                solid_body_collection.Add_Force(Create_Edge_Springs(*surface,linear_stiffness,linear_damping)); // were *2 and *10
                solid_body_collection.Add_Force(Create_Bending_Springs(*surface,bending_stiffness,bending_damping));
                WIND_DRAG_3D<T>* drag=new WIND_DRAG_3D<T>(*surface,solid_body_collection.rigid_body_collection);solid_body_collection.Add_Force(drag);
                drag->Use_Linear_Normal_Viscosity(17);drag->Use_Spatially_Varying_Wind(0,*wind_grid,*wind_V);}
            break;}
        case 36:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            if(triangulated_surface.mesh.elements.m != cloth_triangles){
                LOG::cerr<<"we got "<<triangulated_surface.mesh.elements.m<<" and expected "<<cloth_triangles<<" with side_length="<<side_length<<std::endl;
                PHYSBAM_FATAL_ERROR();}
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0,(T)64));
            T linear_stiffness=40,linear_damping=2;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
            //T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            //solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
            PHYSBAM_DEBUG_PRINT("Spring stiffnesses",linear_stiffness,linear_damping);}
            //if(!no_altitude_springs) solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_surface,stiffness_multiplier*2*4/(1+sqrt((T)2)),damping_multiplier*4));
            //TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(triangulated_surface);
            //bend->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,(T).1);
            //bend->use_force_differential=false;
            //solid_body_collection.Add_Force(bend);
            break;
        case 32:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            if(triangulated_surface.mesh.elements.m != cloth_triangles) PHYSBAM_FATAL_ERROR();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
            LINEAR_SPRINGS<TV>* edge_springs=Create_Edge_Springs(triangulated_surface,stiffness_multiplier*10/(1+sqrt((T)2)),
                damping_multiplier*20); // were *2 and *10
            //edge_springs->restlength*=(T).9;edge_springs->visual_restlength=edge_springs->restlength;
            //edge_springs->Set_Overdamping_Fraction(damping_multiplier*10/(1+sqrt((T)2))); // put it under more tension with cylinder
            solid_body_collection.Add_Force(edge_springs);
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,2/(1+sqrt((T)2)),(T)8));
            ETHER_DRAG<GRID<TV> >* drag=new ETHER_DRAG<GRID<TV> >(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true,0);
            solid_body_collection.Add_Force(drag);}
            break;
        case 6:
        case 7:
        case 8:
        case 9:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            if(test_number==6) solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Quasistatic_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e5,(T).45)));
            break;}
        case 19:
        case 29:{
            for(int i=0;i<((parameter>=2)?1:2);i++){
                TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>(i);
                solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
                solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,stiffness_multiplier*2/(1+sqrt((T)2)),
                    damping_multiplier*2));
                if(fully_implicit){
                    if(use_axial){
                        T axial_bending_stiffness=axial_bending_stiffness_multiplier*2/(1+sqrt((T)2)),axial_bending_damping=axial_bending_damping_multiplier*8;
                        solid_body_collection.Add_Force(Create_Axial_Bending_Springs(triangulated_surface,(T).01,axial_bending_stiffness,axial_bending_damping));}
                    T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
                    solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));}
                else{
                    if(!no_altitude_springs) solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_surface,stiffness_multiplier*2*4/(1+sqrt((T)2)),damping_multiplier*4));
                    TRIANGLE_BENDING_ELEMENTS<T>* bend=Create_Bending_Elements(triangulated_surface);
                    bend->Set_Area_Cutoff_From_Triangulated_Surface(triangulated_surface,(T).1);
                    bend->use_force_differential=false;
                    solid_body_collection.Add_Force(bend);}}
            break;}
        case 20:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            //TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            //solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_surface,new NEO_HOOKEAN_2D<T>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(deformable_body_collection.particles,binding_segment_mesh);
            spring_force->visual_restlength.Fill((T)0);spring_force->Clamp_Restlength(1);
            solid_body_collection.Add_Force(spring_force);
            break;}
        case 21:{
            SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            LINEAR_SPRINGS<TV>* spring_force=Create_Edge_Springs(segmented_curve,(T)1);
            spring_force->Clamp_Restlength((T).05);solid_body_collection.Add_Force(spring_force);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,(T)3,2));
            if(!no_altitude_springs) solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_surface,2*4/(1+sqrt((T)2)),(T)4));
            solid_body_collection.Add_Force(Create_Bending_Elements(triangulated_surface));
            break;}
        case 23:{
            break;}
        case 24:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            // fvm
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            T internal_poissons_ratio=test_24_poissons_ratio<(T).5?test_24_poissons_ratio:0;
            INCOMPRESSIBLE_FINITE_VOLUME<TV,3>* fvm=Create_Incompressible_Finite_Volume(tetrahedralized_volume);
            fvm->mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
            if(test_24_poissons_ratio<(T).5) fvm->disable_projection=true;
            solid_body_collection.Add_Force(fvm);
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e4,internal_poissons_ratio,(T).01,(T).25)));
            break;}
        case 25:
        case 33:
        case 34:{
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++){
                LINEAR_SPRINGS<TV>* spring_force=0;
                if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i)))
                    spring_force=Create_Edge_Springs(*surface,(T)1);
                else if(SEGMENTED_CURVE<TV>* curve=dynamic_cast<SEGMENTED_CURVE<TV>*>(deformable_body_collection.deformable_geometry.structures(i)))
                    spring_force=Create_Edge_Springs(*curve,(T)1);
                if(spring_force){
                    spring_force->Clamp_Restlength((T).05);
                    solid_body_collection.Add_Force(spring_force);}}
            if(test_number==25) solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));}
            break;
        case 26:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,tetrahedralized_volume.mesh,0));
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)1e6,(T).45,(T).01,(T).25),true,(T).1));
            break;}
        case 27:{
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++){
                LINEAR_SPRINGS<TV>* spring_force=0;
                if(TRIANGULATED_SURFACE<T>* surface=dynamic_cast<TRIANGULATED_SURFACE<T>*>(deformable_body_collection.deformable_geometry.structures(i))){
                    spring_force=Create_Edge_Springs(*surface,(T)1e4);
                    spring_force->Clamp_Restlength((T).05);
                    solid_body_collection.Add_Force(spring_force);}}
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));}
            break;
        case 28:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            if(TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>()){
                //solid_body_collection.Add_Force(Create_Edge_Springs(*tetrahedralized_volume,
                //    stiffness_multiplier*10000/(1+sqrt((T)2)),damping_multiplier*10));
                solid_body_collection.Add_Force(Create_Finite_Volume(*tetrahedralized_volume,new ROTATED_LINEAR<T,3>((T)1e4,(T).3,(T).1)));
            }
            if(TRIANGULATED_SURFACE<T>* triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>*>()){
                T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
                solid_body_collection.Add_Force(Create_Edge_Springs(*triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
                T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
                solid_body_collection.Add_Force(Create_Bending_Springs(*triangulated_surface,bending_stiffness,bending_damping));}
            //solid_body_collection.template Find_Force<FINITE_VOLUME_3D<T>&>().strain_measure.Initialize_Rest_State_To_Equilateral_Tetrahedrons((T).1);
            break;}
        case 37:{
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(TETRAHEDRALIZED_VOLUME<T>* structure=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(i)))
                solid_body_collection.Add_Force(Create_Finite_Volume(*structure,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            break;}
        case 38:{
            SEGMENTED_CURVE<TV>& curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            solid_body_collection.Add_Force(Create_Edge_Springs(curve,100/(1+sqrt((T)2)),(T)3));
            solid_body_collection.Add_Force(Create_Segment_Bending_Springs(curve,100/(1+sqrt((T)2)),(T)3));
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 47:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));

            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
            T axial_bending_stiffness=axial_bending_stiffness_multiplier*2/(1+sqrt((T)2)),axial_bending_damping=axial_bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Axial_Bending_Springs(triangulated_surface,(T).01,axial_bending_stiffness,axial_bending_damping));
            break;}
        case 48:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            tests.Add_Gravity();
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            //solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
            solid_body_collection.Add_Force(Create_Axial_Bending_Springs(triangulated_surface,(T).01,bending_stiffness,bending_damping));
            break;}
        case 41:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
            //INCOMPRESSIBLE_FINITE_VOLUME<TV,2>* fvm=Create_Incompressible_Finite_Volume(triangulated_surface,new SPLINE_MODEL<T,2>((T)2e4,(T)0,(T).5,(T)7));
            //fvm->max_cg_iterations=50;
            //solid_body_collection.Add_Force(fvm);
            //solid_body_collection.Add_Force(Create_Finite_Volume(triangulated_surface,new SPLINE_MODEL<T,2>((T)2e4))); // ,(T)0,(T)0.01,(T)0.25)));
            
            T linear_stiffness=stiffness_multiplier*10/(1+sqrt((T)2)),linear_damping=damping_multiplier*15;
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,linear_stiffness,linear_damping)); // were *2 and *10
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
            //solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,(T)1000,(T)4)); // were *2 and *10
            //solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,(T)1000,(T)4)); // were *2 and *10
        }
            break;
        case 42:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>();
            SEGMENTED_CURVE<TV>& segmented_curve=deformable_body_collection.deformable_geometry.template Find_Structure<SEGMENTED_CURVE<TV>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,(T)1000,(T)2)); // were *2 and *10
            LINEAR_SPRINGS<TV>* stick_springs=Create_Edge_Springs(segmented_curve,(T)1000,(T)2);
            stick_springs->visual_restlength.Fill((T)0.01);stick_springs->Clamp_Restlength((T).01);
            solid_body_collection.Add_Force(stick_springs);}
            break;
        case 44:{
            for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++) if(TETRAHEDRALIZED_VOLUME<T>* structure=dynamic_cast<TETRAHEDRALIZED_VOLUME<T>*>(deformable_body_collection.deformable_geometry.structures(i)))
                solid_body_collection.Add_Force(Create_Finite_Volume(*structure,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            break;}
        case 46:{
            TRIANGULATED_SURFACE<T>& triangulated_surface=deformable_body_collection.deformable_geometry.template Find_Structure<TRIANGULATED_SURFACE<T>&>(0);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,triangulated_surface.mesh,0));
            solid_body_collection.Add_Force(Create_Edge_Springs(triangulated_surface,stiffness_multiplier*2/(1+sqrt((T)2)),
                    damping_multiplier*2));
            //if(!no_altitude_springs) solid_body_collection.Add_Force(Create_Altitude_Springs(triangulated_surface,stiffness_multiplier*2*4/(1+sqrt((T)2)),damping_multiplier*4));
            T bending_stiffness=bending_stiffness_multiplier*2/(1+sqrt((T)2)),bending_damping=bending_damping_multiplier*8;
            solid_body_collection.Add_Force(Create_Bending_Springs(triangulated_surface,bending_stiffness,bending_damping));
            break;}
        case 49:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));
            break;}
        case 50:{
            for(int i=0;i<parameter;i++){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
                solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
                solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));}
            break;}
        case 51:{
            for(int i=0;i<parameter;i++){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
                solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
                solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,new NEO_HOOKEAN<T,3>((T)2e5,(T).45,(T).01,(T).25),true,(T).1));}
            break;}
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        switch(test_number){
            case 1:
            case 49:
            case 50:
            case 51:
            case 24:{
                VECTOR<int,3> processes_per_dimension(2,1,1);
                solid_body_collection.deformable_body_collection.mpi_solids->Simple_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,processes_per_dimension);
                //solid_body_collection.deformable_body_collection.mpi_solids->KD_Tree_Partition(solids_parameters.deformable_body_collection,particles.X.array);
                break;}
            case 5:
            case 13:
            case 30:
            case 31:
            case 32:
            case 35:
            case 47:
            case 48:
            case 41:
            case 45:{
                LOG::cout<<"particles.Size()="<<particles.Size()<<std::endl;
                //solid_body_collection.deformable_body_collection.mpi_solids->Simple_Partition(solid_body_collection,particles.X.array,VECTOR<int,3>(4,1,1));
                solid_body_collection.deformable_body_collection.mpi_solids->KD_Tree_Partition(solid_body_collection.deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,ARRAY<TV>(particles.X));
                FILE_UTILITIES::Write_To_File<T>(output_directory+"/particles_of_partition",solid_body_collection.deformable_body_collection.mpi_solids->particles_of_partition);
                break;}
            default:
                LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}}

    if(substitute_springs){
        solid_body_collection.deformable_body_collection.deformables_forces.Delete_Pointers_And_Clean_Memory();
        TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        if(tetrahedralized_volume.mesh.elements.m){
            T linear_stiffness=stiffness_multiplier*1e4,linear_damping=damping_multiplier*1;
            DEFORMABLES_FORCES<TV>* altitude_springs=Create_Altitude_Springs(tetrahedralized_volume,(T)linear_stiffness/(1+sqrt((T)2)),(T)linear_damping);
            solid_body_collection.Add_Force(altitude_springs);
            DEFORMABLES_FORCES<TV>* edge_springs=Create_Edge_Springs(tetrahedralized_volume,(T)linear_stiffness/(1+sqrt((T)2)),(T)linear_damping);
            solid_body_collection.Add_Force(edge_springs);}}

    if(fully_implicit) for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(fully_implicit) for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}
//#####################################################################
// Function Set_Particle_Is_Simulated
//#####################################################################
void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(test_number==25 || test_number==28 || test_number==33){
        // make free particles simulated
        FREE_PARTICLES<TV>& free_particles=deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>&>();
        for(int i=0;i<free_particles.nodes.Size();i++) particle_is_simulated(free_particles.nodes(i))=true;}
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
void Read_Output_Files_Solids(const int frame) PHYSBAM_OVERRIDE
{
    BASE::Read_Output_Files_Solids(frame);
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    dynamic_cast<NEWMARK_EVOLUTION<TV>&>(*solids_evolution).print_matrix=print_matrix;
    LOG::cout<<"Preprocess Frame "<<frame<<std::endl;
    if(test_number==30 && frame==1){
        RANDOM_NUMBERS<T> random;random.Set_Seed(1823);
        T perturbation_size=side_length/number_side_panels*4;
        DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        for(int p=0;p<particles.Size();p++) particles.X(p).y+=random.Get_Uniform_Number((T)0,perturbation_size);}
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(LINEAR_SPRINGS<TV>* linear_springs=solid_body_collection.template Find_Force<LINEAR_SPRINGS<TV>*>())
        linear_springs->Print_Deformation_Statistics();
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE
{
    if(test_number==30){
        if(time>test_30_wind_off){
            INTERPOLATION_CURVE<T,TV> wind_stop;
            wind_stop.Add_Control_Point(test_30_wind_off-(T).3,TV((T).001,1,(T)-.25));
            wind_stop.Add_Control_Point(test_30_wind_off,TV());
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(0);drag.Use_Constant_Wind(0,TV());}
        else if(time>test_30_friction_off){
            WIND_DRAG_3D<T>& drag=solid_body_collection.template Find_Force<WIND_DRAG_3D<T>&>();
            drag.Use_Linear_Normal_Viscosity(20);drag.Use_Constant_Wind(0,TV((T).001,1,(T)-.25));}}
    else if(test_number==32){
        if(time>test_32_wind_drag_on){
            // TODO: take out this ramp of ether drag
            INTERPOLATION_CURVE<T,T> wind_viscosity;
            wind_viscosity.Add_Control_Point(test_32_wind_drag_on,0);
            wind_viscosity.Add_Control_Point(test_32_wind_drag_ramp_off,20);
            ETHER_DRAG<GRID<TV> >& drag=solid_body_collection.template Find_Force<ETHER_DRAG<GRID<TV> >&>();
            drag.Use_Constant_Wind(wind_viscosity.Value(time));}}
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==5 || test_number==41){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;j=0;V(i+m*j)=TV();i=0;j=n-1;V(i+m*j)=TV();}
    else if(test_number==47){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;j=0;V(i+m*j)=TV();i=0;j=n-1;V(i+m*j)=TV();}
        //i=0;for(j=0;j<n;j++) V(i+m*j)=TV();}
    else if(test_number==48){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;for(j=0;j<n;j++) V(i+m*j)=TV();
        i=m-1;for(j=0;j<n;j++) V(i+m*j)=TV();}
    else if(test_number==10){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;j=0;V(i+m*j)=test_10_frame_velocity;i=0;j=n-1;V(i+m*j)=test_10_frame_velocity;}
    else if(test_number==14){
        if(velocity_time<1){
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
            i=0;j=0;V(i+m*j)=TV((T).88,0,0);i=m/3;j=0;V(i+m*j)=TV((T).30,0,0);
            i=(2*m/3);j=0;V(i+m*j)=TV(-(T).30,0,0);i=m-1;j=0;V(i+m*j)=TV(-(T).88,0,0);}
        else if(velocity_time<2){
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
            i=0;j=0;V(i+m*j)=TV(-(T).88,0,0);i=m/3;j=0;V(i+m*j)=TV(-(T).30,0,0);
            i=(2*m/3);j=0;V(i+m*j)=TV((T).30,0,0);i=m-1;j=0;V(i+m*j)=TV((T).88,0,0);}
        else{
            int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
            i=0;j=0;V(i+m*j)=TV();i=m/3;j=0;V(i+m*j)=TV();
            i=(2*m/3);j=0;V(i+m*j)=TV();i=m-1;j=0;V(i+m*j)=TV();}}
    else if(test_number==16 || test_number==17 || test_number==18) V(constrained_particle)=TV();
    else if(test_number==21){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
        i=0;j=0;V(i+m*j)=TV();i=m-1;j=0;V(i+m*j)=TV();
        i=0;j=number_side_panels;V(i+m*j)=TV();i=m-1;j=number_side_panels;V(i+m*j)=TV();
        V.Subset(constrained_particles(0)).Fill(TV());
        V.Subset(constrained_particles(1)).Fill(TV());}
    else if(test_number==22){
        FREE_PARTICLES<TV>& free_particles=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>&>();
        V.Subset(free_particles.nodes).Fill(TV());}
    else if(test_number==25){
        V.Subset(cloth_ramp_particles).Fill(TV());}
    else if(test_number==30){
        if(velocity_time < test_30_constrained_off){
            INTERPOLATION_CURVE<T,T> x;
            x.Add_Control_Point((T)0,(T)5);
            x.Add_Control_Point((T)test_30_constrained_off,(T)0);
            int j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;j=(int)((T).8*n);
            V((int)((T).33*m)+m*j).z=0;V((int)((T).33*m)+m*j).x=x.Value(velocity_time);V((int)((T).66*m)+m*j).z=0;V((int)((T).66*m)+m*j).x=-x.Value(velocity_time);
            j=(int)(.6*n);V((int)(.25*m)+m*j).z=(T)1.5*x.Value(velocity_time);V((int)(.5*m)+m*j).z=-(T)1.5*x.Value(velocity_time);V((int)(.65*m)+m*j).z=(T)1.5*x.Value(velocity_time);
            j=(int)(.5*n);V((int)(.33*m)+m*j).z=0;V((int)(.33*m)+m*j).x=x.Value(velocity_time);}}
    else if(test_number==32){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        if(velocity_time<(T)3.0){
            i=0;j=0;V(i+m*j)*=TV(0,0,(T)1);i=0;j=n-1;V(i+m*j)*=TV(0,0,(T)1);i=m-1;j=0;V(i+m*j)*=TV(0,0,(T)1);i=m-1;j=n-1;V(i+m*j)*=TV(0,0,(T)1);}
        else{
            i=0;j=0;V(i+m*j)*=TV();i=0;j=n-1;V(i+m*j)*=TV();i=m-1;j=0;V(i+m*j)*=TV();i=m-1;j=n-1;V(i+m*j)*=TV();}}
    else if(test_number==36){
        int i=0;int m=(int)(aspect_ratio*number_side_panels)+1;//,n=number_side_panels+1;
        for(int j=0;j<number_side_panels+1;j++) V(i+m*j)*=TV(0,0,0);}
    else if(test_number==38){
        for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
            int index=hair_layout_grid.counts.z*(j+i*hair_layout_grid.counts.y)+1;
            V(index)=TV();}}
    else if(test_number==44) for(int i=0;i<fixed_particles.m;i++) V(fixed_particles(i))=TV();
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==5 || test_number==10 || test_number==41){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;j=0;V(i+m*j)=TV();i=0;j=n-1;V(i+m*j)=TV();}
    if(test_number==47){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;j=0;V(i+m*j)=TV();i=0;j=n-1;V(i+m*j)=TV();}
    //i=0;for(j=0;j<n;j++) V(i+m*j)=TV();}
    if(test_number==48){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        i=0;for(j=0;j<n;j++) V(i+m*j)=TV();
        i=m-1;for(j=0;j<n;j++) V(i+m*j)=TV();}
    if(test_number==14){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
        i=0;j=0;V(i+m*j)=TV();i=(m/3);j=0;V(i+m*j)=TV();
        i=(2*m/3);j=0;V(i+m*j)=TV();i=m-1;j=0;V(i+m*j)=TV();}
    if(test_number==16 || test_number==17 || test_number==18) V(constrained_particle)=TV();
    if(test_number==21){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1;
        i=0;j=0;V(i+m*j)=TV();i=m-1;j=0;V(i+m*j)=TV();
        i=0;j=number_side_panels;V(i+m*j)=TV();i=m-1;j=number_side_panels;V(i+m*j)=TV();
        V.Subset(constrained_particles(0)).Fill(TV());
        V.Subset(constrained_particles(1)).Fill(TV());}
    else if(test_number==22){
        FREE_PARTICLES<TV>& free_particles=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<FREE_PARTICLES<TV>&>();
        V.Subset(free_particles.nodes).Fill(TV());}
    else if(test_number==25){
        V.Subset(cloth_ramp_particles).Fill(TV());}
    else if(test_number==30){
        if(velocity_time < test_30_constrained_off){
            int j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;j=(int)(.8*n);
            V((int)(.33*m)+m*j).z=0;V((int)(.33*m)+m*j).x=(T)0;V((int)(.66*m)+m*j).z=0;V((int)(.66*m)+m*j).x=0;
            j=(int)(.6*n);V((int)(.25*m)+m*j).z=0;V((int)(.5*m)+m*j).z=0;V((int)(.65*m)+m*j).z=0;
            j=(int)(.5*n);V((int)(.33*m)+m*j).z=0;V((int)(.33*m)+m*j).x=(T)0;}}
    else if(test_number==32){
        int i,j;int m=(int)(aspect_ratio*number_side_panels)+1,n=number_side_panels+1;
        if(velocity_time<(T)3.0){
            i=0;j=0;V(i+m*j)*=TV(0,0,(T)1);i=0;j=n-1;V(i+m*j)*=TV(0,0,(T)1);i=m-1;j=0;V(i+m*j)*=TV(0,0,(T)1);i=m-1;j=n-1;V(i+m*j)*=TV(0,0,(T)1);}
        else{
            i=0;j=0;V(i+m*j)*=TV();i=0;j=n-1;V(i+m*j)*=TV();i=m-1;j=0;V(i+m*j)*=TV();i=m-1;j=n-1;V(i+m*j)*=TV();}}
    else if(test_number==36){
        int i=0;int m=(int)(aspect_ratio*number_side_panels)+1;//,n=number_side_panels+1;
        for(int j=0;j<number_side_panels+1;j++) V(i+m*j)*=TV(0,0,0);}
    else if(test_number==38){
        for(int i=0;i<hair_layout_grid.counts.x;i++) for(int j=0;j<hair_layout_grid.counts.y;j++){
            int index=hair_layout_grid.counts.z*(j+i*hair_layout_grid.counts.y)+1;
            V(index)=TV();}}
    else if(test_number==44) for(int i=0;i<fixed_particles.m;i++) V(fixed_particles(i))=TV();
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    if(test_number!=6 && test_number!=7 && test_number!=8 && test_number!=9 && test_number!=44) return;
    if(test_number==44){
        ARRAY<TV>& X_save=deformable_body_rest_positions_array;
        for(int i=0;i<X_save.Size();i++) X(fixed_particles(i))=X_save(i);
        return;}
    ARRAY<TV>& X_save=deformable_body_rest_positions;
    for(int j=0,index=-1;j<mattress_grid.counts.y;j++) for(int k=0;k<mattress_grid.counts.x;k++){
        index++;
        if(test_number==6 || test_number==7 || test_number==8){X(index)=X_save(index)-time*attachment_twist.linear;X(X.Size()+index)=X_save(X.Size()+index)+time*attachment_twist.linear;}
        else{
            ROTATION<TV> orientation,opposite_orientation;
            TV end1_center_of_mass(0,0,(T)1),end2_center_of_mass(0,0,(T)-1);
            T magnitude=attachment_twist.angular.Magnitude();
            if(magnitude>1e-6){
                orientation=ROTATION<TV>(time*magnitude,attachment_twist.angular);
                opposite_orientation=ROTATION<TV>(-time*magnitude,attachment_twist.angular);}
            X(index)=end1_center_of_mass+orientation.Rotate(X_save(index)-end1_center_of_mass)+time*attachment_twist.linear;
            X(X.Size()+index)=end2_center_of_mass+opposite_orientation.Rotate(X_save(X.Size()+index)-end2_center_of_mass)-time*attachment_twist.linear;}}
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    assert(test_number==6 || test_number==7 || test_number==8 || test_number==9);
    for(int j=0,index=0;j<mattress_grid.counts.y;j++) for(int k=0;k<mattress_grid.counts.x;k++){index++;
        X(index)=TV();X(X.Size()+index)=TV();}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==5 && id==int(1)){
        if(time<2) twist.linear=TV();
        else if(time<(T)3.5) twist.linear=TV(1,(T).5,0);
        else if(time<4) twist.linear=TV(0,(T)-1.5,0);
        else twist.linear=TV(-1.5,0,0);}
    else if(test_number==10 && id==int(0)){twist.linear=test_10_frame_velocity;}
    else if(test_number==10 && id==int(1)){twist.linear=-TV((T).88,0,0)+test_10_frame_velocity;}
    else if(test_number==11 && id==int(2)){if(time>=(T)1.6666667) twist.linear=TV((T).88,0,0);}
    else if(test_number==12 && id==int(5)){if(time<(T)1.1) twist.linear=-TV(0,(T)1.0,0);else twist.linear=TV();}
    else if((test_number==13 || test_number==19) && id==int(1)){if(time>(T).75){twist.angular=TV(0,(T)-pi/4,0);}
        //if(time>(T)8.333333){twist.linear=TV((T).5,(T).5,0);}
    }
    else if(test_number==27 && id==int(0)){
        INTERPOLATION_CURVE<T,TV> curve;
        curve.Add_Control_Point((T)0,TV(0,(T)1.5,0));
        curve.Add_Control_Point((T)1.2,TV(0,(T)1.5,(T)1.8));
        curve.Add_Control_Point((T)4,TV(0,(T)1.5,(T)1.8));
        curve.Add_Control_Point((T)5.2,TV(0,(T)0,-1));
        twist.linear=curve.Derivative(time);}
    else if(test_number==30 && id==int(1)){
        //INTERPOLATION_CURVE<T,TV> curve;
        // if(time>test_30_friction_off) rigid_body_2.coefficient_of_friction=(T).3; // TODO: This no longer makes sense here
        //RIGID_BODY<TV>& rigid_body_3=solid_body_collection.rigid_body_collection.Rigid_Body(2);
        //curve.Add_Control_Point((T)4.16,TV(-.5,(T)1.5,(T).7));
        //curve.Add_Control_Point((T)6,TV(-.5,(T)0,(T)-4));
        //twist.linear=curve.Derivative(time);
    }
    else if(test_number==32 && id==int(0)){twist.angular=TV(0,(T)-pi/4,0);}
    else if(test_number==46 && id==int(0)){twist=TWIST<TV>();}
    else return false;
    return true;
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE
{
    if(test_number==5 && id==int(1)){
        if(time<2) frame.t=TV(0,0,(T).5);
        else if(time<(T)3.5) frame.t=TV((time-2),(T).5*(time-2),(T).5);
        else if(time<4) frame.t=TV((T)1.5,(T)((T).75-(T)1.5*(time-(T)3.5)),(T).5);
        else frame.t=TV((T)(1.5-1.5*(time-4)),0,(T).5);}
    else if(test_number==10 && id==int(0)){frame.t=TV((T)0,-(T).25,(T).5)+time*test_10_frame_velocity;}
    else if(test_number==10 && id==int(1)){frame.t=TV((T)2,-(T).25,(T).5)+time*(-TV((T).88,0,0)+test_10_frame_velocity);}
    else if(test_number==11 && id==int(2)){if(time>=(T)1.6666667) frame.t=TV((T)1.2,(T).25,(T).5)+(time-(T)1.666667)*TV((T).88,0,0);}
    else if(test_number==12 && id==int(5)){if(time<(T)1.1) frame.t=TV(0,(T)1.25,0)+time*-TV(0,(T)1.0,0);}
    else if((test_number==13 || test_number==19) && id==int(1)){
        frame.t=TV(0,(T).30,0);
        if(time>(T).75) frame.r=ROTATION<TV>::From_Rotation_Vector(time*TV(0,(T)-pi/4,0));
        //if(time>(T)8.333333){frame.t=TV(0,(T).25,0)+(time-(T)8.333333)*TV((T).5,(T).5,0);}
    }
    else if(test_number==27 && id==int(0)){
        INTERPOLATION_CURVE<T,TV> curve;
        curve.Add_Control_Point((T)0,TV(0,(T)1.5,0));
        curve.Add_Control_Point((T)1.2,TV(0,(T)1.5,(T)1.8));
        curve.Add_Control_Point((T)4,TV(0,(T)1.5,(T)1.8));
        curve.Add_Control_Point((T)5.2,TV(0,(T)0,-1));
        frame.t=curve.Value(time);}
    else if(test_number==30 && id==int(1)){
        //INTERPOLATION_CURVE<T,TV> curve;
        //if(time>test_30_friction_off) rigid_body_2.coefficient_of_friction=(T).3;
        //RIGID_BODY<TV>& rigid_body_3=solid_body_collection.rigid_body_collection.Rigid_Body(2);
        //curve.Add_Control_Point((T)4.16,TV(-.5,(T)1.5,(T).7));
        //curve.Add_Control_Point((T)6,TV(-.5,(T)0,(T)-4));
        //frame.t=curve.Value(time);
    }
    else if(test_number==32 && id==int(0)){frame.r=ROTATION<TV>::From_Rotation_Vector(time*TV(0,(T)-pi/4,0))*ROTATION<TV>((T)pi/2,TV(1,0,0));}
}
//#####################################################################
// Function Bind_Redgreen_Segment_Midpoints
//#####################################################################
void Bind_Redgreen_Segment_Midpoints(RED_GREEN_TRIANGLES<TV>& redgreen)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    redgreen.Initialize_Segment_Index_From_Midpoint_Index();
    ARRAY<int> parents;ARRAY<T> weights;
    for(int s=0;s<redgreen.segment_midpoints.m;s++) if(const int midpoint=redgreen.segment_midpoints(s)){
        redgreen.Unrefined_Parents(midpoint,parents,weights);
        switch(parents.m){
            case 2:binding_list.Add_Binding(new LINEAR_BINDING<TV,2>(particles,midpoint,VECTOR<int,2>(parents),VECTOR<T,2>(weights)));break;
            case 3:binding_list.Add_Binding(new LINEAR_BINDING<TV,3>(particles,midpoint,VECTOR<int,3>(parents),VECTOR<T,3>(weights)));break;
            default: PHYSBAM_FATAL_ERROR("unexpected number of parents in refined triangle mesh");}}
}
//#####################################################################
};
}
#endif
