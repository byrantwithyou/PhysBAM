//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS
//#####################################################################
//    1. Sphere free falling to the ground
//    2. Torus free falling to the ground
//    3. Maggot free falling to the ground
//    4. Large armadillo caving in on itself
//    5. Deformable ball falling on a rigid ball
//    6. Smash test - small boxes, large sphere
//    7. Plane test - ball
//    8. Falling mattress
//    9. Testing rigid objects stuff
//   10. Increasing gravity
//   11. Increasing gravity (individual)
//   16. Smash test - large boxes, small mattress
//   17. Matress, no gravity, random start
//   18. Matress, no gravity, point start
//   23. Big lateral stretch
//   24. Big 6 sides stretch
//   25. Big 8 corners stretch
//   26. Big stretch/bend
//   27. Force inversion
//   28. Taffy test
//   29. Armadillo collapsing and rebounding
//   30. Projectile hitting a wall
//   31. Impact with sphere
//   32  Twisting chain
//   33. Through gears
//   34. Roll-over
//   35. Dancing jello (first attempt)
//   36. Spinning jello
//   37. Two jellos falling and colliding (first high res video)
//   38. Bunch of jellos falling
//   39. One jello falling on a stationary one
//   40. Two jellos running on each other
//   41. Bunch of jellos rolling towards camera
//   42. Bunch of jellos rolling towards camera (Reloaded)
//   43. Through smooth gears
//   44. 2 jellos collision
//   47. Fish past a magnet?
//   49. Hand through a tube?
//   50. Fish through a torus
//   51. Fish through a tube
//   52  Jello's falling one by one on each other
//   53. Stretch tet
//   54. Stretch tet (II)
//   55. Constrained tet
//   56. Size comparison - several shapes
//   57. Two-direction stretch
//   58. Various objects through gears
//   59. New taffy test
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <PhysBAM_Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/SVK_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
namespace PhysBAM{

template<class T_input>
class STANDARD_TESTS:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<VECTOR<T_input,3> > >
{
    typedef T_input T;
    typedef VECTOR<T,3> TV;typedef VECTOR<int,3> TV_INT;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<GRID<TV> > BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::stream_type;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::data_directory;

    std::ofstream svout;
    SOLIDS_STANDARD_TESTS<TV> tests;

    GRID<TV> mattress_grid,mattress_grid2,mattress_grid3,mattress_grid1;
    T attachment_velocity;
    bool semi_implicit;
    bool test_forces;
    bool use_extended_neohookean;
    bool use_extended_neohookean2;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_extended_neohookean_smooth;
    bool use_extended_svk, use_svk;
    bool use_corotated;
    bool use_mooney_rivlin,use_extended_mooney_rivlin;
    bool use_corot_blend;
    bool dump_sv;
    bool with_bunny,with_hand,with_big_arm,gears_of_pain;
    bool override_collisions,override_no_collisions;
    int kinematic_id,kinematic_id2,kinematic_id3,kinematic_id4,kinematic_id5;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve,curve2,curve3,curve4,curve5;
    bool print_matrix;
    int parameter,degrees_incline;
    int fishes,jello_size,number_of_jellos;
    T stiffness_multiplier;
    T damping_multiplier;
    T boxsize;
    T rebound_time,rebound_stiffness;
    bool use_constant_ife;
    bool forces_are_removed;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    ARRAY<TV> jello_centers;
    T stretch;
    T hole;
    bool nobind;
    ARRAY<TV> fish_V;
    T input_cutoff;
    T input_efc;
    T input_poissons_ratio,input_youngs_modulus;
    T input_friction;
    

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),use_extended_neohookean2(false),
        use_extended_neohookean_refined(false),use_extended_neohookean_hyperbola(false),use_extended_neohookean_smooth(false),use_extended_svk(false),
        use_corotated(false),use_corot_blend(false),dump_sv(false),print_matrix(false),use_constant_ife(false),input_cutoff(0),input_efc(0),input_poissons_ratio(-1),input_youngs_modulus(0)
    {
    }

    virtual ~STANDARD_TESTS()
    {}

    // Unused callbacks
    void Post_Initialization() PHYSBAM_OVERRIDE {}
    void Postprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Apply_Constraints(const T dt,const T time) PHYSBAM_OVERRIDE {}
//    void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt) PHYSBAM_OVERRIDE {}
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    //void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    //void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    // void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    // void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
  //  bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
   // void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    parse_args->Add_Option_Argument("-semi_implicit","use semi implicit forces");
    parse_args->Add_Option_Argument("-test_forces","use fully implicit forces");
    parse_args->Add_Option_Argument("-use_svk");
    parse_args->Add_Option_Argument("-use_ext_neo");
    parse_args->Add_Option_Argument("-use_ext_neo2");
    parse_args->Add_Option_Argument("-use_ext_neo_ref");
    parse_args->Add_Option_Argument("-use_ext_neo_hyper");
    parse_args->Add_Option_Argument("-use_ext_neo_smooth");
    parse_args->Add_Option_Argument("-use_ext_svk");
    parse_args->Add_Option_Argument("-use_ext_mooney");
    parse_args->Add_Option_Argument("-use_mooney");
    parse_args->Add_Option_Argument("-use_corotated");
    parse_args->Add_Option_Argument("-use_corot_blend");
    parse_args->Add_Option_Argument("-with_bunny");
    parse_args->Add_Option_Argument("-with_hand");
    parse_args->Add_Option_Argument("-with_big_arm");
    parse_args->Add_Option_Argument("-dump_sv");
    parse_args->Add_Integer_Argument("-parameter",0,"parameter used by multiple tests to change the parameters of the test");
    parse_args->Add_Double_Argument("-stiffen",1,"","stiffness multiplier for various tests");
    parse_args->Add_Double_Argument("-dampen",1,"","damping multiplier for various tests");
    parse_args->Add_Option_Argument("-residuals","print residuals during timestepping");
    parse_args->Add_Option_Argument("-print_energy","print energy statistics");
    parse_args->Add_Double_Argument("-cgsolids",1e-3,"CG tolerance for backward Euler");
    parse_args->Add_Option_Argument("-use_be","use backward euler");
    parse_args->Add_Option_Argument("-print_matrix");
    parse_args->Add_Option_Argument("-project_nullspace","project out nullspace");
    parse_args->Add_Integer_Argument("-projection_iterations",5,"number of iterations used for projection in cg");
    parse_args->Add_Integer_Argument("-solver_iterations",1000,"number of iterations used for solids system");
    parse_args->Add_Option_Argument("-use_constant_ife","use constant extrapolation on inverting finite element fix");
    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-collisions","Does not yet work in all sims, see code for details");
    parse_args->Add_Option_Argument("-no_collisions","Does not yet work in all sims, see code for details");
    parse_args->Add_Double_Argument("-stretch",1,"stretch");
    parse_args->Add_Double_Argument("-hole",.5,"hole");
    parse_args->Add_Double_Argument("-rebound_time",.2,"number of seconds to rebound in test 29");
    parse_args->Add_Double_Argument("-rebound_stiffness",5,"log10 of youngs modulus of final stiffness");
    parse_args->Add_Option_Argument("-nobind");
    parse_args->Add_Double_Argument("-cutoff",.4,"cutoff");
    parse_args->Add_Double_Argument("-efc",20,"efc");
    parse_args->Add_Double_Argument("-poissons_ratio",-1,"poissons_ratio");
    parse_args->Add_Double_Argument("-youngs_modulus",0,"youngs modulus, only for test 41 so far");
    parse_args->Add_Integer_Argument("-jello_size",20,"resolution of each jello cube");
    parse_args->Add_Integer_Argument("-number_of_jellos",12,"number of falling jello cubes in test 41");
    parse_args->Add_Integer_Argument("-degrees_incline",5,"degrees of incline");
    parse_args->Add_Double_Argument("-friction",.3,"amount of friction");
    parse_args->Add_Option_Argument("-gears_of_pain");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d",test_number);
    frame_rate=24;
    parameter=parse_args->Get_Integer_Value("-parameter");
    jello_size=parse_args->Get_Integer_Value("-jello_size");
    
    switch(test_number){
        case 17: case 18: case 24: case 25: case 27: case 10: case 11: case 23: case 57:
            if(!parameter) parameter=10;
            mattress_grid=GRID<TV>(parameter+1,parameter+1,parameter+1,(T)-1,(T)1,(T)-1,(T)1,(T)-1,(T)1);
            break;
        case 34:
            mattress_grid=GRID<TV>(13,13,13,(T)-2,(T)2,(T)-2,(T)2,(T)-2,(T)2);
            break;
        case 26:
            mattress_grid=GRID<TV>(40,5,5,(T)-4,(T)4,(T)-.5,(T).5,(T)-.5,(T).5);
            break;
        case 28:
            mattress_grid=GRID<TV>(80,10,10,(T)-8,(T)8,(T)-.5,(T).5,(T)-.5,(T).5);
            break;
        case 16:
            mattress_grid=GRID<TV>(11,6,11,(T)-1,(T)1,(T)-.5,(T).5,(T)-1,(T)1);
            break;
        case 35: case 36: case 41:
            mattress_grid=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            mattress_grid1=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            mattress_grid2=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.012,(T)0.012,(T)-0.012,(T)0.012,(T)-0.012,(T)0.012);
            mattress_grid3=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.015,(T)0.015,(T)-0.015,(T)0.015,(T)-0.015,(T)0.015);
            break;
        case 37: case 39: case 40: case 38: case 44:
            mattress_grid=GRID<TV>(40,40,40,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            break;
        case 42: case 52:
            mattress_grid=GRID<TV>(20,20,20,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            break;
    	default:
            mattress_grid=GRID<TV>(20,10,20,(T)-1,(T)1,(T)-.5,(T).5,(T)-1,(T)1);
    }

    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    solids_parameters.use_trapezoidal_rule_for_velocities=!parse_args->Get_Option_Value("-use_be");
    solids_parameters.use_rigid_deformable_contact=false;
    solid_body_collection.deformable_body_collection.soft_bindings.use_gauss_seidel_for_impulse_based_collisions=true;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.triangle_collision_parameters.use_gauss_jacobi=true;
    solids_parameters.triangle_collision_parameters.repulsions_limiter_fraction=1;
    solids_parameters.triangle_collision_parameters.collisions_final_repulsion_limiter_fraction=.1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    stretch=(T)parse_args->Get_Double_Value("-stretch");
    input_friction=(T)parse_args->Get_Double_Value("-friction");
    test_forces=parse_args->Is_Value_Set("-test_forces");
    use_extended_neohookean=parse_args->Is_Value_Set("-use_ext_neo");
    use_extended_neohookean2=parse_args->Is_Value_Set("-use_ext_neo2");
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_ref");
    use_extended_neohookean_hyperbola=parse_args->Is_Value_Set("-use_ext_neo_hyper");
    use_extended_mooney_rivlin=parse_args->Is_Value_Set("-use_ext_mooney");
    use_mooney_rivlin=parse_args->Is_Value_Set("-use_mooney");
    use_extended_neohookean_smooth=parse_args->Is_Value_Set("-use_ext_neo_smooth");
    use_svk=parse_args->Is_Value_Set("-use_svk");
    use_extended_svk=parse_args->Is_Value_Set("-use_ext_svk");
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    use_constant_ife=parse_args->Get_Option_Value("-use_constant_ife");
    solids_parameters.implicit_solve_parameters.test_system=parse_args->Is_Value_Set("-test_system");
    override_collisions=parse_args->Is_Value_Set("-collisions");
    override_no_collisions=parse_args->Is_Value_Set("-no_collisions")&&(!override_collisions);
    hole=(T)parse_args->Get_Double_Value("-hole");
    rebound_stiffness=(T)parse_args->Get_Double_Value("-rebound_stiffness");
    rebound_time=(T)parse_args->Get_Double_Value("-rebound_time");
    with_bunny=parse_args->Is_Value_Set("-with_bunny");
    with_hand=parse_args->Is_Value_Set("-with_hand");
    with_big_arm=parse_args->Is_Value_Set("-with_big_arm");
    number_of_jellos=parse_args->Get_Integer_Value("-number_of_jellos");
    degrees_incline=parse_args->Get_Integer_Value("-degrees_incline");
    gears_of_pain=parse_args->Is_Value_Set("-gears_of_pain");
    
    semi_implicit=parse_args->Is_Value_Set("-semi_implicit");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
    nobind=parse_args->Is_Value_Set("-nobind");
    if(parse_args->Is_Value_Set("-cutoff")) input_cutoff=(T)parse_args->Get_Double_Value("-cutoff");
    if(parse_args->Is_Value_Set("-efc")) input_efc=(T)parse_args->Get_Double_Value("-efc");
    if(parse_args->Is_Value_Set("-poissons_ratio")) input_poissons_ratio=(T)parse_args->Get_Double_Value("-poissons_ratio");
    if(parse_args->Is_Value_Set("-youngs_modulus")) input_youngs_modulus=(T)parse_args->Get_Double_Value("-youngs_modulus");
    
    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 7:
        case 8:
        case 9:
        case 10:
        case 11:
        case 16:
        case 17:
        case 18:
        case 56:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            break;
        case 5:
        case 6:
            frame_rate=24;
            last_frame=(int)(3*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            break;
        case 24:
        case 25:
        case 26:
        case 27: case 23: case 53: case 54: case 55: case 57:
            attachment_velocity = 0.2;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=2000;
            break;
        case 28:
            attachment_velocity = 0.4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            //solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=1000;
        case 29:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            //if (with_hand || with_bunny)
           // {
                //solids_parameters.triangle_collision_parameters.perform_self_collision=true;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
            //}
            frame_rate=120;
            last_frame=600;
            break;
        case 30:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            frame_rate=150;
            break;
        case 31:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            last_frame=5000;
            frame_rate=240;
            break;
        case 32:
            // solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            last_frame = 4000;
            break;
        case 33:
        case 43:
        case 58:
            solids_parameters.cfl=(T)10;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            frame_rate=120;
            last_frame=10*120;
            break;
        case 34:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            last_frame = 400;
            break;
        case 35: case 36: //case 41:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            // solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            // solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            last_frame = 1000;
            break;
        case 37:
        case 38:
        case 39:
        case 40:
        case 42:
        case 44:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            if (override_no_collisions) solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            frame_rate=960;
            last_frame=1000;
            break;
        case 41:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            if (override_no_collisions) solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            frame_rate=240;
            last_frame=2000;
            break;
        case 52:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            if (override_no_collisions) solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            frame_rate=360;
            last_frame=4000;
            break;
        case 47:
            frame_rate=24;
            last_frame=480;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            break;
        case 49:
            frame_rate=24;
            last_frame=480;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            break;
        case 48:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;//(int)(200*frame_rate);
            solids_parameters.cfl=(T)1;
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            break;
        case 50:
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            break;
        case 51:
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            break;
        default:
            LOG::cerr<<"Parsing: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    solid_body_collection.Print_Residuals(parse_args->Get_Option_Value("-residuals"));
    solid_body_collection.print_energy=parse_args->Get_Option_Value("-print_energy");
}
void Parse_Late_Options() PHYSBAM_OVERRIDE {BASE::Parse_Late_Options();}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
void Get_Initial_Data()
{
    bool automatically_add_to_collision_structures=true;
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    T density=TV::dimension==1?1:TV::dimension==2?100:1000;
    switch(test_number){
        case 1: case 7:{
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
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 4: case 29:{
            //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0))),true,true,density,.005);
            
            if (with_hand)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0),ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,1.0);
            else if(with_bunny)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0))),true,true,density,1.0);
            else if(with_big_arm)
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_380K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005);
            else
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T).4,(T)0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005);
            //            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)2.3,0))),true,true,density,5.0);
//            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
           // tests.Create_Mattress(GRID<TV>(TV_INT(13,13,13),RANGE<TV>(TV(-1,1,-1),TV(1,3,1))),true,0);
            
            tests.Add_Ground();

            break;}
        case 5:{
            RIGID_BODY<TV>& tmp_sphere=tests.Add_Rigid_Body("sphere",(T)1.0,(T).5);
            //RIGID_BODY<TV>& tmp_sphere=tests.Add_Analytic_Box(TV(1,1,1));
            tmp_sphere.X()=TV(0,(T).25,0);
            tmp_sphere.is_static=true;
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(.3,(T)6,0))),true,true,density,.5);
            tests.Add_Ground();
            break;}
        case 6:{
            //RIGID_BODY<TV>& tmp_sphere=tests.Add_Rigid_Body("sphere",(T)1.0,(T).5);
            RIGID_BODY<TV>& bottom_box=tests.Add_Analytic_Box(TV(1,1,1));
            RIGID_BODY<TV>& top_box=tests.Add_Analytic_Box(TV(1,1,1));
            bottom_box.X()=TV(0,(T)0,0);
            top_box.X()=TV(0,(T)2,0);
            bottom_box.is_static=true;
            top_box.is_static=false;
            kinematic_id=top_box.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(top_box.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,(T)2,0)));
            curve.Add_Control_Point(1,FRAME<TV>(TV(0,(T)1,0)));
            curve.Add_Control_Point(2,FRAME<TV>(TV(0,(T)1,0)));
            curve.Add_Control_Point(3,FRAME<TV>(TV(0,(T)2,0)));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)1,0))),true,true,density,.4);
            tests.Add_Ground();
            break;}
        case 8:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            tests.Add_Ground();
            break;}
        case 9:{
            RIGID_BODY<TV>& box1=tests.Add_Rigid_Body("cylinder",(T)1.0,(T).5);
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Cylinder(10,1);
            box1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(1,0,0));//ROTATION<TV>::From_Euler_Angles(-(T)pi/10,7*(T)pi/16,0);
            box1.X()=TV(0,-6,0);
            box2.X()=TV(0,6,0);
            box1.is_static=false;
            box2.is_static=false;
            break;
        }
        case 10:{
            last_frame=1200;
            for(int i=1;i<=7;i++){
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(4*i,1,0)));
                tests.Create_Mattress(mattress_grid,true,&initial_state);}
            tests.Add_Ground();
            break;}
        case 11:{
            last_frame=1200;
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,1,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            tests.Add_Ground();
            break;}
        case 16: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(20,20,20));
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(20,20,20));
            box1.X()=TV(0,-11,0);
            box2.X()=TV(0,11,0);
            box1.is_static=true;
            box2.is_static=false;
            kinematic_id=box2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,11,0)));
            curve.Add_Control_Point(5,FRAME<TV>(TV(0,8,0)));
            curve.Add_Control_Point(6,FRAME<TV>(TV(0,8,0)));
            curve.Add_Control_Point(11,FRAME<TV>(TV(0,11,0)));
            last_frame=250;
            break;}
        case 17:
        case 18:
        case 24:
        case 25:
        case 26: case 23: case 57:
        case 27:{
            tests.Create_Mattress(mattress_grid,true,0);
            break;}
        case 28: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            RIGID_BODY<TV>& cylinder1=tests.Add_Analytic_Cylinder(10,1);
            cylinder1.X()=TV(0,6,0);
            cylinder1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(1,0,0));
            RIGID_BODY<TV>& cylinder2=tests.Add_Analytic_Cylinder(10,1);
            cylinder2.X()=TV(0,-6,0);
            cylinder2.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,0,1));
            cylinder1.is_static=false;
            cylinder2.is_static=false;
            kinematic_id= cylinder1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder1.particle_index)=true;
            kinematic_id2=cylinder2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder2.particle_index)=true;
            T period = 10.0; T start = 0.0; T radius = 6.0; int pts = 75; T freq = .25;
            for (int ind=0; ind <= pts; ind++){
                curve.Add_Control_Point(start+freq*ind,FRAME<TV>(TV(radius*sin(2.0*pi*freq*ind/period),radius*cos(2.0*pi*freq*ind/period),0)));
                curve2.Add_Control_Point(start+freq*ind,FRAME<TV>(TV(-radius*sin(2.0*pi*freq*ind/period),-radius*cos(2.0*pi*freq*ind/period),0)));
            }
            last_frame=1000;
            break;}
        case 30: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)-90,(T)6,(T)11))),true,true,density, 1.0);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(20,20,20));
            box1.X()=TV(19,10,10);
            box1.is_static=true;
            tests.Add_Ground(1.0);
            last_frame=480;
            break;}
        case 31: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0))),true,true,density,.005);
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(0.125,density,5);
            sphere.is_static=false;
            tests.Add_Ground();
            kinematic_id=sphere.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(sphere.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,0.5,-10)));
            curve.Add_Control_Point(20,FRAME<TV>(TV(0,0.5,90)));

            RIGID_BODY<TV>& ball1=tests.Add_Analytic_Sphere(0.065,density,6);
            ball1.X() = TV(0.27,0.645,-0.12);
            RIGID_BODY<TV>& ball2=tests.Add_Analytic_Sphere(0.065,density,6);
            ball2.X() = TV(-0.275,0.605,-0.18);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(0.4,0.06,0.4));
            box1.X() = TV(-0.024,0.03,0.1);
            RIGID_BODY<TV>& cylinder1=tests.Add_Analytic_Cylinder(0.4,0.06,196);
            cylinder1.X() = TV(-0.024,0,0.1)+TV(0.2,0,0);
            RIGID_BODY<TV>& cylinder2=tests.Add_Analytic_Cylinder(0.4,0.06,196);
            cylinder2.X() = TV(-0.024,0,0.1)+TV(-0.2,0,0);
            RIGID_BODY<TV>& cylinder3=tests.Add_Analytic_Cylinder(0.4,0.06,196);
            cylinder3.X() = TV(-0.024,0,0.1)+TV(0,0,0.2);
            cylinder3.Rotation() = ROTATION<TV>((T)pi/2,TV(0,1,0));
            RIGID_BODY<TV>& cylinder4=tests.Add_Analytic_Cylinder(0.4,0.06,196);
            cylinder4.X() = TV(-0.024,0,0.1)+TV(0,0,-0.2);
            cylinder4.Rotation() = ROTATION<TV>((T)pi/2,TV(0,1,0));
            RIGID_BODY<TV>& sphere1=tests.Add_Analytic_Sphere(0.06,density,6);
            sphere1.X() = TV(-0.024,0,0.1)+TV(0.2,0,0.2);
            RIGID_BODY<TV>& sphere2=tests.Add_Analytic_Sphere(0.06,density,6);
            sphere2.X() = TV(-0.024,0,0.1)+TV(0.2,0,-0.2);
            RIGID_BODY<TV>& sphere3=tests.Add_Analytic_Sphere(0.06,density,6);
            sphere3.X() = TV(-0.024,0,0.1)+TV(-0.2,0,0.2);
            RIGID_BODY<TV>& sphere4=tests.Add_Analytic_Sphere(0.06,density,6);
            sphere4.X() = TV(-0.024,0,0.1)+TV(-0.2,0,-0.2);

            ball1.is_static = true;
            ball2.is_static = true;
            box1.is_static = true;
            cylinder1.is_static = true;
            cylinder2.is_static = true;
            cylinder3.is_static = true;
            cylinder4.is_static = true;
            sphere1.is_static = true;
            sphere2.is_static = true;
            sphere3.is_static = true;
            sphere4.is_static = true;
            
            break;}

        case 32:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-2.8,0,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(2.8,0,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-5.6,0,0))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(5.6,0,0))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-8.4,0,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(8.4,0,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);

            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((T).5,(T)2,32,64);
            torus1.is_static=false;
            torus1.coefficient_of_friction = 0.05;

            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((T).5,(T)2,32,64);
            torus2.is_static=false;
            torus2.coefficient_of_friction = 0.05;

            kinematic_id=torus1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(torus1.particle_index)=true;
            kinematic_id2=torus2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(torus2.particle_index)=true;

            T x_start = 11.2;
            T x_stop  = 15;
            T t_stop  = 2;
            T t_rot   = 8;

            curve.Add_Control_Point(0,FRAME<TV>(TV(x_start,0,0),ROTATION<TV>((T)0,TV(1,0,0))));
            curve.Add_Control_Point(t_stop,FRAME<TV>(TV(x_stop,0,0),ROTATION<TV>((T)0,TV(1,0,0))));

            for (int i=0; i<47; i++) curve.Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(x_stop,0,0),ROTATION<TV>((T)pi/2.0*i,TV(1,0,0))));

            curve2.Add_Control_Point(0,FRAME<TV>(TV(-x_start,0,0),ROTATION<TV>((T)0,TV(1,0,0))));
            curve2.Add_Control_Point(t_stop,FRAME<TV>(TV(-x_stop,0,0),ROTATION<TV>((T)0,TV(1,0,0))));

            for (int i=0; i<47; i++) curve2.Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(-x_stop,0,0),ROTATION<TV>((T)-pi/2.0*i,TV(1,0,0))));

            break;}
        case 33:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0.9,0),ROTATION<TV>((T)pi*0.525,TV(1,0,0))*ROTATION<TV>((T)pi/2,TV(0,1,0)))),true,true,density,0.06);

            T scale = 0.3;

            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);

            gear1.coefficient_of_friction = 0.1;
            gear2.coefficient_of_friction = 0.1;

            kinematic_id=gear1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear1.particle_index)=true;
            kinematic_id2=gear2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear2.particle_index)=true;

            T angular_velocity = 1;

            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-(T).4*scale,1.5*scale,-.75*scale),ROTATION<TV>(-i,TV(0,0,1))));
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV((T).4*scale,1.5*scale,-.75*scale),ROTATION<TV>(i,TV(0,0,1))));}

            tests.Add_Ground();
            break;}
            
        case 58:{
            T scale = 0.3;

            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3*scale,0),ROTATION<TV>(T(pi/2),TV(0,1,0))*ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,.35);
            if(!gears_of_pain){ tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4.5*scale,-1.2*scale),ROTATION<TV>((T)pi*0.525,TV(1,0,0))*ROTATION<TV>(0*(T)pi/2,TV(0,1,0)))),true,true,density,0.06);            
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.7*scale,-3.0*scale))),true,true,density,.25);
            }
           
            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            
            gear1.coefficient_of_friction = 0.1;
            gear2.coefficient_of_friction = 0.1;
            
            kinematic_id=gear1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear1.particle_index)=true;
            kinematic_id2=gear2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear2.particle_index)=true;
            
            T angular_velocity = 1; T gear_dx;
            if (gears_of_pain) gear_dx=(T).377; else gear_dx=(T).4;
            
            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,-.75*scale),ROTATION<TV>(-i,TV(0,0,1)))); //.4 is default, gears barely touch at .375
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,-.75*scale),ROTATION<TV>(i,TV(0,0,1))));}

            RIGID_BODY<TV>& box0=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale));
            
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale));
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale));
            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1.5*scale,.06*scale);
            box0.X()=TV(0,4.0*scale,-3.0*scale);

            if(!gears_of_pain){
                RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale));
                box1.X()=TV(0,1.2*scale,1.0*scale);
                //box1.Rotation()=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
                box1.is_static=true;
            }
            
            box2.X()=TV(0,1.2*scale,-1.0*scale);
            //box2.Rotation()=ROTATION<TV>(-(T)pi/4.0,TV(1,0,0));
            box3.X()=TV(0,3.0*scale,-1.4*scale);

            box3.Rotation()=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
            box0.is_static=false; //Will move later
            box2.is_static=true;
            box3.is_static=false;
            
            box3.coefficient_of_friction = .0;
            kinematic_id3=box3.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box3.particle_index)=true; 
            curve3.Add_Control_Point(0,FRAME<TV>(TV(0,3.0*scale,-1.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));
            curve3.Add_Control_Point(1,FRAME<TV>(TV(0,3.0*scale,-1.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));
            curve3.Add_Control_Point(2,FRAME<TV>(TV(0,3.0*scale,-3.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));

            box0.coefficient_of_friction = .05;
            kinematic_id5=box0.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box0.particle_index)=true; 
            curve5.Add_Control_Point(0,FRAME<TV>(TV(0,4.0*scale,-3.0*scale)));
            curve5.Add_Control_Point(1.0,FRAME<TV>(TV(0,4.0*scale,-3.0*scale),ROTATION<TV>(0*(T)pi/12.0,TV(1,0,0))));
            curve5.Add_Control_Point(1.1,FRAME<TV>(TV(0,4.0*scale,-3.0*scale),ROTATION<TV>((T)pi/12.0,TV(1,0,0))));
            
            cylinder.X()=TV(0,5.0*scale,0*scale);
            cylinder.is_static=false;
            kinematic_id4=cylinder.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder.particle_index)=true;
            curve4.Add_Control_Point(0,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curve4.Add_Control_Point(1,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curve4.Add_Control_Point(2,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curve4.Add_Control_Point(3,FRAME<TV>(TV(0,2.0*scale,0*scale)));

            
            tests.Add_Ground();
            break;}

        case 34:{
            T radius = 8.0;
            T velocity = 7;

            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder((T)32,radius,64);
            cylinder.is_static=false;
            cylinder.coefficient_of_friction = 1e8;
            kinematic_id=cylinder.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder.particle_index)=true;
            for (int i=0; i<128; i++) curve.Add_Control_Point(i,FRAME<TV>(TV(-25+i*velocity,radius*1.05,0),ROTATION<TV>(-i*velocity/radius,TV(0,0,1))));

            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,2,0)));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,3,0))),true,true,density,3);
            // tests.Create_Mattress(mattress_grid,true,&initial_state);
            tests.Add_Ground(1e8);
            break;}
        case 35:{
             RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(8,7,.5),ROTATION<TV>(T(pi/3),TV(1,0,1))));
             RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(8,1.2,0),ROTATION<TV>(T(pi/6),TV(0,1,0))));
            RIGID_BODY_STATE<TV> initial_state3(FRAME<TV>(TV(-.1,6,0),ROTATION<TV>(T(pi/4),TV(1,1,1))));
            RIGID_BODY_STATE<TV> initial_state4(FRAME<TV>(TV(0,1.5,0),ROTATION<TV>(T(0),TV(1,1,1))));

            tests.Create_Mattress(mattress_grid1,true,&initial_state1);
            tests.Create_Mattress(mattress_grid2,true,&initial_state2);
            tests.Create_Mattress(mattress_grid2,true,&initial_state3);
            tests.Create_Mattress(mattress_grid3,true,&initial_state4);
            tests.Add_Ground();
            break;
        }
        case 36:{
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(0,8,0),ROTATION<TV>(T(0),TV(1,0,1))));
            //RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(8,1.2,0),ROTATION<TV>(T(pi/6),TV(0,1,0))));
            //RIGID_BODY_STATE<TV> initial_state3(FRAME<TV>(TV(-.1,6,0),ROTATION<TV>(T(pi/4),TV(1,1,1))));
            //RIGID_BODY_STATE<TV> initial_state4(FRAME<TV>(TV(0,1.5,0),ROTATION<TV>(T(0),TV(1,1,1))));

            tests.Create_Mattress(mattress_grid1,true,&initial_state1);
            //tests.Create_Mattress(mattress_grid2,true,&initial_state2);
            //tests.Create_Mattress(mattress_grid2,true,&initial_state3);
            //tests.Create_Mattress(mattress_grid3,true,&initial_state4);
            tests.Add_Ground();
            break;
        }
        case 37:
        {
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(1.7,5.2,1.6),ROTATION<TV>(T(pi/4),TV(1.3,0.3,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(0,2.2,0),ROTATION<TV>(T(pi/7),TV(0.5,2,0.3))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground();
            break;
        }
        case 38:
        {
            number_of_jellos = 13;

            for (int i=1; i<=number_of_jellos; i++)
            {
                jello_centers.Append(TV(5.3*sin(277*i),15+2*cos(123*i),5.3*cos(297*i)));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(i),ROTATION<TV>(10*sin(178*i),TV(sin(145*i),cos(345*i),cos(478*i)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state);
            }
            tests.Add_Ground();
            break;
        }
        case 39:
        {
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(0.07,0.3,0.09),ROTATION<TV>(T(pi/0.103),TV(1.35,0.785,1.675))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(0,0.01,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground();
            break;
        }
        case 40:
        {
            jello_centers.Append(TV(-0.266,0.022,0.013)); 
            jello_centers.Append(TV(0.266, 0.029,-0.013));
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(1),ROTATION<TV>(T(pi/0.13),TV(1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(2),ROTATION<TV>(T(pi/0.076),TV(0.7,1,0.1))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground();
            break;
        }
        case 41:
      /*  {

            jello_centers.Append(TV(-0.054,0.042,0.013)); 
            jello_centers.Append(TV(0.031, 0.029,-0.013));
            jello_centers.Append(TV(-0.074,0.037,-0.043)); 
            jello_centers.Append(TV(0.052, 0.059,0.036));
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(1),ROTATION<TV>(-T(pi/0.11),TV(1.3,-1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(2),ROTATION<TV>(T(pi/0.076),TV(0.7,1,0.1))));
            RIGID_BODY_STATE<TV> initial_state3(FRAME<TV>(jello_centers(3),ROTATION<TV>(T(pi/0.13),TV(-1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state4(FRAME<TV>(jello_centers(4),ROTATION<TV>(-T(pi/0.046),TV(0.7,.1,-0.2))));
            tests.Create_Mattress(mattress_grid1,true,&initial_state1);
            tests.Create_Mattress(mattress_grid2,true,&initial_state2);
            tests.Create_Mattress(mattress_grid1,true,&initial_state3);
            tests.Create_Mattress(mattress_grid3,true,&initial_state4);
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(0.1);
            inclined_floor.Rotation()=ROTATION<TV>((T)pi/10,TV(1,0,0));
            break;
        }*/
        {
            //Todo: At some point, use Oriented_Box to do a tighter collision check. C-Re, C-e, R where R is the rotation matrix and e is the edge
            RANDOM_NUMBERS<T> random;
            T max_jello_size = .03;
            T bound = .3;
            TV new_center; T new_rotate;
            bool stuck=false;
            RIGID_BODY_STATE<TV> initial_state;
            
            random.Set_Seed(123);
            //break;} 
            for (int i=1; i<=number_of_jellos; i++){
                do {
                    random.Fill_Uniform(new_center,-bound,bound);
                   // new_center = TV(random.Get_Uniform_Number(-bound,bound),random.Get_Uniform_Number((T).5*bound,(T)1*bound),random.Get_Uniform_Number(-bound,bound));
                    stuck=false;
                    new_center.y = (T).5*(new_center.y + 1.5*bound);
                    for (int j=1; j<i&&(!stuck); j++){
                        if((new_center-jello_centers(j)).Magnitude()<=(T)4*max_jello_size*max_jello_size) stuck=true;
                       // if ((new_center.x-jello_centers(j).x)*(new_center.x-jello_centers(j).x)+(new_center.y-jello_centers(j).y)*(new_center.y-jello_centers(j).y)+(new_center.z-jello_centers(j).z)*(new_center.z-jello_centers(j).z)<=(T)4*max_jello_size*max_jello_size) stuck=true;
                }}while(stuck);
                jello_centers.Append(new_center);
                random.Fill_Uniform(new_center,-bound,bound);
//                     new_center = TV(random.Get_Uniform_Number(-bound,bound),random.Get_Uniform_Number(-bound,bound),random.Get_Uniform_Number(-bound,bound));
                new_rotate = random.Get_Uniform_Number(-(T)pi,(T)pi);
                RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(i),ROTATION<TV>(new_rotate,new_center)));
                if (i % 4 ==0) {tests.Create_Mattress(mattress_grid3,true,&initial_state1);}
                else if (i % 4 ==1) {tests.Create_Mattress(mattress_grid2,true,&initial_state1);}
                else {tests.Create_Mattress(mattress_grid1,true,&initial_state1);}
            }
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Rotation()=ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            break;
        }
        case 42:
        {
            int count = 0;
            for (int i=1; i<=3; i++)
            for (int j=1; j<=3; j++)
            for (int k=1; k<=3; k++)
            {
                count++;
                jello_centers.Append(TV(-100+i*5,j*5+3,k*5+sin(75*count)));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(count),ROTATION<TV>(10*sin(178*count),TV(sin(145*count),cos(345*count),cos(478*count)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state);
            }
            tests.Add_Ground();
            break;
        }
        case 52:
        {
            /*RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(0.07,0.035,0.08,128);
            shell.coefficient_of_friction = 0.3;
            shell.X()=TV(0,0.035,0);
            shell.Rotation()=ROTATION<TV>((T)pi/2.0,TV(1,0,0));
            shell.is_static=true;*/

            int count = 0;
            for (int i=1; i<=27; i++)
            {
                count++;
                jello_centers.Append(TV(-500+i*5,-1000,0));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(count),ROTATION<TV>(10*sin(178*i),TV(sin(145*i),cos(345*i),cos(478*i)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state);
            }
            
            int number_of_boxes = 128;
            T height = 0.08;
            T width = 0.02;
            T radius = 0.035;

            ARRAY<RIGID_BODY<TV>*> boxes;
            
            for (int i=1; i<=number_of_boxes; i++)
            {
                T phi = 2*pi*(i-1)/number_of_boxes;

                boxes.Append(&tests.Add_Analytic_Box(TV(width,height,width)));
                boxes(i)->X() = TV((radius+width/2)*sin(phi),height/2,(radius+width/2)*cos(phi));
                boxes(i)->Rotation() = ROTATION<TV>((T)phi,TV(0,1,0));
                boxes(i)->is_static = true;
                boxes(i)->coefficient_of_friction = 0.3;
            }

            /*RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box4=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box5=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box6=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box7=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));
            RIGID_BODY<TV>& box8=tests.Add_Analytic_Box(TV(0.07,0.07,0.07));

            box1.X()=TV(0.07,0.035,0);
            box2.X()=TV(-0.07,0.035,0);
            box3.X()=TV(0,0.035,0.07);
            box4.X()=TV(0,0.035,-0.07);
            
            box5.X()=TV(0.07/sqrt(2),0.035,0.07/sqrt(2));
            box6.X()=TV(-0.07/sqrt(2),0.035,0.07/sqrt(2));
            box7.X()=TV(-0.07/sqrt(2),0.035,-0.07/sqrt(2));
            box8.X()=TV(0.07/sqrt(2),0.035,-0.07/sqrt(2));

            box5.Rotation() = ROTATION<TV>((T)pi/4.0,TV(0,1,0));
            box6.Rotation() = ROTATION<TV>((T)pi/4.0,TV(0,1,0));
            box7.Rotation() = ROTATION<TV>((T)pi/4.0,TV(0,1,0));
            box8.Rotation() = ROTATION<TV>((T)pi/4.0,TV(0,1,0));

            box1.is_static=true;
            box2.is_static=true;
            box3.is_static=true;
            box4.is_static=true;
            box5.is_static=true;
            box6.is_static=true;
            box7.is_static=true;
            box8.is_static=true;*/

            tests.Add_Ground();
            break;
        }
        case 44:
        {
            jello_centers.Append(TV(-10,3.2,0));
            jello_centers.Append(TV(10,3.5,-1));
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(1),ROTATION<TV>(T(pi/4),TV(1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(2),ROTATION<TV>(T(pi/5),TV(0.7,1,0.1))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground();
            break;
        }
        case 50:
        {
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((T)1.4,(T)1.6,32,64);
            torus1.is_static=true;
            torus1.coefficient_of_friction = 0.00;
            torus1.X()=TV(10,10,0);
            torus1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            /*RIGID_BODY<TV>& torus2=tests.Add_Analytic_Shell((T).9,(T)1,(T).2,64);
            torus2.is_static=true;
            torus2.coefficient_of_friction = 0.05;
            torus2.X()=TV(0,3,0);
            torus2.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,0,1));*/


            ///tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_21K.tet",
               //                                 RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,0,0))),true,true,density);
            //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/torus_thin_24K.tet",
              //                                  RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-4.8,2,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);
            //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/SubmarineFluidMesh.tet",
              //                                  RIGID_BODY_STATE<TV>(FRAME<TV>(TV(4.8,4,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);
           // tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/snake_22K.tet",
             //                                   RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-9.6,0,0))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density);
            //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/shell-2-5K.tet",
              //                                  RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-14.4,0,0),ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density);

            tests.Add_Ground(1);
            break;
        }
        case 47:
        {
            fishes=5;
            for (int i=0; i<fishes; i++)
            {
              //  jello_centers.Append(TV(-500+i*5,-1000,0));
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                    RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-10*i,10,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density);

            }
            boxsize=(T)2;
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(3*boxsize,boxsize,boxsize));
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(3*boxsize,boxsize,boxsize));
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(boxsize,boxsize,boxsize));
            RIGID_BODY<TV>& box4=tests.Add_Analytic_Box(TV(boxsize,boxsize,boxsize));

            
            box1.X()=TV(0,.5*boxsize,boxsize);
            box2.X()=TV(0,.5*boxsize,-boxsize);
            box3.X()=TV(-boxsize,.5*boxsize,0);
            box4.X()=TV(boxsize,.5*boxsize,0);
            
            
            box1.is_static=true;
            box2.is_static=true;
            box3.is_static=true;
            box4.is_static=true;
            tests.Add_Ground();
            break;
        }
        case 49:
        {
            TV start(10,10,0);
            T outer=(T)1,inner=(T).5,length=10;
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus1.is_static=true;
            torus1.coefficient_of_friction = 0.05;
            torus1.X()=start;
            torus1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus2.is_static=true;
            torus2.coefficient_of_friction = 0.05;
            torus2.X()=start+TV(length,0,0);
            torus2.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(length,outer,inner,64);
            shell.is_static=true;
            shell.coefficient_of_friction = 0.05;
            shell.X()=start+TV(length/2,0,0);
            shell.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density);
            tests.Add_Ground(1);
            break;
        }
        case 51:
        {
            TV start(10,10,0);
            T outer=(T)2.5,inner=hole,length=10;
            RIGID_BODY<TV>& torus1=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus1.is_static=true;
            torus1.coefficient_of_friction = 0.05;
            torus1.X()=start;
            torus1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            torus2.is_static=true;
            torus2.coefficient_of_friction = 0.05;
            torus2.X()=start+TV(length,0,0);
            torus2.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(length,outer,inner,64);
            shell.is_static=true;
            shell.coefficient_of_friction = 0.05;
            shell.X()=start+TV(length/2,0,0);
            shell.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density);
            tests.Add_Ground(1);
            break;
        }
        case 43:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0),ROTATION<TV>(T(0),TV(0,1,0)))),true,true,density);
            int cogs=40;
            T radius=.8;
            T sep=radius+.003;
            T cog_rad=.035;
            RIGID_BODY<TV>& gear1=tests.Add_Analytic_Smooth_Gear(TV(radius,cog_rad,1),cogs,8);
            RIGID_BODY<TV>& gear2=tests.Add_Analytic_Smooth_Gear(TV(radius,cog_rad,1),cogs,8);
            gear1.coefficient_of_friction = 2;
            gear2.coefficient_of_friction = 2;
            kinematic_id=gear1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear1.particle_index)=true;
            kinematic_id2=gear2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear2.particle_index)=true;
            T angular_velocity = 1;

            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-sep,1.5,0),ROTATION<TV>(-i,TV(0,0,1))*ROTATION<TV>(pi/(2*cogs),TV(0,0,1))));
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(sep,1.5,0),ROTATION<TV>(i,TV(0,0,1))*ROTATION<TV>(pi/(2*cogs),TV(0,0,1))));}

            tests.Add_Ground();
            break;}
        case 53:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            int a=particles.array_collection->Add_Element();
            int b=particles.array_collection->Add_Element();
            int c=particles.array_collection->Add_Element();
            int d=particles.array_collection->Add_Element();
            particles.X(a)=TV(1,0,0);
            particles.X(b)=TV(0,-1,0);
            particles.X(c)=TV(0,.5,.5);
            particles.X(d)=TV(0,.5,-.5);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            particles.mass(d)=1;
            tv->mesh.elements.Append(VECTOR<int,4>(a,b,c,d));
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(tv);
            break;}
        case 54: case 55:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            int a=particles.array_collection->Add_Element();
            int b=particles.array_collection->Add_Element();
            int c=particles.array_collection->Add_Element();
            int d=particles.array_collection->Add_Element();
            particles.X(a)=TV(1,0,0);
            particles.X(b)=TV(0,0,0);
            particles.X(c)=TV(0,0,1);
            particles.X(d)=TV(0,1,0);
            particles.mass(a)=1;
            particles.mass(b)=1;
            particles.mass(c)=1;
            particles.mass(d)=1;
            tv->mesh.elements.Append(VECTOR<int,4>(a,b,c,d));
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(tv);
            break;}
        case 56:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)5,0))),true,true,density,.5);
            tests.Add_Ground();
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(5,5,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.2);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(10,5,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density,.005);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(15,5,0),ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,1.0);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(20,(T)5,0))),true,true,density,1.0);
            break;
        }
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

    // correct number nodes
    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++) deformable_body_collection.deformable_geometry.structures(i)->Update_Number_Nodes();

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
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    Get_Initial_Data();

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 8:
        case 16:
        case 5:
        case 6:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            break;}
        case 4: {
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01);
            break;}
        case 7:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            for(int i=1; i<=deformable_body_collection.particles.X.m; i++) deformable_body_collection.particles.X(i).y=3;
            break;}
        case 9: case 56:{break;}
        case 10:{
            bool* bools[7]={&use_corotated,0,&use_constant_ife,&use_corot_blend,&use_extended_neohookean,&use_extended_neohookean_smooth,&use_extended_neohookean_hyperbola};
            for(int i=1;i<=7;i++){
                if(bools[i-1]) *bools[i-1]=true;
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(i);
                Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
                if(bools[i-1]) *bools[i-1]=false;}
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 11:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)5e4,(T).45,(T).01);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 17:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,-1,1);
            break;}
        case 18:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            RANDOM_NUMBERS<T> rand;
            rand.Fill_Uniform(particles.X,-0,0);
            break;}
        case 23:
        case 24:
        case 25:
        case 26:
        case 27:
        case 57:
        case 28:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            if(test_number==57){
                int m=mattress_grid.counts.x;
                int n=mattress_grid.counts.y;
                int p=mattress_grid.counts.z;
                for(int i=1;i<=m;i++) for(int j=1;j<=n;j++){constrained_particles.Append(i+m*(j-1));constrained_particles.Append(i+m*(j-1)+(p-1)*m*n);}
                for(int i=1;i<=m;i++) for(int k=1;k<=p;k++){constrained_particles.Append(i+m*n*(k-1));constrained_particles.Append(i+m*(n-1)+m*n*(k-1));}
                constrained_velocities=particles.X.Subset(constrained_particles)*attachment_velocity;
                constrained_velocities.template Project<T,&TV::x>().Fill(0);}
            break;}
        case 29:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)0e2,(T).45,(T).01);
            forces_are_removed=true;
            break;}

        case 30:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).40,(T).01);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            for(int i=1; i<=deformable_body_collection.particles.X.m; i++){ deformable_body_collection.particles.V(i).x=(T)60;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
            break;}
        case 31:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4,(T).4,(T).001);

            // RIGID_BODY<TV>& box1=tests.Add_Analytic_Sphere(0.065,density,5);
            // box1.X() = TV(0.27,0.645,-0.12);
            // RIGID_BODY<TV>& box2=tests.Add_Analytic_Sphere(0.065,density,5);
            // box2.X() = TV(-0.275,0.605,-0.18);
            // RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(0.5,0.06,0.5));
            // box3.X() = TV(-0.024,0.03,0.1);

            for (int i=1; i<=particles.X.m; i++)
            {
                T y   = particles.X(i).y;
                TV v1 = particles.X(i) - TV(0.27,0.645,-0.12);
                TV v2 = particles.X(i) - TV(-0.275,0.605,-0.18);
                
                if ((y<=0.06) || (sqr(v1.x) + sqr(v1.y) + sqr(v1.z) <= sqr(0.065)) || (sqr(v2.x) + sqr(v2.y) + sqr(v2.z) <= sqr(0.065)))
                {
                    constrained_particles.Append(i);
                }
            }

            break;} 
        case 32:{
            T youngs_modulus = 7e5;
            T poissons_ratio = .4;
            T damping = 0.05;
            T g=0.8;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
            Add_Constitutive_Model(tetrahedralized_volume3,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume4=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(4);
            Add_Constitutive_Model(tetrahedralized_volume4,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume5=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(5);
            Add_Constitutive_Model(tetrahedralized_volume5,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume6=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(6);
            Add_Constitutive_Model(tetrahedralized_volume6,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume7=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(7);
            Add_Constitutive_Model(tetrahedralized_volume7,youngs_modulus,poissons_ratio,damping);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=g;
            break;}
        case 43:
        case 33:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,1e4,0.4,0.005);

            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 58:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,1e4,0.4,0.005);
            if(!gears_of_pain){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
                Add_Constitutive_Model(tetrahedralized_volume2,5e4,0.4,0.005);
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
                Add_Constitutive_Model(tetrahedralized_volume3,1e4,0.4,0.005);}
             
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 34:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .4;
            T damping = 0.1;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 35:{
            T youngs_modulus = 2e5;
            T poissons_ratio = .45;
            T damping = 0.01;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
            Add_Constitutive_Model(tetrahedralized_volume3,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume4=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(4);
            Add_Constitutive_Model(tetrahedralized_volume4,youngs_modulus,poissons_ratio,damping);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 36:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01);
            //solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            for(int i=1; i<=deformable_body_collection.particles.X.m/2; i++){ deformable_body_collection.particles.V(i).x=(T)1;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
            break;}
        case 37:{
            T youngs_modulus = 2.25e5;
            T poissons_ratio = .4;
            T damping = 0;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping,0.4,50);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping,0.4,50);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=1;i<=m;i++)
            for(int j=1;j<=n;j++)
            for(int ij=1;ij<=mn;ij++)
            {
                particles.V(i+m*(j-1)+m*n*(ij-1)) = TV(-2*sin(j/(T)n),-cos(2*ij/(T)mn)*2,-2*sin(3*i/(T)m));
                particles.V(i+m*(j-1)+m*n*(ij-1)+m*n*mn) = TV(1*sin(2*ij/(T)mn),0.5*cos(3*i/(T)m)*2,1*sin(j/(T)n));
            }
            break;}
        case 38:{
            T youngs_modulus = 4e5;
            T poissons_ratio = .4;
            T damping = 0.001;
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for (int k=1; k<=jello_centers.m; k++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping,0.4,50);
                for(int i=1;i<=m;i++)
                for(int j=1;j<=n;j++)
                for(int ij=1;ij<=mn;ij++)
                {
                    particles.V(i+m*(j-1)+m*n*(ij-1)+(k-1)*m*n*mn) = TV(-3*sin(sin(127*k)*j/(T)n),-5*cos(sin(384*k)*4*ij/(T)mn)*2+5*sin(385*k),-6*sin(3*i*sin(457*k)/(T)m));
                }
            }
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 39:{
            T youngs_modulus = 1e4;
            T poissons_ratio = .4;
            T damping = 0.001;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping,0.4,50);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping,0.4,50);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=1;i<=m;i++)
                for(int j=1;j<=n;j++)
                    for(int ij=1;ij<=mn;ij++)
                    {
                        particles.V(i+m*(j-1)+m*n*(ij-1)) = TV(-0.3*sin(j/(T)n)-0.15,-cos(2*ij/(T)mn)*0.3,-0.15*sin(3*i/(T)m)-0.3);
                    }
            break;}
        case 40:{
            T youngs_modulus = 1e4;
            T poissons_ratio = .4;
            T damping = 0.001;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping,0.4,50);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping,0.4,50);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=1;i<=m;i++)
                for(int j=1;j<=n;j++)
                    for(int ij=1;ij<=mn;ij++)
                    {
                        particles.V(i+m*(j-1)+m*n*(ij-1)) = TV(0.32*cos(6*ij/(T)mn)+2.15,-0.28*cos(4*i/(T)m)+0.55,-0.31*sin(4*j/(T)n));
                        particles.V(i+m*(j-1)+m*n*(ij-1)+m*n*mn) = TV(0.27*sin(6*j/(T)n)-2.17,0.32*cos(6*ij/(T)mn)+0.53,0.25*sin(4*i/(T)m));
                    }
            for (int i=1; i<=m*n*mn; i++)
            {
                int index = i;
                particles.V(index) += TV(particles.X(index).y-jello_centers(1).y,-(particles.X(index).x-jello_centers(1).x),0)*51;
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(1).z,-(particles.X(index).y-jello_centers(1).y))*(-13);
                index+=m*n*mn;
                particles.V(index) += TV(particles.X(index).y-jello_centers(2).y,-(particles.X(index).x-jello_centers(2).x),0)*(-39);
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(2).z,-(particles.X(index).y-jello_centers(2).y))*23;
            }
            break;}
        case 41:{
            T youngs_modulus = 1e3;
            T poissons_ratio = .4;
            T damping = 0.001;

            for (int k=1; k<=number_of_jellos; k++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping,0.4,50);
            }
            
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;
        }
        case 42:{
            T youngs_modulus = 3e5;
            T poissons_ratio = .4;
            T damping = 0.001;
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for (int k=1; k<=jello_centers.m; k++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping,0.4,50);
                for(int i=1;i<=m;i++)
                for(int j=1;j<=n;j++)
                for(int ij=1;ij<=mn;ij++)
                {
                    particles.V(i+m*(j-1)+m*n*(ij-1)+(k-1)*m*n*mn) = TV(-2*sin(sin(127*k)*j/(T)n)+17+4*sin(181*k),-3*cos(sin(384*k)*4*ij/(T)mn)*2+(2+4*sin(461*k)),1.5*sin(3*i*sin(457*k)/(T)m));
                }
                for (int i=1; i<=m*n*mn; i++)
                {
                    int index = i+(k-1)*m*n*mn;
                    particles.V(index) += TV(particles.X(index).y-jello_centers(k).y,-(particles.X(index).x-jello_centers(k).x),0)*(6+2*sin(453*k));
                    particles.V(index) += TV(0,particles.X(index).z-jello_centers(k).z,-(particles.X(index).y-jello_centers(k).y))*(cos(413*k));
                }
            }
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 52:{
            T youngs_modulus = 1e4;
            T poissons_ratio = .4;
            T damping = 0.001;
            for (int k=1; k<=jello_centers.m; k++)
            {
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(k);
                Add_Constitutive_Model(tetrahedralized_volume,youngs_modulus,poissons_ratio,damping,0.4,50);
            }
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            break;}
        case 44:{
            T youngs_modulus = 4e5;
            T poissons_ratio = .4;
            T damping = 0.001;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping,0.4,50);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping,0.4,50);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            int m=mattress_grid.counts.x;
            int n=mattress_grid.counts.y;
            int mn=mattress_grid.counts.z;
            for(int i=1;i<=m;i++)
                for(int j=1;j<=n;j++)
                    for(int ij=1;ij<=mn;ij++)
                    {
                        particles.V(i+m*(j-1)+m*n*(ij-1)) = TV(-2*sin(j/(T)n)+19,-cos(2*ij/(T)mn)*2+1,-2*sin(3*i/(T)m));
                        particles.V(i+m*(j-1)+m*n*(ij-1)+m*n*mn) = TV(1*sin(2*ij/(T)mn)-21,0.5*cos(3*i/(T)m)*2+1.2,1*sin(j/(T)n));
                    }
            for (int i=1; i<=m*n*mn; i++)
            {
                int index = i;
                particles.V(index) += TV(particles.X(index).y-jello_centers(1).y,-(particles.X(index).x-jello_centers(1).x),0)*9;
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(1).z,-(particles.X(index).y-jello_centers(1).y))*2;
                index+=m*n*mn;
                particles.V(index) += TV(particles.X(index).y-jello_centers(2).y,-(particles.X(index).x-jello_centers(2).x),0)*(-12);
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(2).z,-(particles.X(index).y-jello_centers(2).y))*(-1);
            }
            break;}
        case 50:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .45;
            T damping = 0.01;
           // T g=0.8;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            for(int i=1;i<=particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);
            /*for (int i=1; i<=m*n*mn; i++)
            {
                particles.V(i) += TV(particles.X(i).y-2.4,-(particles.X(i).x+30),0)*10;
            }*/
            /*TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
            Add_Constitutive_Model(tetrahedralized_volume2,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
            Add_Constitutive_Model(tetrahedralized_volume3,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume4=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(4);
            Add_Constitutive_Model(tetrahedralized_volume4,youngs_modulus,poissons_ratio,damping);
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume5=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(5);
            Add_Constitutive_Model(tetrahedralized_volume5,youngs_modulus,poissons_ratio,damping);*/
            /*TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume6=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(6);
            Add_Constitutive_Model(tetrahedralized_volume6,youngs_modulus,poissons_ratio,damping);*/
           // TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume7=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(7);
            //Add_Constitutive_Model(tetrahedralized_volume7,youngs_modulus,poissons_ratio,damping);
           // solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
           // solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=g;

            break;}
        case 51:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .45;
            T damping = 0.01;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            for(int i=1;i<=particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);
            break;}
        case 55:
        case 54:
        case 53:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .45;
            T damping = 0.01;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            if(test_number==55) particles.X(1).x=stretch;
            break;}
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        VECTOR<int,3> processes_per_dimension(2,1,1);
        deformable_body_collection.mpi_solids->Simple_Partition(deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,processes_per_dimension);}

    for(int i=1;i<=deformable_body_collection.deformable_geometry.structures.m;i++){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision)
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}

    if(parse_args->Is_Value_Set("-solver_iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-solver_iterations");
    if(parse_args->Is_Value_Set("-cgsolids")) solids_parameters.implicit_solve_parameters.cg_tolerance=(T)parse_args->Get_Double_Value("-cgsolids");

    if(!semi_implicit) for(int i=1;i<=solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=1;i<=solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=1;i<=deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}

//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    T final_time=50;
    if(test_number==24){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y = velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        TV velocity_z = velocity_time<final_time?TV(0,0,attachment_velocity):TV();
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int j=n/3+1;j<=2*n/3+1;j++){V(i+m*(j-1))=-velocity_z;V(i+m*(j-1)+(mn-1)*m*n)=velocity_z;}
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int ij=mn/3+1;ij<=2*mn/3+1;ij++){V(i+m*n*(ij-1))=-velocity_y;V(i+m*(n-1)+m*n*(ij-1))=velocity_y;}
        for(int ij=mn/3+1;ij<=2*mn/3+1;ij++)for(int j=n/3+1;j<=2*n/3+1;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
    }
    if(test_number==25){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        V(1)=         velocity_time<final_time?TV(-attachment_velocity,-attachment_velocity,-attachment_velocity):TV();
        V(m)=         velocity_time<final_time?TV( attachment_velocity,-attachment_velocity,-attachment_velocity):TV();
        V(m*(n-1)+1)= velocity_time<final_time?TV(-attachment_velocity, attachment_velocity,-attachment_velocity):TV();
        V(m*n)=       velocity_time<final_time?TV( attachment_velocity, attachment_velocity,-attachment_velocity):TV();
        V((mn-1)*n*m+1)=         velocity_time<final_time?TV(-attachment_velocity,-attachment_velocity,attachment_velocity):TV();
        V((mn-1)*n*m+m)=         velocity_time<final_time?TV( attachment_velocity,-attachment_velocity,attachment_velocity):TV();
        V((mn-1)*n*m+m*(n-1)+1)= velocity_time<final_time?TV(-attachment_velocity, attachment_velocity,attachment_velocity):TV();
        V(mn*n*m)=               velocity_time<final_time?TV( attachment_velocity, attachment_velocity,attachment_velocity):TV();
    }
    if(test_number==26){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y = velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
        for(int i=3*m/7+1;i<=4*m/7+1;i++)for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(i+m*(j-1)+m*n*(ij-1))=-velocity_y;V(i+m*(j-1)+m*n*(ij-1))=-velocity_y;}
    }
    if(test_number==27){
        final_time=70;
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=velocity_x;V(m+m*(j-1)+m*n*(ij-1))=-velocity_x;}
    }
    if(test_number==23){
        final_time=70;
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
    }
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y = velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
    }
    if(test_number==53){V(1)=TV(1,0,0);V(2).x=0;V(3).x=0;V(4).x=0;}
    if(test_number==54){V(1)=TV(1,0,0);V(2)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==55){V(2)=V(1)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==31)
    {
        int number_of_constrained_particles = constrained_particles.m;
        for (int i=1; i<=number_of_constrained_particles; i++)
            V(constrained_particles(i))=TV();
    }
    if(test_number==57) V.Subset(constrained_particles)=constrained_velocities;
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    if(test_number==24){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int j=n/3+1;j<=2*n/3+1;j++){V(i+m*(j-1))=TV();V(i+m*(j-1)+(mn-1)*m*n)=TV();}
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int ij=mn/3+1;ij<=2*mn/3+1;ij++){V(i+m*n*(ij-1))=TV();V(i+m*(n-1)+m*n*(ij-1))=TV();}
        for(int ij=mn/3+1;ij<=2*mn/3+1;ij++)for(int j=n/3+1;j<=2*n/3+1;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==25){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        V(1)=TV();
        V(m)=TV();
        V(m*(n-1)+1)=TV();
        V(m*n)=TV();
        V((mn-1)*n*m+1)=TV();
        V((mn-1)*n*m+m)=TV();
        V((mn-1)*n*m+m*(n-1)+1)=TV();
        V(mn*n*m)=TV();
    }
    if(test_number==26){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
        for(int i=3*m/7+1;i<=4*m/7+1;i++)for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(i+m*(j-1)+m*n*(ij-1))=TV();V(i+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==23){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==27){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
  /*  if(test_number==50){
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        int n=particles.array_collection->Size();
        for(int i=1; i <=n; i++)
        {
            if(externally_forced[i] &&V(i).x<=0){V(i)=TV();}
        }

    }*/
    if(test_number==53){V(1)=TV();V(2).x=0;V(3).x=0;V(4).x=0;}
    if(test_number==54 || test_number==55){V(1)=V(2)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==31)
    {
        int number_of_constrained_particles = constrained_particles.m;
        for (int i=1; i<=number_of_constrained_particles; i++)
            V(constrained_particles(i))=TV();
    }
    if(test_number==24){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int j=n/3+1;j<=2*n/3+1;j++){V(i+m*(j-1))=TV();V(i+m*(j-1)+(mn-1)*m*n)=TV();}
        for(int i=m/3+1;i<=2*m/3+1;i++)for(int ij=mn/3+1;ij<=2*mn/3+1;ij++){V(i+m*n*(ij-1))=TV();V(i+m*(n-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==57) V.Subset(constrained_particles).Fill(TV());
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
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if(id==kinematic_id) frame=curve.Value(time);
    if(id==kinematic_id2) frame=curve2.Value(time);
    if(id==kinematic_id3) frame=curve3.Value(time);
    if(id==kinematic_id4) frame=curve4.Value(time);
    if(id==kinematic_id5) frame=curve5.Value(time);
    
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(id==kinematic_id) twist=curve.Derivative(time);
    if(id==kinematic_id2) twist=curve2.Derivative(time);
    if(id==kinematic_id3) twist=curve3.Derivative(time);
    if(id==kinematic_id4) twist=curve4.Derivative(time);
    if(id==kinematic_id5) twist=curve5.Derivative(time);
    return false;
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_forces){
        solid_body_collection.deformable_body_collection.Test_Energy(time);
        solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}
    if(test_number==10 || test_number==11) solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=10*time;
    if(test_number==33)
    {
        TETRAHEDRALIZED_VOLUME<T>& tet_volume =  solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        FINITE_VOLUME<TV,3>& fvm = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        
        int number_of_vertices = solid_body_collection.deformable_body_collection.collisions.check_collision.m;
        for (int i=1; i<=number_of_vertices; i++) solid_body_collection.deformable_body_collection.collisions.check_collision(i)=false;
        
        for(int t=1;t<=fvm.Fe_hat.m;t++)
            if(fvm.Fe_hat(t).x11<3)
                for (int i=1; i<=4; i++) solid_body_collection.deformable_body_collection.collisions.check_collision(tet_volume.mesh.elements(t)(i))=true;
    }
    if(test_number==58)
    {
        solid_body_collection.deformable_body_collection.collisions.check_collision.Fill(true);
        for(int f=1;FINITE_VOLUME<TV,3>* fvm = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>*>(f);f++)
            for(int t=1;t<=fvm->Fe_hat.m;t++)
                if(fvm->Fe_hat(t).x11>=300)
                    solid_body_collection.deformable_body_collection.collisions.check_collision.Subset(fvm->strain_measure.mesh_object.mesh.elements(t)).Fill(false);
    }
    if(test_number==31)
    {
        TETRAHEDRALIZED_VOLUME<T>& tet_volume = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        FINITE_VOLUME<TV,3>& fvm = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

        int number_of_vertices = solid_body_collection.deformable_body_collection.collisions.check_collision.m;
        for (int i=1; i<=number_of_vertices; i++) solid_body_collection.deformable_body_collection.collisions.check_collision(i)=true;

        for(int t=1;t<=fvm.Fe_hat.m;t++)
            if(fvm.Fe_hat(t).x11>=3)
                for (int i=1; i<=4; i++)
                    if (solid_body_collection.deformable_body_collection.particles.X(tet_volume.mesh.elements(t)(i)).y <= 0.06)
                        solid_body_collection.deformable_body_collection.collisions.check_collision(tet_volume.mesh.elements(t)(i))=false;

        int number_of_constrained_particles = constrained_particles.m;
        for (int i=1; i<=number_of_constrained_particles; i++)
            solid_body_collection.deformable_body_collection.collisions.check_collision(constrained_particles(i))=false;
    }
    if(test_number==51) fish_V=solid_body_collection.deformable_body_collection.particles.V;
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{   if(test_number==29 && time > .1){
        T critical=(T)1.0;
        T critical2=(T)1.0+rebound_time;
        T critical3=(T)1.0+rebound_time+.4;
        T start_young=(T)0; T end_young=(T)rebound_stiffness;
    T pois = (T).45; if(input_poissons_ratio!=-1) pois=input_poissons_ratio;
                DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(time<critical) forces_are_removed=true;
    //if (forces_are_removed){ LOG::cout << "Hey look " << time << std::endl;}
    if(time>critical && forces_are_removed){
        int n=deformable_body_collection.particles.array_collection->Size(); LOG::cout << "Hey look at me " << n << std::endl;
        for (int i=1; i <= n; i++){deformable_body_collection.particles.X(i).y = 10.0;}
                    forces_are_removed=false;
    }
        if(time>critical && time<critical2) {
            FINITE_VOLUME<TV,3>& fv = deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
            CONSTITUTIVE_MODEL<T,3>& icm = fv.constitutive_model;
            T young = pow(10.0,start_young + (time-critical)/(critical2-critical)*(end_young-start_young));
            icm.Update_Lame_Constants(young,pois,(T).01);
            forces_are_removed=false;
        }
    if(time>critical && rebound_time < 1e-6){
        FINITE_VOLUME<TV,3>& fv = deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        CONSTITUTIVE_MODEL<T,3>& icm = fv.constitutive_model;
        icm.Update_Lame_Constants(pow(10.0,end_young),pois,(T).01); 
        forces_are_removed=false;
    }
    if(time>critical3){
        FINITE_VOLUME<TV,3>& fv = deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        CONSTITUTIVE_MODEL<T,3>& icm = fv.constitutive_model;
        T young = pow(10.0,end_young-(T)1.5);
        icm.Update_Lame_Constants(young,pois,(T).01);
        solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;        
    }
    }
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    
    //if (true)
    //{
        FINITE_VOLUME<TV,3>& force_field = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        ARRAY<DIAGONAL_MATRIX<T,3> >& sv = force_field.Fe_hat;
        T Jmin = (T)1; T s1=(T)0; T s2=(T)0; T s3=(T)0;

        for (int i=1; i<=sv.m; i++)
        {
            if (Jmin > sv(i).x11*sv(i).x22*sv(i).x33){ s1=sv(i).x11; s2=sv(i).x22; s3=sv(i).x33; Jmin = sv(i).x11*sv(i).x22*sv(i).x33;}
        }
    LOG::cout<<"Minimum determinant "<<Jmin << " " << s1 << " " << s2 << " " << s3 <<std::endl;
    //}
    if(test_number==51) for(int i=1;i<=particles.X.m;i++) if(particles.V(i).x>6) particles.V(i).x=6;
    T min_volume=FLT_MAX;
    for(int v=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(v);v++){
        for(int i=1;i<=tetrahedralized_volume->mesh.elements.m;i++){
            T vol=tetrahedralized_volume->Signed_Size(i);
            if(vol<min_volume) min_volume=vol;}}
    LOG::cout<<"Minimum tet volume: "<<min_volume<<std::endl;
    if(test_number==29)
        LOG::cout << "Self collisions enabled = " << solids_parameters.triangle_collision_parameters.perform_self_collision << " " << time << std::endl;
    LOG::cout<<"Minimum tet volume: "<<min_volume<<std::endl;
    
}
//#####################################################################
// Function Bind_Intersecting_Particles
//#####################################################################
void Bind_Intersecting_Particles()
{
    if(nobind) return;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    bool added_binding=false;
    for(int i=1;i<=rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){int b=rigid_body_collection.static_and_kinematic_rigid_bodies(i);
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(b);
        for(int p=1;p<=particles.X.m;p++){
            if(!binding_list.Binding(p) && rigid_body.Implicit_Geometry_Lazy_Inside(particles.X(p),1e-3)){
                added_binding=true;
                LOG::cout<<"adding binding: "<<b<<"  "<<p<<std::endl;
                binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,p,rigid_body_collection,b,rigid_body.Object_Space_Point(particles.X(p))));}}}

    if(added_binding) solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame)
{
    dynamic_cast<NEWMARK_EVOLUTION<TV>&>(*solids_evolution).print_matrix=print_matrix;
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    static int first_time=1;
    if (dump_sv)
    {
        std::string output_file = STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d/SV_%d",test_number,frame);
        svout.open(output_file.c_str());
    }

    if(test_number==32 && frame==1100) Bind_Intersecting_Particles();

    if(test_number==52)
    {
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
     
        for (int k=1; k<=jello_centers.m; k++) if (frame==(k-1)*90+1)
        {
            TV center = TV();
            for (int i=1; i<=m*n*mn; i++)
            {
                int index = i+(k-1)*m*n*mn;
                center += particles.X(index);
            }
            center /= m*n*mn;
            for (int i=1; i<=m*n*mn; i++)
            {
                int index = i+(k-1)*m*n*mn;
                particles.X(index) += TV(0,0.2,0)-center;
            }
            for(int i=1;i<=m;i++)
            for(int j=1;j<=n;j++)
            for(int ij=1;ij<=mn;ij++)
            {
                particles.V(i+m*(j-1)+m*n*(ij-1)+(k-1)*m*n*mn) = 0.1*TV(-sin(sin(187*k)*5*i/(T)m),-cos(cos(217*k)*6*j/(T)n),sin(5*ij*sin(471*k)/(T)mn));
            }
        }
    }
    if(test_number==32 && this->restart && first_time && !nobind)
    {
        first_time=0;
        int bind1[]={33187,33237,33289,33476,33609,33617,33621,33628,33632,33633,33635,33638,33639,33643,33644,33649,33798,33799,33817,33819,33832,35723,35724,
                     35725,35728,35733,35734,35736,35737,35738,35753,35756,35757,36138,36148,36149,36159,36161,36162,36168,36177,36188,36193,36194,36195,36196,
                     36197,36198,36199,36201,36202,36203,36204,36205,36206,36207,36209,36210,36211,36212,36213,36214,36221,36232,36233,36235,36237,36238,36241,
                     36242,36248,36249,36250,36251,36262,36270,36273,36281,36282,36283,36285,36336,36337,36338,36339,36363,36367,36370,36371,36575,36576,36760,
                     36765,36766,36775,36776,36823,37806,37807,37808,37811,37813,37825,37826,37827,37828,37830,37831,37832,37835,37836,37837,37838,37839,37840,
                     37841,37846,37854,37867,37868,37869,37870,37871,37872,37874,37875,37877,37878,37879,37881,37882,37883,37884,37916,37918,37919,37920,38220};
        int bind2[]={27664,27716,27718,27827,27839,27872,27876,27880,27882,27884,27892,27895,27907,27913,28221,28240,28314,29095,29182,29379,29380,29381,29382,
                     29413,29468,29469,29470,29471,29472,29473,29474,29482,29488,29511,29573,29574,29577,29578,29595,29597,29599,29602,29609,29611,29621,29630,
                     29633,29634,29636,29637,29645,29646,29648,29714,29715,29717,29718,29719,29720,29724,29725,29727,29728,29731,29732,29742,29743,29744,29745,
                     29810,29822,29823,29835,29836,29837,29838,31573,31588,31619,31620,31630,31633,31636,31637,31670,31671,31672,31673,31674,31676,31679,31680,
                     31682,31687,31689,31693,31694,31695,31696,31697,31698,31699,31700,31737,31738,31739,31741,31742,31744,31747,31748,31750,31949,31999,32760,
                     32761,32762,32763,32764,32766,32767,32768,32770,32771,32772,32773,32795,32796,32797,32798,32799,32815,32816,32817,32818,32819,32823,32824,
                     32826,32831,32832};
        for(size_t i=0;i<sizeof(bind1)/sizeof(*bind1);i++) Add_Debug_Particle(particles.X(bind1[i]),TV(1,0,0));
        for(size_t i=0;i<sizeof(bind2)/sizeof(*bind2);i++) Add_Debug_Particle(particles.X(bind2[i]),TV(0,1,0));

        RIGID_BODY<TV>& torus1=solid_body_collection.rigid_body_collection.Rigid_Body(1);
        RIGID_BODY<TV>& torus2=solid_body_collection.rigid_body_collection.Rigid_Body(2);
        for(size_t i=0;i<sizeof(bind1)/sizeof(*bind1);i++)
            binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,bind1[i],rigid_body_collection,1,torus1.Object_Space_Point(particles.X(bind1[i]))));
        for(size_t i=0;i<sizeof(bind2)/sizeof(*bind2);i++)
            binding_list.Add_Binding(new RIGID_BODY_BINDING<TV>(particles,bind2[i],rigid_body_collection,2,torus2.Object_Space_Point(particles.X(bind2[i]))));
        solid_body_collection.Update_Simulated_Particles();
    }
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
void Add_External_Forces(ARRAY_VIEW<TV> F,const T time) PHYSBAM_OVERRIDE
{
    T v0=4,v1=6;
    if(test_number==50 || test_number==51){
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
        ARRAY<bool> use(particles.X.m);
        use.Subset(externally_forced).Fill(true);
        for(int p=1; p<=particles.X.m; p++)
        {
            T height=particles.X(p).x;
            if(!use(p) && height<10) continue;
            T force_multiplier=sqr(max((T)1,time-(T)1));
            if(height>8.8) force_multiplier*=2;
//            if(height>22) continue;
            if(height>14+time) continue;
            T vm=fish_V(p).Magnitude();
            if(vm>v1) continue;
            if(vm>v0) force_multiplier*=(vm-v0)/(v1-v0);
            F(p)+=TV(1.0*force_multiplier*(14+time-height),0,0);
//            Add_Debug_Particle(particles.X(p),TV(0,0,1));
        }
    }
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,T stiffness,T poissons_ratio,T damping, T cutoff = 0.4, T efc = 20)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* icm=0;
    if(input_efc) efc=input_efc;
    if(input_cutoff) cutoff=input_cutoff;
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;
    
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean2) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_svk) icm=new ST_VENANT_KIRCHHOFF<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_extended_svk) icm=new SVK_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,.6,efc);
    else if(use_extended_mooney_rivlin) icm=new MOONEY_RIVLIN_3D_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,efc);
    else if(use_mooney_rivlin) icm=new MOONEY_RIVLIN_3D2<T>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_extended_neohookean_hyperbola) icm=new NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.1);
    else if(use_extended_neohookean_smooth) icm=new NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.1);
    else if(use_corotated) icm=new COROTATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_blend) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else{
        NEO_HOOKEAN<T,3>* nh=new NEO_HOOKEAN<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
        icm=nh;
        nh->use_constant_ife=use_constant_ife;}
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,icm));
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(dump_sv) svout.close();
}
};
}
#endif
