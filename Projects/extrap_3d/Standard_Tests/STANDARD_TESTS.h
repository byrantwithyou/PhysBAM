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
//   48. Compress with sphere
//   49. See-saw?
//   50. Fish through a torus
//   51. Fish through a tube
//   52  Jello's falling one by one on each other
//   53. Stretch tet
//   54. Stretch tet (II)
//   55. Constrained tet
//   56. Size comparison - several shapes
//   57. Two-direction stretch
//   58. Various objects through gears
//   59. Random armadillo
//   77. Squeeze in a box
//  100. Primary contour field
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
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GEN_NEO_HOOKEAN_ENERGY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/GENERAL_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_J_INTERP_ENERGY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/RC_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/RC2_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ST_VENANT_KIRCHHOFF.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/SVK_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Bindings/RIGID_BODY_BINDING.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/ETHER_DRAG.h>
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
    bool use_extended_neohookean3;
    bool use_int_j_neo;
    bool use_rc_ext;
    bool use_rc2_ext;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_extended_neohookean_smooth;
    bool use_extended_svk, use_svk;
    bool use_corotated;
    bool use_corotated_fixed;
    bool use_mooney_rivlin,use_extended_mooney_rivlin;
    bool use_corot_blend;
    bool dump_sv;
    bool with_bunny,with_hand,with_big_arm,gears_of_pain;
    bool override_collisions,override_no_collisions;
    int kinematic_id,kinematic_id2,kinematic_id3,kinematic_id4,kinematic_id5,kinematic_id6,kinematic_id7,kinematic_id8;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve,curve2,curve3,curve4,curve5,curve6,curve7,curve8;
    bool print_matrix;
    int parameter;
    int fishes,jello_size,number_of_jellos,seed_input;
    T stiffness_multiplier;
    T damping_multiplier;
    T boxsize;
    T degrees_wedge,degrees_incline;
    T rebound_time,rebound_stiffness,rebound_drop;
    bool use_constant_ife;
    bool forces_are_removed,self_collision_flipped,sloped_floor;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    ARRAY<TV> jello_centers;
    T stretch,plateau,repulsion_thickness;
    T hole;
    bool nobind;
    ARRAY<TV> fish_V;
    T input_cutoff;
    T input_efc;
    T input_poissons_ratio,input_youngs_modulus;
    T input_friction,stretch_cutoff;
    T J_min,J_max,la_min;
    T hand_scale;
    bool test_model_only;
    T ether_drag;
    ARRAY<ARRAY<VECTOR<T,2> > > contrail;
    T sigma_range;
    int image_size;
    ARRAY<VECTOR<VECTOR<T,2>,2> > contour_segments;
    ARRAY<int> stuck_particles;
    bool pin_corners;
    RANDOM_NUMBERS<T> rand;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),use_extended_neohookean2(false),
        use_extended_neohookean_refined(false),use_extended_neohookean_hyperbola(false),use_extended_neohookean_smooth(false),use_extended_svk(false),
        use_corotated(false),use_corotated_fixed(false),use_corot_blend(false),dump_sv(false),print_matrix(false),use_constant_ife(false),input_cutoff(FLT_MAX),input_efc(FLT_MAX),
        input_poissons_ratio(-1),input_youngs_modulus(0),J_min(0),J_max((T).1),la_min(0),test_model_only(false),ether_drag(0),pin_corners(true)
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
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {if((test_number==60 || test_number==17 || test_number==18) && time<1e-5) dt=std::min(dt,3e-7);}
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
    parse_args->Add_Option_Argument("-use_ext_neo3");
    parse_args->Add_Option_Argument("-use_int_j_neo");
    parse_args->Add_Option_Argument("-use_rc_ext");
    parse_args->Add_Option_Argument("-use_rc2_ext");
    parse_args->Add_Option_Argument("-use_ext_neo_ref");
    parse_args->Add_Option_Argument("-use_ext_neo_hyper");
    parse_args->Add_Option_Argument("-use_ext_neo_smooth");
    parse_args->Add_Option_Argument("-use_ext_svk");
    parse_args->Add_Option_Argument("-use_ext_mooney");
    parse_args->Add_Option_Argument("-use_mooney");
    parse_args->Add_Option_Argument("-use_corotated");
    parse_args->Add_Option_Argument("-use_corotated_fixed");
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
    parse_args->Add_Integer_Argument("-seed",1234,"random seed to use");
    parse_args->Add_Integer_Argument("-solver_iterations",1000,"number of iterations used for solids system");
    parse_args->Add_Option_Argument("-use_constant_ife","use constant extrapolation on inverting finite element fix");
    parse_args->Add_Option_Argument("-test_system");
    parse_args->Add_Option_Argument("-collisions","Does not yet work in all sims, see code for details");
    parse_args->Add_Option_Argument("-no_collisions","Does not yet work in all sims, see code for details");
    parse_args->Add_Double_Argument("-stretch",1,"stretch");
    parse_args->Add_Double_Argument("-hole",.5,"hole");
    parse_args->Add_Double_Argument("-rebound_time",.2,"number of seconds to rebound in test 29");
    parse_args->Add_Double_Argument("-rebound_stiffness",5,"log10 of youngs modulus of final stiffness");
    parse_args->Add_Double_Argument("-rebound_drop",1.5,"log10 of youngs modulus of dropoff of final stiffness");
    parse_args->Add_Option_Argument("-nobind");
    parse_args->Add_Double_Argument("-cutoff",.4,"cutoff");
    parse_args->Add_Double_Argument("-repulsion_thickness",1e-4,"repulsion thickness");
    parse_args->Add_Double_Argument("-efc",20,"efc");
    parse_args->Add_Double_Argument("-poissons_ratio",-1,"poissons_ratio");
    parse_args->Add_Double_Argument("-youngs_modulus",0,"youngs modulus, only for test 41 so far");
    parse_args->Add_Integer_Argument("-jello_size",20,"resolution of each jello cube");
    parse_args->Add_Integer_Argument("-number_of_jellos",12,"number of falling jello cubes in test 41");
    parse_args->Add_Double_Argument("-degrees_incline",5.0,"degrees of incline");
    parse_args->Add_Double_Argument("-degrees_wedge",2.0,"degrees of side incline");
    parse_args->Add_Double_Argument("-friction",.3,"amount of friction");
    parse_args->Add_Option_Argument("-gears_of_pain");
    parse_args->Add_Option_Argument("-sloped_floor");
    parse_args->Add_Double_Argument("-ja",.6,"J to stop interpolating");
    parse_args->Add_Double_Argument("-jb",.4,"J to start interpolating");
    parse_args->Add_Double_Argument("-lb",-.1,"final poisson's ratio");
    parse_args->Add_Double_Argument("-hand_scale",.8,"hand scale on test 58");
    parse_args->Add_Option_Argument("-test_model_only");
    parse_args->Add_Double_Argument("-stretch_cutoff",300,"Stretch cutoff on test 58");
    parse_args->Add_Double_Argument("-ether_drag",0,"Ether drag");
    parse_args->Add_Integer_Argument("-image_size",500,"image size for plots");
    parse_args->Add_Double_Argument("-sigma_range",3,"sigma range for plots");
    parse_args->Add_Option_Argument("-pin_corners");
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
    plateau=(T)0;
    parameter=parse_args->Get_Integer_Value("-parameter");
    jello_size=parse_args->Get_Integer_Value("-jello_size");
    
    switch(test_number){
        case 17: case 18: case 24: case 25: case 27: case 10: case 11: case 23: case 57: case 77: case 80:
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
            mattress_grid2=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.016,(T)0.016,(T)-0.016,(T)0.016,(T)-0.016,(T)0.016);
            mattress_grid3=GRID<TV>(jello_size,jello_size,jello_size,(T)-0.0125,(T)0.0125,(T)-0.0125,(T)0.0125,(T)-0.0125,(T)0.0125);
            break;
        case 37: case 39: case 40: case 38: case 44:
            mattress_grid=GRID<TV>(2,2,2,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            break;
        case 42: case 52:
            mattress_grid=GRID<TV>(20,20,20,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01,(T)-0.01,(T)0.01);
            break;
    	default:{
            if(!parameter) parameter=10;            
        
            mattress_grid=GRID<TV>(2*parameter,parameter,2*parameter,(T)-1,(T)1,(T)-.5,(T).5,(T)-1,(T)1);
        }
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
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    stiffness_multiplier=(T)parse_args->Get_Double_Value("-stiffen");
    damping_multiplier=(T)parse_args->Get_Double_Value("-dampen");
    stretch=(T)parse_args->Get_Double_Value("-stretch");
    input_friction=(T)parse_args->Get_Double_Value("-friction");
    test_forces=parse_args->Is_Value_Set("-test_forces");
    use_extended_neohookean=parse_args->Is_Value_Set("-use_ext_neo");
    use_extended_neohookean2=parse_args->Is_Value_Set("-use_ext_neo2");
    use_extended_neohookean3=parse_args->Is_Value_Set("-use_ext_neo3");
    use_int_j_neo=parse_args->Is_Value_Set("-use_int_j_neo");
    use_rc_ext=parse_args->Is_Value_Set("-use_rc_ext");
    use_rc2_ext=parse_args->Is_Value_Set("-use_rc2_ext");
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_ref");
    use_extended_neohookean_hyperbola=parse_args->Is_Value_Set("-use_ext_neo_hyper");
    use_extended_mooney_rivlin=parse_args->Is_Value_Set("-use_ext_mooney");
    use_mooney_rivlin=parse_args->Is_Value_Set("-use_mooney");
    use_extended_neohookean_smooth=parse_args->Is_Value_Set("-use_ext_neo_smooth");
    use_svk=parse_args->Is_Value_Set("-use_svk");
    use_extended_svk=parse_args->Is_Value_Set("-use_ext_svk");
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_corotated_fixed=parse_args->Is_Value_Set("-use_corotated_fixed");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    use_constant_ife=parse_args->Get_Option_Value("-use_constant_ife");
    solids_parameters.implicit_solve_parameters.test_system=parse_args->Is_Value_Set("-test_system");
    override_collisions=parse_args->Is_Value_Set("-collisions");
    override_no_collisions=parse_args->Is_Value_Set("-no_collisions")&&(!override_collisions);
    hole=(T)parse_args->Get_Double_Value("-hole");
    rebound_stiffness=(T)parse_args->Get_Double_Value("-rebound_stiffness");
    rebound_time=(T)parse_args->Get_Double_Value("-rebound_time");
    rebound_drop=(T)parse_args->Get_Double_Value("-rebound_drop");
    with_bunny=parse_args->Is_Value_Set("-with_bunny");
    with_hand=parse_args->Is_Value_Set("-with_hand");
    with_big_arm=parse_args->Is_Value_Set("-with_big_arm");
    number_of_jellos=parse_args->Get_Integer_Value("-number_of_jellos");
    degrees_incline=parse_args->Get_Double_Value("-degrees_incline");
    degrees_wedge=parse_args->Get_Double_Value("-degrees_wedge");
    gears_of_pain=parse_args->Is_Value_Set("-gears_of_pain");
    sloped_floor=parse_args->Is_Value_Set("-sloped_floor");
    J_min=(T)parse_args->Get_Double_Value("-ja");
    J_max=(T)parse_args->Get_Double_Value("-jb");
    la_min=(T)parse_args->Get_Double_Value("-lb");
    hand_scale=(T)parse_args->Get_Double_Value("-hand_scale");
    test_model_only=parse_args->Get_Option_Value("-test_model_only");
    ether_drag=(T)parse_args->Get_Double_Value("-ether_drag");
    repulsion_thickness=(T)parse_args->Get_Double_Value("-repulsion_thickness");
    sigma_range=(T)parse_args->Get_Double_Value("-sigma_range");
    image_size=parse_args->Get_Integer_Value("-image_size");
    stretch_cutoff=parse_args->Get_Double_Value("-stretch_cutoff");
    pin_corners=parse_args->Get_Option_Value("-pin_corners");
    
    semi_implicit=parse_args->Is_Value_Set("-semi_implicit");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
    nobind=parse_args->Is_Value_Set("-nobind");
    if(parse_args->Is_Value_Set("-cutoff")) input_cutoff=(T)parse_args->Get_Double_Value("-cutoff");
    if(parse_args->Is_Value_Set("-efc")) input_efc=(T)parse_args->Get_Double_Value("-efc");
    if(parse_args->Is_Value_Set("-poissons_ratio")) input_poissons_ratio=(T)parse_args->Get_Double_Value("-poissons_ratio");
    if(parse_args->Is_Value_Set("-youngs_modulus")) input_youngs_modulus=(T)parse_args->Get_Double_Value("-youngs_modulus");
    if(parse_args->Is_Value_Set("-seed")) rand.Set_Seed(parse_args->Get_Integer_Value("-seed"));

    switch(test_number){
        case 1:
        case 2:
        case 3:
       // case 4:
        case 7:
        case 8:
        case 80:
        case 9:
        case 10:
        case 11:
        case 16:
        case 17:
        case 18:
        case 56:
        case 77:
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
        case 27: case 23: case 53: case 54: case 55: case 57: case 100: case 48:
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
        case 29: case 4:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-5;

            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;//This gets turned off later then back on
            //std::cout << "rame collisions are " << override_collisions << std::endl;
            self_collision_flipped=false;
            //}
            frame_rate=120;
            last_frame=420;
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
        case 59:
        case 60:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 1e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=false;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=3000;
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
            solids_parameters.cfl=(T)10;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 3e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            frame_rate=120;
            last_frame=10*300;
            break;
        case 43:
        case 58:
            solids_parameters.cfl=(T)10;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = repulsion_thickness;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=true;
            solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
            frame_rate=120;
            last_frame=1500;
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
            solids_parameters.cfl=(T)5;
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 2e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;
            
            if (override_no_collisions){
                solids_parameters.triangle_collision_parameters.perform_self_collision=false;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;                
            }
            frame_rate=600;
            last_frame=1000;
            break;
        case 52:
            solids_parameters.triangle_collision_parameters.collisions_repulsion_thickness = 2e-4;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=true;
            solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=true;

            if (override_no_collisions){
                solids_parameters.triangle_collision_parameters.perform_self_collision=false;
                solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
                solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;                
            }
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
    
    //solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
    //solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
    //solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;

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
        case 80:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_20K.tet.gz", // adaptive_torus_float.tet
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(15,5,0),ROTATION<TV>(-T(pi/2),TV(1,0,0)))),true,true,density,1.0);
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
            for(int i=0;i<7;i++){
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
        case 77: {
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,0,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);

            RIGID_BODY<TV>& box_bottom=tests.Add_Analytic_Box(TV(6,2,6));
            RIGID_BODY<TV>& box_side_1=tests.Add_Analytic_Box(TV(2,6,6));
            RIGID_BODY<TV>& box_side_2=tests.Add_Analytic_Box(TV(2,6,6));
            RIGID_BODY<TV>& box_side_3=tests.Add_Analytic_Box(TV(6,6,2));
            RIGID_BODY<TV>& box_side_4=tests.Add_Analytic_Box(TV(6,6,2));
            RIGID_BODY<TV>& box_top=tests.Add_Analytic_Box(TV(6,2,6));
                        
            box_bottom.X()=TV(0,-2,0);
            box_side_1.X()=TV(-2,0,0);
            box_side_2.X()=TV(2,0,0);
            box_side_3.X()=TV(0,0,-2);
            box_side_4.X()=TV(0,0,2);

            box_bottom.is_static=true;
            box_side_1.is_static=true;
            box_side_2.is_static=true;
            box_side_3.is_static=true;
            box_side_4.is_static=true;
            box_top.is_static=false;

            box_bottom.coefficient_of_friction=0;
            box_side_1.coefficient_of_friction=0;
            box_side_2.coefficient_of_friction=0;
            box_side_3.coefficient_of_friction=0;
            box_side_4.coefficient_of_friction=0;
            box_top.coefficient_of_friction=0;

            kinematic_id=box_top.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box_top.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,2,0)));
            curve.Add_Control_Point(10,FRAME<TV>(TV(0,0,0)));
            last_frame=300;
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
        case 59: case 60: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)0.373,(T)0))),true,true,density,.005);
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
            T y_translate=0.1;

            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(-0.2,2.5-y_translate,-0.04),ROTATION<TV>(T(-pi/2.*0.9),TV(0,0,1))*ROTATION<TV>(T(pi/2.),TV(0,1,0)))),true,true,density,0.013);
            
            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375,1);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375,1);
            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1.5,.15,24);
 
            gear1.coefficient_of_friction = 1;
            gear2.coefficient_of_friction = 1;
            cylinder.coefficient_of_friction = 0;
 
            kinematic_id=gear1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear1.particle_index)=true;
            kinematic_id2=gear2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear2.particle_index)=true;
            kinematic_id3=cylinder.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder.particle_index)=true;
             
            T angular_velocity = 1;
            T scale = 1;
 
            curve3.Add_Control_Point(0,FRAME<TV>(TV(0,4-y_translate,0),ROTATION<TV>(0,TV(0,0,1))));
            curve3.Add_Control_Point(2,FRAME<TV>(TV(0,4-y_translate,0),ROTATION<TV>(0,TV(0,0,1))));
            curve3.Add_Control_Point(4,FRAME<TV>(TV(0,2-y_translate,0),ROTATION<TV>(0,TV(0,0,1))));
 
            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-(T).4*scale,1.5*scale-y_translate,-.75*scale),ROTATION<TV>(-i,TV(0,0,1))));
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV((T).4*scale,1.5*scale-y_translate,-.75*scale),ROTATION<TV>(i,TV(0,0,1))));}

            tests.Add_Ground();
            break;}
            
        case 58:{
            T scale = 0.75;
            
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3*scale,0),ROTATION<TV>(T(pi/2),TV(0,1,0))*ROTATION<TV>(T(pi/2),TV(1,0,0)))),true,true,density,hand_scale);
            if(!gears_of_pain){ //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/fish_42K.tet",
                //  RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,4.5*scale,-1.2*scale),ROTATION<TV>((T)pi*0.525,TV(1,0,0))*ROTATION<TV>(0*(T)pi/2,TV(0,1,0)))),true,true,density,0.06);            
                //tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/bunny.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.7*scale,-3.0*scale))),true,true,density,.25);
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.6*scale,-0.0*scale),ROTATION<TV>(T(-pi/2),TV(0,0,1))*ROTATION<TV>(T(pi/2),TV(0,1,0)))),true,true,density,.0065);                
            }else{
                tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)4.6*scale,-0.0*scale),ROTATION<TV>(T(-pi/2),TV(0,0,1))*ROTATION<TV>(T(pi/2),TV(0,1,0)))),true,true,density,.005);      
            }
            
            
            RIGID_BODY<TV>& gear1=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            RIGID_BODY<TV>& gear2=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            
            gear1.coefficient_of_friction = input_friction;
            gear2.coefficient_of_friction = input_friction;
            RIGID_BODY<TV>& gear3=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            RIGID_BODY<TV>& gear4=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            
            gear3.coefficient_of_friction = input_friction;
            gear4.coefficient_of_friction = input_friction;
            RIGID_BODY<TV>& gear5=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            RIGID_BODY<TV>& gear6=tests.Add_Rigid_Body("gear",.375*scale,1.0*scale);
            
            gear5.coefficient_of_friction = input_friction;
            gear6.coefficient_of_friction = input_friction;
            
            kinematic_id=gear1.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear1.particle_index)=true;
            kinematic_id2=gear2.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear2.particle_index)=true;
            kinematic_id3=gear3.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear3.particle_index)=true;
            kinematic_id4=gear4.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear4.particle_index)=true;
            kinematic_id5=gear5.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear5.particle_index)=true;
            kinematic_id6=gear6.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(gear6.particle_index)=true; 
            
            T angular_velocity = 1; T gear_dx;
            if (gears_of_pain) gear_dx=(T).377; else gear_dx=(T).4;
            
            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,-.75*scale),ROTATION<TV>(-i,TV(0,0,1)))); //.4 is default, gears barely touch at .375
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,-.75*scale),ROTATION<TV>(i,TV(0,0,1))));
                curve3.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,-1.75*scale),ROTATION<TV>(-i,TV(0,0,1)))); //.4 is default, gears barely touch at .375
                curve4.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,-1.75*scale),ROTATION<TV>(i,TV(0,0,1))));
                curve5.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-gear_dx*scale,1.5*scale,.25*scale),ROTATION<TV>(-i,TV(0,0,1)))); //.4 is default, gears barely touch at .375
                curve6.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(gear_dx*scale,1.5*scale,.25*scale),ROTATION<TV>(i,TV(0,0,1))));            
            }
            
            RIGID_BODY<TV>& box0=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale));
            
            //  RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale));
            //  RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(2.0*scale,.1*scale,2.0*scale));
            RIGID_BODY<TV>& cylinder=tests.Add_Analytic_Cylinder(1.5*scale,.12*scale);
            box0.X()=TV(0,4.0*scale,-0.0*scale);
            
            /*  if(!gears_of_pain){
             RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(2.0*scale,2.0*scale,.1*scale));
             box1.X()=TV(0,1.2*scale,0.831*scale);
             //box1.Rotation()=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
             box1.is_static=true;
             }*/
            
            //box2.X()=TV(0,1.2*scale,-0.831*scale);
            // box2.Rotation()=ROTATION<TV>(-(T)pi/4.0,TV(1,0,0));
            // box3.X()=TV(0,3.0*scale,-1.4*scale);
            
            // box3.Rotation()=ROTATION<TV>((T)pi/4.0,TV(1,0,0));
            box0.is_static=false; //Will move later
            // box2.is_static=true;
            // box3.is_static=false;
            
            /* box3.coefficient_of_friction = .5;
             kinematic_id3=box3.particle_index;
             rigid_body_collection.rigid_body_particle.kinematic(box3.particle_index)=true; 
             curve3.Add_Control_Point(0,FRAME<TV>(TV(0,3.0*scale,-1.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));
             curve3.Add_Control_Point(.1,FRAME<TV>(TV(0,3.0*scale,-1.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));
             curve3.Add_Control_Point(.11,FRAME<TV>(TV(0,3.0*scale,-3.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));*/
            
            box0.coefficient_of_friction = .05;
            kinematic_id8=box0.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box0.particle_index)=true; 
            curve8.Add_Control_Point(0,FRAME<TV>(TV(0,4.0*scale,-0.0*scale)));
            curve8.Add_Control_Point(1.5+hand_scale,FRAME<TV>(TV(0,4.0*scale,-0.0*scale),ROTATION<TV>(0*(T)pi/2.0,TV(1,0,0))));
            curve8.Add_Control_Point(1.53+hand_scale,FRAME<TV>(TV(0,3.1*scale,-1.0*scale),ROTATION<TV>((T)pi/2.0,TV(1,0,0))));
            
            T drop_time;
            if(gears_of_pain) drop_time=(T)17; else drop_time=(T)17;
            cylinder.X()=TV(0,5.0*scale,0*scale);
            cylinder.is_static=false;
            kinematic_id7=cylinder.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(cylinder.particle_index)=true;
            curve7.Add_Control_Point(0,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curve7.Add_Control_Point(drop_time,FRAME<TV>(TV(0,5.0*scale,0*scale)));
            curve7.Add_Control_Point(drop_time+(T)1,FRAME<TV>(TV(0,2.0*scale,0*scale)));
            
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);            
            inclined_floor.X()=TV(0,.0*scale,0);
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
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(0),ROTATION<TV>(T(pi/0.13),TV(1.3,1.5,0.7))));
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(jello_centers(1),ROTATION<TV>(T(pi/0.076),TV(0.7,1,0.1))));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground();
            break;
        }
        case 41:

        {
            //Todo: At some point, use Oriented_Box to do a tighter collision check. C-Re, C-e, R where R is the rotation matrix and e is the edge
            T max_jello_size = .036;//maximum edge length
            T bound = .2;
            TV new_center; T new_rotate;
            T board_height = .6-.1*bound; plateau=board_height;
            bool stuck=false;
            RIGID_BODY_STATE<TV> initial_state;
            
            //break;} 
            for (int i=1; i<=number_of_jellos; i++){
                do {
                    rand.Fill_Uniform(new_center,-bound,bound);
                   // new_center = TV(rand.Get_Uniform_Number(-bound,bound),rand.Get_Uniform_Number((T).5*bound,(T)1*bound),rand.Get_Uniform_Number(-bound,bound));
                    stuck=false;
                    new_center.x = .6*(new_center.x);
                    new_center.z = 1.1*(new_center.z+bound)-bound;
                    if (new_center.z < -.5*bound) new_center.y = (T).5*(new_center.y + 2.0*bound+3.0*(new_center.z+bound));
                    else new_center.y = .3*(new_center.y+bound) + .5*max_jello_size + board_height+.2*bound;
                    for (int j=1; j<i&&(!stuck); j++){
                        //LOG::cout << i << " " << j << " " << new_center << " " << jello_centers(j) << " " << (new_center-jello_centers(j)).Magnitude() << " " << (T)8*max_jello_size*max_jello_size << std::endl;
                        if((new_center-jello_centers(j)).Magnitude()<=(T)2*max_jello_size) stuck=true;
                       // if ((new_center.x-jello_centers(j).x)*(new_center.x-jello_centers(j).x)+(new_center.y-jello_centers(j).y)*(new_center.y-jello_centers(j).y)+(new_center.z-jello_centers(j).z)*(new_center.z-jello_centers(j).z)<=(T)4*max_jello_size*max_jello_size) stuck=true;
                }}while(stuck);
                jello_centers.Append(new_center);
                rand.Fill_Uniform(new_center,-bound,bound);
//                     new_center = TV(rand.Get_Uniform_Number(-bound,bound),rand.Get_Uniform_Number(-bound,bound),rand.Get_Uniform_Number(-bound,bound));
                new_rotate = rand.Get_Uniform_Number(-(T)pi,(T)pi);
                RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(jello_centers(i),ROTATION<TV>(new_rotate,new_center)));
                if (i % 5 ==0) {tests.Create_Mattress(mattress_grid3,true,&initial_state1);}
                else if (i % 5 ==1) {tests.Create_Mattress(mattress_grid2,true,&initial_state1);}
                else {tests.Create_Mattress(mattress_grid1,true,&initial_state1);}
            }
            if(sloped_floor){
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Rotation()=ROTATION<TV>((T)pi*degrees_wedge/(T)180,TV(0,0,1))*ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            RIGID_BODY<TV>& inclined_floor2=tests.Add_Ground(input_friction);
            inclined_floor2.Rotation()=ROTATION<TV>(-(T)pi*degrees_wedge/(T)180,TV(0,0,1))*ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            }
            else
            {
            RIGID_BODY<TV>& inclined_floor=tests.Add_Ground(input_friction);
            inclined_floor.Rotation()=ROTATION<TV>((T)pi*degrees_incline/(T)180,TV(1,0,0));
            }
            
            T dy = .3*bound*sin((T)pi*degrees_incline/(T)180);
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.4*bound));            
            RIGID_BODY<TV>& box2=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box3=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box4=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box5=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box6=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box7=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.2*bound));            
            RIGID_BODY<TV>& box8=tests.Add_Analytic_Box(TV(2.2*bound,.25*bound,.4*bound));            
            box1.X()=TV((T)0,board_height,-.4*bound);
            box2.X()=TV((T)0*bound,board_height,-0.1*bound);
            box3.X()=TV((T)0*bound,board_height,0.1*bound);
            box4.X()=TV((T)0*bound,board_height,0.3*bound);
            box5.X()=TV((T)0*bound,board_height,0.5*bound);
            box6.X()=TV((T)0*bound,board_height,0.7*bound);
            box7.X()=TV((T)0*bound,board_height,0.9*bound);
            box8.X()=TV((T)0*bound,board_height,1.2*bound);

            curve.Add_Control_Point(0,FRAME<TV>(box1.X()));
            curve.Add_Control_Point(.2,FRAME<TV>(box1.X()));
            //curve.Add_Control_Point(.21,FRAME<TV>(box1.X()-TV(0,.2*board_height,0)));
            //curve.Add_Control_Point(.23,FRAME<TV>(box1.X()-TV(-2.5*bound,.2*board_height,0)));
            curve2.Add_Control_Point(0,FRAME<TV>(box2.X()));
            curve2.Add_Control_Point(.2,FRAME<TV>(box2.X()));
            curve2.Add_Control_Point(.27,FRAME<TV>(box2.X()-TV(0,dy,0)));
            //curve2.Add_Control_Point(.31,FRAME<TV>(box2.X()-TV(0,dy+.2*board_height,0)));
            //curve2.Add_Control_Point(.33,FRAME<TV>(box2.X()-TV(-2.5*bound,dy+.2*board_height,0)));
            curve3.Add_Control_Point(0,FRAME<TV>(box3.X()));
            curve3.Add_Control_Point(.2,FRAME<TV>(box3.X()));
            curve3.Add_Control_Point(.33,FRAME<TV>(box3.X()-TV(0,2.0*dy,0)));
            //curve3.Add_Control_Point(.41,FRAME<TV>(box3.X()-TV(0,2.0*dy+.2*board_height,0)));
            //curve3.Add_Control_Point(.43,FRAME<TV>(box3.X()-TV(-2.5*bound,2.0*dy+.2*board_height,0)));
            curve4.Add_Control_Point(0,FRAME<TV>(box4.X()));
            curve4.Add_Control_Point(.2,FRAME<TV>(box4.X()));
            curve4.Add_Control_Point(.40,FRAME<TV>(box4.X()-TV(0,3.0*dy,0)));
            //curve4.Add_Control_Point(.51,FRAME<TV>(box4.X()-TV(0,3.0*dy+.2*board_height,0)));
            //curve4.Add_Control_Point(.53,FRAME<TV>(box4.X()-TV(-2.5*bound,3.0*dy+.2*board_height,0)));
            curve5.Add_Control_Point(0,FRAME<TV>(box5.X()));
            curve5.Add_Control_Point(.2,FRAME<TV>(box5.X()));
            curve5.Add_Control_Point(.47,FRAME<TV>(box5.X()-TV(0,4.0*dy,0)));
            //curve5.Add_Control_Point(.61,FRAME<TV>(box5.X()-TV(0,4.0*dy+.2*board_height,0)));
            //curve5.Add_Control_Point(.63,FRAME<TV>(box5.X()-TV(-2.5*bound,4.0*dy+.2*board_height,0)));
            curve6.Add_Control_Point(0,FRAME<TV>(box6.X()));
            curve6.Add_Control_Point(.2,FRAME<TV>(box6.X()));
            curve6.Add_Control_Point(.53,FRAME<TV>(box6.X()-TV(0,5.0*dy,0)));
            //curve6.Add_Control_Point(.71,FRAME<TV>(box6.X()-TV(0,5.0*dy+.2*board_height,0)));
            //curve6.Add_Control_Point(.73,FRAME<TV>(box6.X()-TV(-2.5*bound,5.0*dy+.2*board_height,0)));
            curve7.Add_Control_Point(0,FRAME<TV>(box7.X()));
            curve7.Add_Control_Point(.2,FRAME<TV>(box7.X()));
            curve7.Add_Control_Point(.60,FRAME<TV>(box7.X()-TV(0,6.0*dy,0)));
            //curve7.Add_Control_Point(.81,FRAME<TV>(box7.X()-TV(0,6.0*dy+.2*board_height,0)));
            //curve7.Add_Control_Point(.83,FRAME<TV>(box7.X()-TV(-2.5*bound,6.0*dy+.2*board_height,0)));
            curve8.Add_Control_Point(0,FRAME<TV>(box8.X()));
            curve8.Add_Control_Point(.2,FRAME<TV>(box8.X()));
            curve8.Add_Control_Point(.67,FRAME<TV>(box8.X()-TV(0,7.0*dy,0)));
            //curve8.Add_Control_Point(.91,FRAME<TV>(box8.X()-TV(0,7.0*dy+.2*board_height,0)));
            //curve8.Add_Control_Point(.93,FRAME<TV>(box8.X()-TV(-2.5*bound,7.0*dy+.2*board_height,0)));
    /*        curve2.Add_Control_Point(0,FRAME<TV>(box2.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve2.Add_Control_Point(.3,FRAME<TV>(box2.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve2.Add_Control_Point(.35,FRAME<TV>(box2.X(),ROTATION<TV>(-(T)pi/2.0,TV(0,0,1))));
            curve3.Add_Control_Point(.0,FRAME<TV>(box3.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve3.Add_Control_Point(.4,FRAME<TV>(box3.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve3.Add_Control_Point(.45,FRAME<TV>(box3.X(),ROTATION<TV>(-(T)pi/2.0,TV(0,0,1))));
            curve4.Add_Control_Point(.0,FRAME<TV>(box4.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve4.Add_Control_Point(.5,FRAME<TV>(box4.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve4.Add_Control_Point(.55,FRAME<TV>(box4.X(),ROTATION<TV>(-(T)pi/2.0,TV(0,0,1))));
            curve5.Add_Control_Point(.0,FRAME<TV>(box5.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve5.Add_Control_Point(.6,FRAME<TV>(box5.X(),ROTATION<TV>(0*(T)pi/2.0,TV(0,0,1))));
            curve5.Add_Control_Point(.65,FRAME<TV>(box5.X(),ROTATION<TV>(-(T)pi/2.0,TV(0,0,1))));*/
            
               //         curve3.Add_Control_Point(0,FRAME<TV>(TV(0,3.0*scale,-1.4*scale),ROTATION<TV>((T)pi/4.0,TV(1,0,0))));
            box1.is_static=false;
            kinematic_id=box1.particle_index;
            box2.is_static=false;
            kinematic_id2=box2.particle_index;
            box3.is_static=false;
            kinematic_id3=box3.particle_index;
            box4.is_static=false;
            kinematic_id4=box4.particle_index;
            box5.is_static=false;
            kinematic_id5=box5.particle_index;
            box6.is_static=false;
            kinematic_id6=box6.particle_index;
            box7.is_static=false;
            kinematic_id7=box7.particle_index;
            box8.is_static=false;
            kinematic_id8=box8.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(box1.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box2.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box3.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box4.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box5.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box6.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box7.particle_index)=true;
            rigid_body_collection.rigid_body_particle.kinematic(box8.particle_index)=true;
            
             
            
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
            T flat_bottom_radius=0.02;
            T rounded_innterior_radius=0.02;
            T vertical_straight_length=0.12;
            T thickness=0.01;
            T coefficient_of_friction=0.3;
            RIGID_BODY<TV>& rounding=tests.Add_Analytic_Bowl(flat_bottom_radius,rounded_innterior_radius,thickness,128,32);
            rounding.coefficient_of_friction = coefficient_of_friction;
            rounding.Rotation()=ROTATION<TV>((T)pi,TV(1,0,0));
            rounding.X()=TV(0,rounded_innterior_radius+thickness,0);
            rounding.is_static=true;

            RIGID_BODY<TV>& walls=tests.Add_Analytic_Shell(vertical_straight_length,flat_bottom_radius+rounded_innterior_radius+thickness,flat_bottom_radius+rounded_innterior_radius,128);
            walls.coefficient_of_friction = coefficient_of_friction;
            walls.Rotation()=ROTATION<TV>((T)pi/2,TV(1,0,0));
            walls.X()=TV(0,vertical_straight_length/2,0);
            walls.is_static=true;

            RIGID_BODY<TV>& bottom=tests.Add_Analytic_Cylinder(thickness,flat_bottom_radius+rounded_innterior_radius+thickness,128);
            bottom.coefficient_of_friction = coefficient_of_friction;
            bottom.Rotation()=ROTATION<TV>((T)pi/2,TV(1,0,0));
            bottom.X()=TV(0,thickness/2,0);
            bottom.is_static=true;

            int count = 0;
            for (int i=1; i<=25; i++)
            {
                count++;
                jello_centers.Append(TV(i*0.05,50,0));
                RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(jello_centers(count),ROTATION<TV>(10*sin(178*i),TV(sin(145*i),cos(345*i),cos(478*i)))));
                tests.Create_Mattress(mattress_grid,true,&initial_state);
            }
            
            // int number_of_boxes = 128;
            // T height = 0.08;
            // T width = 0.02;
            // T radius = 0.035;

            // ARRAY<RIGID_BODY<TV>*> boxes;
            
            // for (int i=1; i<=number_of_boxes; i++)
            // {
                // T phi = 2*pi*(i-1)/number_of_boxes;

                // boxes.Append(&tests.Add_Analytic_Box(TV(width,height,width)));
                // boxes(i)->X() = TV((radius+width/2)*sin(phi),height/2,(radius+width/2)*cos(phi));
                // boxes(i)->Rotation() = ROTATION<TV>((T)phi,TV(0,1,0));
                // boxes(i)->is_static = true;
                // boxes(i)->coefficient_of_friction = 0.3;
            // }

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
            //torus1.is_static=true;
            torus1.coefficient_of_friction = 0.05;
            torus1.X()=start;
            torus1.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            last_frame=240;
            RIGID_BODY<TV>& torus2=tests.Add_Analytic_Torus((outer-inner)/2,(outer+inner)/2,32,64);
            //torus2.is_static=true;
            torus2.coefficient_of_friction = 0.05;
            torus2.X()=start+TV(length,0,0);
            torus2.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
           // Add_Joint(1,2,new POINT_JOINT<TV>,FRAME<TV>(TV(15,10,0))); // ropes

            /*last_frame=240;
            RIGID_BODY<TV>& shell=tests.Add_Analytic_Shell(length,outer,inner,64);
            shell.is_static=true;
            shell.coefficient_of_friction = 0.05;
            shell.X()=start+TV(length/2,0,0);
            shell.Rotation()=ROTATION<TV>((T)pi/2.0,TV(0,1,0));
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/hand_30k.tet",
                                                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,10,0),ROTATION<TV>(T(pi),TV(0,1,0)))),true,true,density);*/
            
            //Add_Rigid_Body("plank",(T).125,(T).5,0,FRAME<TV>(TV(0,3,0)),"seesaw");
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
        case 48:{
            tests.Create_Mattress(GRID<TV>(TV_INT()+(parameter?parameter:10)+1,RANGE<TV>::Centered_Box()));
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
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(25,5,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            break;}
        case 100:{
            TETRAHEDRALIZED_VOLUME<T>* tv=TETRAHEDRALIZED_VOLUME<T>::Create(particles);
            particles.array_collection->Add_Elements(4);
            particles.X(1)=TV(0,1,0);
/*            particles.X(2)=TV(1/sqrt(3.),0,0);
            particles.X(3)=TV(-1/sqrt(12),0,.5);
            particles.X(4)=TV(-1/sqrt(12),0,-.5);*/
            particles.X(2)=ROTATION<TV>(pi*1./12.,TV(0,1,0)).Rotate(TV(1/sqrt(3.),0,0));
            particles.X(3)=TV(particles.X(2).z,0,particles.X(2).x);
            particles.X(4)=TV(-1/sqrt(6.),0,-1/sqrt(6.)+1e-2);
            tv->mesh.elements.Append(VECTOR<int,4>(1,2,3,4));
            particles.mass.Fill(1);
            solid_body_collection.deformable_body_collection.deformable_geometry.Add_Structure(tv);
            contrail.Resize(1);
            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.deformable_geometry.structures);

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
    PARTICLES<TV>& particles=deformable_body_collection.particles;

    Get_Initial_Data();

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 8:
        case 80:
        case 16:
        case 5:
        case 6:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            if(test_number!=80) solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            if(test_number==80) particles.X.template Project<T,&TV::x>()*=-(T).97;
            if(test_number==80) particles.X.template Project<T,&TV::y>()*=(T).98;
            break;}
        case 77:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
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
            for(int i=0;i<deformable_body_collection.particles.X.m;i++) deformable_body_collection.particles.X(i).y=3;
            break;}
        case 9: case 56:{break;}
        case 10:{
            bool* bools[7]={&use_corotated,0,&use_constant_ife,&use_corot_blend,&use_extended_neohookean,&use_extended_neohookean_smooth,&use_extended_neohookean_hyperbola};
            for(int i=0;i<7;i++){
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
            int s=parameter+1,s2=s*s,s3=s*s2;
            if(pin_corners){
                stuck_particles.Append(1);
                stuck_particles.Append(s);
                stuck_particles.Append(s2);
                stuck_particles.Append(s3);
                stuck_particles.Append(s2-s+1);
                stuck_particles.Append(s3-s2+1);
                stuck_particles.Append(s3-s+1);
                stuck_particles.Append(s3-s2+s);}
            ARRAY<TV> OX(particles.X.Subset(stuck_particles));
            rand.Fill_Uniform(particles.X,-1,1);
            particles.X.Subset(stuck_particles)=OX;
            break;}
        case 18:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
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
                for(int i=0;i<m;i++) for(int j=0;j<n;j++){constrained_particles.Append(i+m*(j-1));constrained_particles.Append(i+m*(j-1)+(p-1)*m*n);}
                for(int i=0;i<m;i++) for(int k=0;k<p;k++){constrained_particles.Append(i+m*n*(k-1));constrained_particles.Append(i+m*(n-1)+m*n*(k-1));}
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
            for(int i=0;i<deformable_body_collection.particles.X.m;i++){ deformable_body_collection.particles.V(i).x=(T)60;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
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
        case 59:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4,(T).3,(T).001);
            
            ARRAY<TV> rand_particles(particles.X.m);
            rand.Fill_Uniform(rand_particles,-0.3,.3);

            for (int i=1; i<=particles.X.m; i++)
            {
                T y   = particles.X(i).y;
                TV v1 = particles.X(i) - TV(0.27,0.645,-0.12);
                TV v2 = particles.X(i) - TV(-0.275,0.605,-0.18);
                
                if ((y<=0.06) || (sqr(v1.x) + sqr(v1.y) + sqr(v1.z) <= sqr(0.065)) || (sqr(v2.x) + sqr(v2.y) + sqr(v2.z) <= sqr(0.065)))
                {
                    constrained_particles.Append(i);
                }
                else
                {
                    particles.X(i)=rand_particles(i)+TV(0,0.3,0)+TV(-0.024,0.06,0.1);
                }
            }

            break;} 
        case 60:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e4,(T).3,(T).001);
            rand.Fill_Uniform(particles.X,-0.3,.3);
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
            Add_Constitutive_Model(tetrahedralized_volume,1e4,0.35,0.005);

            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=0;
            break;}
        case 58:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,2e4,0.4,0.005);
            if(!gears_of_pain){
                TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume2=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(2);
                Add_Constitutive_Model(tetrahedralized_volume2,3e4,0.4,0.005);
                //TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume3=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(3);
                //Add_Constitutive_Model(tetrahedralized_volume3,1e4,0.4,0.005);
            }
             
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
            for(int i=0;i<deformable_body_collection.particles.X.m/2;i++){ deformable_body_collection.particles.V(i).x=(T)1;deformable_body_collection.particles.V(i).y=0;deformable_body_collection.particles.V(i).z=0;}
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
            for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
            for(int ij=0;ij<mn;ij++)
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
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++)
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
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
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
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
                    {
                        particles.V(i+m*j+m*n*ij) = TV(0.32*cos(6*ij/(T)mn)+2.15,-0.28*cos(4*i/(T)m)+0.55,-0.31*sin(4*j/(T)n));
                        particles.V(i+m*j+m*n*ij+m*n*mn) = TV(0.27*sin(6*j/(T)n)-2.17,0.32*cos(6*ij/(T)mn)+0.53,0.25*sin(4*i/(T)m));
                    }
            for (int i=0; i<m*n*mn; i++)
            {
                int index = i;
                particles.V(index) += TV(particles.X(index).y-jello_centers(0).y,-(particles.X(index).x-jello_centers(0).x),0)*51;
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(0).z,-(particles.X(index).y-jello_centers(0).y))*(-13);
                index+=m*n*mn;
                particles.V(index) += TV(particles.X(index).y-jello_centers(1).y,-(particles.X(index).x-jello_centers(1).x),0)*(-39);
                particles.V(index) += TV(0,particles.X(index).z-jello_centers(1).z,-(particles.X(index).y-jello_centers(1).y))*23;
            }
            break;}
        case 41:{
            T youngs_modulus = 1e3;
            T poissons_ratio = .4;
            T damping = 0.001;

            for (int k=0; k<number_of_jellos; k++)
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
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++)
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
            for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                    for(int ij=0;ij<mn;ij++)
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
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);

            break;}
        case 51:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .45;
            T damping = 0.01;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            for(int i=0;i<particles.X.m;i++)
                if(particles.X(i).x>=1.5)
                    externally_forced.Append(i);
            break;}
        case 55:
        case 54:
        case 48:
        case 53:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .45;
            T damping = 0.01;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>(1);
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            if(test_number==55) particles.X(1).x=stretch;
            break;}
        case 100:{
            T youngs_modulus = 1e5;
            T poissons_ratio = .4;
            T damping = 0.1;
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume1=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume1,youngs_modulus,poissons_ratio,damping);
            break;}
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    if(ether_drag) solid_body_collection.Add_Force(new ETHER_DRAG<GRID<TV> >(particles,solid_body_collection.rigid_body_collection,true,true,ether_drag,0));

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        VECTOR<int,3> processes_per_dimension(2,1,1);
        deformable_body_collection.mpi_solids->Simple_Partition(deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,processes_per_dimension);}

    for(int i=0;i<deformable_body_collection.deformable_geometry.structures.m;i++){
        deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.deformable_geometry.structures(i));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision)
            solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(i));}

    if(parse_args->Is_Value_Set("-solver_iterations")) solids_parameters.implicit_solve_parameters.cg_iterations=parse_args->Get_Integer_Value("-solver_iterations");
    if(parse_args->Is_Value_Set("-cgsolids")) solids_parameters.implicit_solve_parameters.cg_tolerance=(T)parse_args->Get_Double_Value("-cgsolids");

    if(!semi_implicit) for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    if(!semi_implicit) for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}

//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
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
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
        for(int i=3*m/7+1;i<=4*m/7+1;i++)for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(i+m*(j-1)+m*n*(ij-1))=-velocity_y;V(i+m*(j-1)+m*n*(ij-1))=-velocity_y;}
    }
    if(test_number==27){
        final_time=70;
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=velocity_x;V(m+m*(j-1)+m*n*(ij-1))=-velocity_x;}
    }
    if(test_number==23){
        final_time=35;
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=-(T)(velocity_time<=25)*velocity_x;V(m+m*(j-1)+m*n*(ij-1))=(T)(velocity_time<=25)*velocity_x;}
    }
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y = velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
    }
    if(test_number==53){V(1)=TV(1,0,0);V(2).x=0;V(3).x=0;V(4).x=0;}
    if(test_number==54){V(1)=TV(1,0,0);V(2)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==55){V(2)=V(1)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==31||test_number==59)
    {
        int number_of_constrained_particles = constrained_particles.m;
        for (int i=1; i<=number_of_constrained_particles; i++)
            V(constrained_particles(i))=TV();
    }
    if(test_number==57) V.Subset(constrained_particles)=constrained_velocities;
    if(test_number==100){V(1)=TV(0,(particles.X(1).y<5-1e-4)*.5,0);V(2).y=0;V(3).y=0;V(4).y=0;}
    if(test_number==48){
        for(int i=0;i<stuck_particles.m;i++){
            int p=stuck_particles(i);
            V(p)=V(p).Projected_Orthogonal_To_Unit_Direction(particles.X(p).Normalized());}}
    if(test_number==17) V.Subset(stuck_particles).Fill(TV());
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
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
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
        for(int i=3*m/7+1;i<=4*m/7+1;i++)for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(i+m*(j-1)+m*n*(ij-1))=TV();V(i+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==23){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==27){
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        for(int ij=0;ij<mn;ij++)for(int j=0;j<n;j++){V(1+m*(j-1)+m*n*(ij-1))=TV();V(m+m*(j-1)+m*n*(ij-1))=TV();}
    }
    if(test_number==53){V(1)=TV();V(2).x=0;V(3).x=0;V(4).x=0;}
    if(test_number==54 || test_number==55){V(1)=V(2)=TV();V(3).x=0;V(4).x=0;}
    if(test_number==31||test_number==59)
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
    if(test_number==100){V(1)=TV();V(2).y=0;V(3).y=0;V(4).y=0;}
    if(test_number==48){
        for(int i=0;i<stuck_particles.m;i++){
            int p=stuck_particles(i);
            V(p)=V(p).Projected_Orthogonal_To_Unit_Direction(particles.X(p).Normalized());}}
    if(test_number==17) V.Subset(stuck_particles).Fill(TV());
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
    if(test_number==41)
    {
        if(id==kinematic_id) 
        {if(time >= .2) frame = FRAME<TV>(TV(1,plateau,0)); else frame=curve.Value(time);}
        if(id==kinematic_id2) 
        {if(time >= .27) frame = FRAME<TV>(TV(2,plateau,0)); else frame=curve2.Value(time);}
        if(id==kinematic_id3) 
        {if(time >= .33) frame = FRAME<TV>(TV(3,plateau,0)); else frame=curve3.Value(time);}
        if(id==kinematic_id4) 
        {if(time >= .4) frame = FRAME<TV>(TV(4,plateau,0)); else frame=curve4.Value(time);}
        if(id==kinematic_id5) 
        {if(time >= .47) frame = FRAME<TV>(TV(5,plateau,0)); else frame=curve5.Value(time);}
        if(id==kinematic_id6) 
        {if(time >= .53) frame = FRAME<TV>(TV(6,plateau,0)); else frame=curve6.Value(time);}
        if(id==kinematic_id7) 
        {if(time >= .6) frame = FRAME<TV>(TV(7,plateau,0)); else frame=curve7.Value(time);}
        if(id==kinematic_id8) 
        {if(time >= .67) frame = FRAME<TV>(TV(8,plateau,0)); else frame=curve8.Value(time);}
        return;
    }
    if(id==kinematic_id) frame=curve.Value(time);
    if(id==kinematic_id2) frame=curve2.Value(time);
    if(id==kinematic_id3) frame=curve3.Value(time);
    if(id==kinematic_id4) frame=curve4.Value(time);
    if(id==kinematic_id5) frame=curve5.Value(time);
    if(id==kinematic_id6) frame=curve6.Value(time);
    if(id==kinematic_id7) frame=curve7.Value(time);
    if(id==kinematic_id8) frame=curve8.Value(time);
    
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
    if(id==kinematic_id6) twist=curve6.Derivative(time);
    if(id==kinematic_id7) twist=curve7.Derivative(time);
    if(id==kinematic_id8) twist=curve8.Derivative(time);
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
        
        for(int t=0;t<fvm.Fe_hat.m;t++)
            if(fvm.Fe_hat(t).x11<3)
                for (int i=0; i<4; i++) solid_body_collection.deformable_body_collection.collisions.check_collision(tet_volume.mesh.elements(t)(i))=true;

        if (time >=1 ) solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=9.8;
    }
    if(test_number==58)
    {
        solid_body_collection.deformable_body_collection.collisions.check_collision.Fill(true);
        for(int f=1;FINITE_VOLUME<TV,3>* fvm = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>*>(f);f++)
            for(int t=0;t<fvm->Fe_hat.m;t++){
               // LOG::cout << "booya " << stretch_cutoff << std::endl;
                if(fvm->Fe_hat(t).x11>=stretch_cutoff)
                    solid_body_collection.deformable_body_collection.collisions.check_collision.Subset(fvm->strain_measure.mesh_object.mesh.elements(t)).Fill(false);}
    }
    if(test_number==29){    
        std::cout << "rame!" <<      solids_parameters.triangle_collision_parameters.perform_self_collision << std::endl;   
        //solids_parameters.triangle_collision_parameters.perform_self_collision=self_collision_flipped;
    }
    if(test_number==31)
    {
        // TETRAHEDRALIZED_VOLUME<T>& tet_volume = solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
        // FINITE_VOLUME<TV,3>& fvm = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();

        // int number_of_vertices = solid_body_collection.deformable_body_collection.collisions.check_collision.m;
        // for (int i=1; i<=number_of_vertices; i++) solid_body_collection.deformable_body_collection.collisions.check_collision(i)=true;

        // for(int t=0;t<fvm.Fe_hat.m;t++)
            // if(fvm.Fe_hat(t).x11>=3)
                // for (int i=1; i<=4; i++)
                    // if (solid_body_collection.deformable_body_collection.particles.X(tet_volume.mesh.elements(t)(i)).y <= 0.06)
                        // solid_body_collection.deformable_body_collection.collisions.check_collision(tet_volume.mesh.elements(t)(i))=false;

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
{   if(test_number==29 && time > .01){
    
    if(solids_parameters.triangle_collision_parameters.perform_self_collision && time<=1.3){
        
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=override_collisions;
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=override_collisions;
        solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;//This gets turned off later then back on

    }

    T critical=(T)1.0;
        T critical2=(T)1.0+rebound_time;
        T critical3=(T)1.0+rebound_time+.3;

    //std::cout << "4th frame Frame" << critical3 <<  std::endl;
    T start_young=(T)0; T end_young=(T)rebound_stiffness;
    T pois = (T).45; if(input_poissons_ratio!=-1) pois=input_poissons_ratio;
                DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(time<critical) forces_are_removed=true;
    //if (forces_are_removed){ LOG::cout << "Hey look " << time << std::endl;}
    if(time>critical && forces_are_removed){
       //  int n=deformable_body_collection.particles.array_collection->Size(); LOG::cout << "Hey look at me " << n << std::endl;
        //for (int i=1; i <= n; i++){deformable_body_collection.particles.X(i).y = 10.0;}
                    forces_are_removed=false;
        
        }
    //std::cout << "3rd critical time frame Frame" << (1>0) << " " << time << " " << critical3 << " " << (time>critical3) << " " << std::endl;

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
    //if(solids_parameters.triangle_collision_parameters.perform_self_collision) std::cout << "rame Hooray!" << std::endl;
    //if(solids_parameters.triangle_collision_parameters.perform_self_collision) std::cout << "rame oo-rah!" << std::endl;
    if(time>critical3 && self_collision_flipped==false){
        self_collision_flipped=true;
        //std::cout << "3rd critical time reached frame Frame" << std::endl;
        FINITE_VOLUME<TV,3>& fv = deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        CONSTITUTIVE_MODEL<T,3>& icm = fv.constitutive_model;
        T young = pow(10.0,end_young-rebound_drop);
        icm.Update_Lame_Constants(young,pois,(T).01);
        solids_parameters.triangle_collision_parameters.perform_self_collision=override_collisions;        
       // solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.deformable_geometry.structures(1));
        if(solids_parameters.triangle_collision_parameters.perform_self_collision) std::cout << "rame oh,rad!" << std::endl;
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

        for(int i=0;i<sv.m;i++){
            if (Jmin > sv(i).x11*sv(i).x22*sv(i).x33){ s1=sv(i).x11; s2=sv(i).x22; s3=sv(i).x33; Jmin = sv(i).x11*sv(i).x22*sv(i).x33;}}
    LOG::cout<<"Minimum determinant "<<Jmin << " " << s1 << " " << s2 << " " << s3 <<std::endl;
    //}
    if(test_number==51) for(int i=0;i<particles.X.m;i++) if(particles.V(i).x>6) particles.V(i).x=6;
    T min_volume=FLT_MAX;
    for(int v=1;TETRAHEDRALIZED_VOLUME<T>* tetrahedralized_volume=solid_body_collection.deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>*>(v);v++){
        for(int i=0;i<tetrahedralized_volume->mesh.elements.m;i++){
            T vol=tetrahedralized_volume->Signed_Size(i);
            if(vol<min_volume) min_volume=vol;}}
    LOG::cout<<"Minimum tet volume: "<<min_volume<<std::endl;
    if(test_number==29)
        LOG::cout << "Self collisions enabled = " << solids_parameters.triangle_collision_parameters.perform_self_collision << " " << time << std::endl;
    if(dump_sv){
        for(int f=1;FINITE_VOLUME<TV,3>* force_field=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>*>(f);f++){
            ARRAY<DIAGONAL_MATRIX<T,3> >& sv = force_field->Fe_hat;
            for(int i=0;i<sv.m;i++){
                svout << sv(i).x11 << " " << sv(i).x22 << " " << sv(i).x33 << std::endl;
                Add_Debug_Particle(sv(i).To_Vector(),TV(1,1,0));
                Debug_Particle_Set_Attribute<TV>(ATTRIBUTE_ID_V,-force_field->isotropic_model->P_From_Strain(sv(i),1,i).To_Vector());}}}
    if(test_number==48){
        stuck_particles.Remove_All();
        T r=std::max((T)0,hole-time*stretch);
        if(time>1.2*hole/stretch) r=std::min((time-1.2*hole/stretch)*stretch,hole);
        for(int i=0;i<particles.X.m;i++)
            if(particles.X(i).Magnitude()>r){
                stuck_particles.Append(i);
                particles.X(i)*=r/particles.X(i).Magnitude();}}
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
    for(int i=0;i<rigid_body_collection.static_and_kinematic_rigid_bodies.m;i++){int b=rigid_body_collection.static_and_kinematic_rigid_bodies(i);
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(b);
        for(int p=0;p<particles.X.m;p++){
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
    if (dump_sv)
    {
        std::string output_file = STRING_UTILITIES::string_sprintf("%s/SV_%d",output_directory.c_str(),frame);
        svout.open(output_file.c_str());
    }

    if(test_number==32 && frame==1100) Bind_Intersecting_Particles();

    if(test_number==52)
    {
        PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
     
        for (int k=1; k<=jello_centers.m; k++) if (frame==(k-1)*60+1)
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
            for(int i=0;i<m;i++)
            for(int j=0;j<n;j++)
            for(int ij=0;ij<mn;ij++)
            {
                particles.V(i+m*(j-1)+m*n*(ij-1)+(k-1)*m*n*mn) = 0.1*TV(-sin(sin(187*k)*5*i/(T)m),-cos(cos(217*k)*6*j/(T)n),sin(5*ij*sin(471*k)/(T)mn));
            }

            for (int k_other=k+1; k_other<=jello_centers.m; k_other++)
            {
                TV center = TV();
                for (int i=1; i<=m*n*mn; i++)
                {
                    int index = i+(k_other-1)*m*n*mn;
                    center += particles.X(index);
                }
                center /= m*n*mn;
                for (int i=1; i<=m*n*mn; i++)
                {
                    int index = i+(k_other-1)*m*n*mn;
                    particles.X(index).y += 50-center.y;
                }
                for(int i=0;i<m;i++)
                for(int j=0;j<n;j++)
                for(int ij=0;ij<mn;ij++)
                {
                    particles.V(i+m*(j-1)+m*n*(ij-1)+(k_other-1)*m*n*mn) = TV();
                }
            }        
        }
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
        for(int p=0;p<particles.X.m;p++)
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
    if(input_efc!=FLT_MAX) efc=input_efc;
    if(input_cutoff!=FLT_MAX) cutoff=input_cutoff;
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;
    
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean2) icm=new NEO_HOOKEAN_EXTRAPOLATED2<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean3) icm=new GENERAL_EXTRAPOLATED<T,3>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_int_j_neo) icm=new GENERAL_EXTRAPOLATED<T,3>(*new NEO_J_INTERP_ENERGY<T>(J_min,J_max,la_min*stiffness*stiffness_multiplier),
        stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc_ext) icm=new RC_EXTRAPOLATED<T,3>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_rc2_ext) icm=new RC2_EXTRAPOLATED<T,3>(*new GEN_NEO_HOOKEAN_ENERGY<T>,stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_svk) icm=new ST_VENANT_KIRCHHOFF<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_extended_svk) icm=new SVK_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,.6,efc);
    else if(use_extended_mooney_rivlin) icm=new MOONEY_RIVLIN_3D_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_mooney_rivlin) icm=new MOONEY_RIVLIN_3D2<T>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_extended_neohookean_hyperbola) icm=new NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc);
    else if(use_extended_neohookean_smooth) icm=new NEO_HOOKEAN_EXTRAPOLATED_SMOOTH<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.1);
    else if(use_corotated) icm=new COROTATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corotated_fixed) icm=new COROTATED_FIXED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_blend) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else{
        NEO_HOOKEAN<T,3>* nh=new NEO_HOOKEAN<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff);
        icm=nh;
        nh->use_constant_ife=use_constant_ife;}
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,icm));
    if(test_model_only) Test_Model(*icm);
}

//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,T a0, T a1, TV da0, TV da1, TV df, T e)
{
    T av=TV::Dot_Product(da1+da0,df)/2/e;
    T dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av, dif, fabs(av-dif));
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,TV a0, TV a1, const MATRIX<T,3>& da0, const MATRIX<T,3>& da1, TV df, T e)
{
    TV av=(da1+da0)*df/2/e;
    TV dif=(a1-a0)/e;
    char buff[1000];
    sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
    LOG::cout<<buff;
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
void Test_Model_Helper(const char* str,const MATRIX<T,3>& a0, const MATRIX<T,3>& a1, const VECTOR<SYMMETRIC_MATRIX<T,3>,3>& da0, const VECTOR<SYMMETRIC_MATRIX<T,3>,3>& da1, TV df, T e)
{
    for(int i=0;i<TV::m;i++){
        TV av=(da1(i)+da0(i))*df/2/e;
        TV dif=(a1.Transposed().Column(i)-a0.Transposed().Column(i))/e;
        char buff[1000];
        sprintf(buff, "============ test ============ %s %8.5f %8.5f (%8.5f)\n", str, av.Magnitude(), dif.Magnitude(), (av-dif).Magnitude());
        LOG::cout<<buff;}
}
//#####################################################################
// Function Test_Model_Helper
//#####################################################################
template<class RC> void 
Test_Model_Helper(ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* icm, TV &f, TV &df, T e)
{
    RC* rc=dynamic_cast<RC*>(icm);
    if(!rc) return;
    if(f.Min()>0) rc->base.Test(f.x,f.y,f.z,0);
    int simplex=0;
    if(f.Product()>rc->extrapolation_cutoff) return;
    typename RC::HELPER h0;
    if(!h0.Compute_E(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,rc->extrapolation_cutoff,f,simplex)) return;
    h0.Compute_dE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f,simplex);
    h0.Compute_ddE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f,simplex);
    typename RC::HELPER h1;
    if(!h1.Compute_E(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,rc->extrapolation_cutoff,f+df,simplex)) return;
    h1.Compute_dE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f+df,simplex);
    h1.Compute_ddE(rc->base,rc->extra_force_coefficient*rc->youngs_modulus,f+df,simplex);
#define XX(k) Test_Model_Helper(#k,h0.k, h1.k, h0.d##k, h1.d##k, df, e);Test_Model_Helper(#k,h0.d##k, h1.d##k, h0.dd##k, h1.dd##k, df, e);
    XX(m);
    XX(h);
    XX(phi);
    XX(E);
    XX(z);
    XX(xi);
    XX(s);
    XX(Q);
    XX(u);
    XX(g);
 }
//#####################################################################
// Function Test_Model
//#####################################################################
void Test_Model(ISOTROPIC_CONSTITUTIVE_MODEL<T,3>& icm)
{
    for(int i=0;i<20;i++){
        TV f;
        rand.Fill_Uniform(f,0,2);
        T e=1e-5;
        TV df;
        rand.Fill_Uniform(df,-e,e);
        f=f.Sorted().Reversed();
        if(rand.Get_Uniform_Integer(0,1)==1) f(2)=-f(2);
        LOG::cout<<f<<std::endl;
        icm.Test(DIAGONAL_MATRIX<T,3>(f),1);
        Test_Model_Helper<RC_EXTRAPOLATED<T,3> >(&icm,f,df,e);}
    if(RC2_EXTRAPOLATED<T,3>* rc2=dynamic_cast<RC2_EXTRAPOLATED<T,3>*>(&icm))
        rc2->Test_Model();
    exit(0);
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    if(dump_sv) svout.close();
    if(test_number==100) Plot_Contour_Landscape(frame);
}
//#####################################################################
// Function Contour_Crossing
//#####################################################################
T Contour_Crossing(const TV& g0,const TV& v0,const TV& g1,const TV& v1)
{
    T a=TV::Dot_Product(g0,v0);
    T b=TV::Dot_Product(g1,v1);
    if(!a) return 0;
    if(!b) return 1;
    if(TV::Dot_Product(v0,v1)<0) b=-b;
    if((a>0) == (b>0)) return -1;
    return a/(a-b);
}
//#####################################################################
// Function Add_Primary_Contour_Segments
//#####################################################################
void Add_Primary_Contour_Segments(ISOTROPIC_CONSTITUTIVE_MODEL<T,3>& icm)
{
    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    ARRAY<TV,VECTOR<int,2> > evec(1,image_size,1,image_size);
    ARRAY<TV,VECTOR<int,2> > grad(1,image_size,1,image_size);
    for(int i=0;i<image_size;i++)
        for(int j=0;j<image_size;j++){
            T x=(2*i-image_size)*sigma_range/image_size+1e-5;
            T y=(2*j-image_size)*sigma_range/image_size;
            TV g=icm.P_From_Strain(DIAGONAL_MATRIX<T,3>(particles.X(1).y,x,y),1,1).To_Vector();
            DIAGONALIZED_ISOTROPIC_STRESS_DERIVATIVE<T,3> disd;
            icm.Isotropic_Stress_Derivative(DIAGONAL_MATRIX<T,3>(particles.X(1).y,x,y),disd,1);
            SYMMETRIC_MATRIX<T,3> H(disd.x1111,disd.x2211,disd.x3311,disd.x2222,disd.x3322,disd.x3333);
            DIAGONAL_MATRIX<T,3> ev;
            MATRIX<T,3> eigenvectors;
            H.Fast_Solve_Eigenproblem(ev,eigenvectors);
            evec(VECTOR<int,2>(i,j))=eigenvectors.Column(ev.To_Vector().Arg_Abs_Max());
            grad(VECTOR<int,2>(i,j))=g;}

    for(int i=1;i<image_size;i++)
        for(int j=1;j<image_size;j++){
            TV g00=grad(VECTOR<int,2>(i,j)),g01=grad(VECTOR<int,2>(i,j+1)),g10=grad(VECTOR<int,2>(i+1,j)),g11=grad(VECTOR<int,2>(i+1,j+1));
            TV v00=evec(VECTOR<int,2>(i,j)),v01=evec(VECTOR<int,2>(i,j+1)),v10=evec(VECTOR<int,2>(i+1,j)),v11=evec(VECTOR<int,2>(i+1,j+1));
            T cx0=Contour_Crossing(g00,v00,g10,v10);
            T cx1=Contour_Crossing(g01,v01,g11,v11);
            T c0x=Contour_Crossing(g00,v00,g01,v01);
            T c1x=Contour_Crossing(g10,v10,g11,v11);
            int n=(cx0>=0)+(cx1>=0)+(c0x>=0)+(c1x>=0);
            if(n<2) continue;
            VECTOR<T,2> X00((2*i-image_size)*sigma_range/image_size+1e-5,(2*j-image_size)*sigma_range/image_size);
            VECTOR<T,2> X01((2*i-image_size)*sigma_range/image_size+1e-5,(2*(j+1)-image_size)*sigma_range/image_size);
            VECTOR<T,2> X10((2*(i+1)-image_size)*sigma_range/image_size+1e-5,(2*j-image_size)*sigma_range/image_size);
            VECTOR<T,2> X11((2*(i+1)-image_size)*sigma_range/image_size+1e-5,(2*(j+1)-image_size)*sigma_range/image_size);
            if(X00.Product()>2) continue;
            VECTOR<T,2> Yx0=X00+(X10-X00)*cx0,Yx1=X01+(X11-X01)*cx1,Y0x=X00+(X01-X00)*c0x,Y1x=X10+(X11-X10)*c1x;
            if(cx0>=0 && cx1>=0){
                contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Yx0,Yx1));
                if(n==3 && c0x>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y0x,(Yx0+Yx1)/2));
                if(n==3 && c1x>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y1x,(Yx0+Yx1)/2));}
            if(c0x>=0 && c1x>=0){
                contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y0x,Y1x));
                if(n==3 && cx0>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Yx0,(Y0x+Y1x)/2));
                if(n==3 && cx1>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Yx1,(Y0x+Y1x)/2));}
            if(n>2) continue;
            if(c0x>=0 && cx0>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y0x,Yx0));
            if(c0x>=0 && cx1>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y0x,Yx1));
            if(c1x>=0 && cx0>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y1x,Yx0));
            if(c1x>=0 && cx1>=0) contour_segments.Append(VECTOR<VECTOR<T,2>,2>(Y1x,Yx1));}
}
//#####################################################################
// Function Plot_Contour_Landscape
//#####################################################################
void Plot_Contour_Landscape(int frame)
{
    char buff[1000];
    sprintf(buff, "%s/data", output_directory.c_str());
    FILE_UTILITIES::Create_Directory(buff);
    sprintf(buff, "%s/data/%03d.txt", output_directory.c_str(), frame);
    std::ofstream out(buff);

    PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    FINITE_VOLUME<TV,3>& fv=solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* icm=fv.isotropic_model;
    bool is_neo=dynamic_cast<NEO_HOOKEAN<T,3>*>(icm);
    T min=-image_size/2,max=image_size/2;
    for(int i=min;i<max;i++)
        for(int j=min;j<max;j++){
            if(is_neo && (i<0 || j<0)) continue;
            VECTOR<T,2> X(2*sigma_range*(i+.5)/image_size+1.1e-5,2*sigma_range*(j+.5)/image_size+1.2e-5);
            DIAGONAL_MATRIX<T,3> F(X.Insert(particles.X(1).y,1)),P=icm->P_From_Strain(F,1,0);
            out<<"p "<<X<<" "<<-P.To_Vector().Remove_Index(1).Normalized()<<std::endl;}

    for(int i=0;i<fv.Fe_hat.m;i++){
        contrail(i).Append(fv.Fe_hat(i).To_Vector().Remove_Index(1));
        out<<"c "<<contrail(i)<<std::endl;}

    contour_segments.Remove_All();
    image_size*=10;
    Add_Primary_Contour_Segments(*icm);
    image_size/=10;

    for(int i=0;i<contour_segments.m;i++){
        if(is_neo && (contour_segments(i).x.Min()<.2 || contour_segments(i).y.Min()<.2)) continue;
        out<<"u "<<contour_segments(i)<<std::endl;}

    out<<"t "<<particles.X.Subset(fv.strain_measure.mesh.elements.Flattened())<<std::endl;

    out<<"s "<<particles.X(1).y<<std::endl;
    out<<"m "<<icm->constant_lambda<<"  "<<icm->constant_mu<<std::endl;
}
};
}
#endif
