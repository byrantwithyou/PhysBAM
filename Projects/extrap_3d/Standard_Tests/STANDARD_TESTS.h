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
//   38. Many jollos falling on the ground
//   39. One jello falling on a stationary one
//   40. Two jellos running on each other
//   41. Bunch of jellos rolling towards camera
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Collisions/COLLISION_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_SMOOTH.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/SVK_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D2.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/MOONEY_RIVLIN_3D_EXTRAPOLATED.h>
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
    typedef VECTOR<T,3> TV;
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
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_extended_neohookean_smooth;
    bool use_extended_svk;
    bool use_corotated;
    bool use_mooney_rivlin,use_extended_mooney_rivlin;
    bool use_corot_blend;
    bool dump_sv;
    int kinematic_id,kinematic_id2,kinematic_id3;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve,curve2,curve3;
    bool print_matrix;
    int parameter;
    T stiffness_multiplier;
    T damping_multiplier;
    bool use_constant_ife;
    bool forces_are_removed;
    
    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),
        use_extended_neohookean_refined(false),use_extended_neohookean_hyperbola(false),use_extended_neohookean_smooth(false),use_extended_svk(false),
        use_corotated(false),use_corot_blend(false),dump_sv(false),print_matrix(false),use_constant_ife(false)
    {
    }

    virtual ~STANDARD_TESTS()
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
    void Set_External_Positions(ARRAY_VIEW<TV> X,ARRAY_VIEW<ROTATION<TV> > rotation,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    //void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {if(dump_sv)svout.close();} //modified 22 Nov 2011
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
    parse_args->Add_Option_Argument("-use_ext_neo");
    parse_args->Add_Option_Argument("-use_ext_neo_ref");
    parse_args->Add_Option_Argument("-use_ext_neo_hyper");
    parse_args->Add_Option_Argument("-use_ext_neo_smooth");
    parse_args->Add_Option_Argument("-use_ext_svk");
    parse_args->Add_Option_Argument("-use_ext_mooney");
    parse_args->Add_Option_Argument("-use_mooney");
    parse_args->Add_Option_Argument("-use_corotated");
    parse_args->Add_Option_Argument("-use_corot_blend");
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

    switch(test_number){
        case 17: case 18: case 24: case 25: case 27: case 10: case 11:
            mattress_grid=GRID<TV>(10,10,10,(T)-1,(T)1,(T)-1,(T)1,(T)-1,(T)1);
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
        case 35: case 36:
            mattress_grid1=GRID<TV>(10,10,10,(T)-1.0,(T)1.0,(T)-1.0,(T)1.0,(T)-1.0,(T)1.0);
            mattress_grid2=GRID<TV>(12,12,12,(T)-1.2,(T)1.2,(T)-1.2,(T)1.2,(T)-1.2,(T)1.2);
            mattress_grid3=GRID<TV>(15,15,15,(T)-1.5,(T)1.5,(T)-1.5,(T)1.5,(T)-1.5,(T)1.5);
            break;
        case 37: case 38: case 39: case 40: case 41:
            mattress_grid=GRID<TV>(40,40,40,(T)-1.0,(T)1.0,(T)-1.0,(T)1.0,(T)-1.0,(T)1.0);
            break;
    	default:
            mattress_grid=GRID<TV>(20,10,20,(T)-1,(T)1,(T)-.5,(T).5,(T)-1,(T)1);
    }

    print_matrix=parse_args->Is_Value_Set("-print_matrix");
    parameter=parse_args->Get_Integer_Value("-parameter");
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
    test_forces=parse_args->Is_Value_Set("-test_forces");
    use_extended_neohookean=parse_args->Is_Value_Set("-use_ext_neo");
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_ref");
    use_extended_neohookean_hyperbola=parse_args->Is_Value_Set("-use_ext_neo_hyper");
    use_extended_mooney_rivlin=parse_args->Is_Value_Set("-use_ext_mooney");
    use_mooney_rivlin=parse_args->Is_Value_Set("-use_mooney");
    use_extended_neohookean_smooth=parse_args->Is_Value_Set("-use_ext_neo_smooth");
    use_extended_svk=parse_args->Is_Value_Set("-use_ext_svk");
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    use_constant_ife=parse_args->Get_Option_Value("-use_constant_ife");
    solids_parameters.implicit_solve_parameters.test_system=parse_args->Is_Value_Set("-test_system");
    
    semi_implicit=parse_args->Is_Value_Set("-semi_implicit");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");

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
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
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
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            frame_rate=120;
            last_frame=10*120;
            last_frame=1000;
            break;
        case 34:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            last_frame = 400;
            break;
        case 35: case 36:
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
        case 41:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.triangle_collision_parameters.perform_self_collision=true;
            frame_rate=60;
            last_frame=500;
            break;
        case 24:
        case 25:
        case 26:
        case 27:
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
        case 29: case 30:
            solids_parameters.cfl=(T)5;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            frame_rate=40;
            break;
        case 31:
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            last_frame=1500;
            frame_rate=60;
            break;
        case 5:
        case 6:
            frame_rate=60;
            last_frame=(int)(3*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            break;
        case 48:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_iterations=100000;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;//(int)(200*frame_rate);
            solids_parameters.cfl=(T)1;
            solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure=false;
            break;    
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

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
        case 33:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_110K.tet",
                RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,2.2,0),ROTATION<TV>(T(pi/2),TV(0,0,0)))),true,true,density,0.008);
            
            
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

            curve3.Add_Control_Point(0,FRAME<TV>(TV(0,2.75,0),ROTATION<TV>(0,TV(0,0,1))));
            curve3.Add_Control_Point(0.5,FRAME<TV>(TV(0,2.75,0),ROTATION<TV>(0,TV(0,0,1))));
            curve3.Add_Control_Point(2,FRAME<TV>(TV(0,2,0),ROTATION<TV>(0,TV(0,0,1))));

            for (int i=0; i<60; i++){
                curve.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV(-(T).4,1.5,-.75),ROTATION<TV>(-i,TV(0,0,1))));
                curve2.Add_Control_Point(i/angular_velocity,FRAME<TV>(TV((T).4,1.5,-.75),ROTATION<TV>(i,TV(0,0,1))));}

            tests.Add_Ground();
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
            
            for (int i=0; i<64; i++) curve.Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(x_stop,0,0),ROTATION<TV>((T)pi/2.0*i,TV(1,0,0))));
          
            curve2.Add_Control_Point(0,FRAME<TV>(TV(-x_start,0,0),ROTATION<TV>((T)0,TV(1,0,0))));
            curve2.Add_Control_Point(t_stop,FRAME<TV>(TV(-x_stop,0,0),ROTATION<TV>((T)0,TV(1,0,0))));
            
            for (int i=0; i<64; i++) curve2.Add_Control_Point(t_rot+i*2,FRAME<TV>(TV(-x_stop,0,0),ROTATION<TV>((T)-pi/2.0*i,TV(1,0,0))));
            
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
            tests.Add_Ground(1);
            break;
        }
        case 39:
        {
            RIGID_BODY_STATE<TV> initial_state1(FRAME<TV>(TV(1,8,1.5),ROTATION<TV>(T(pi/4),TV(1.3,0.3,0.7)))); 
            RIGID_BODY_STATE<TV> initial_state2(FRAME<TV>(TV(0,1,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state1);
            tests.Create_Mattress(mattress_grid,true,&initial_state2);
            tests.Add_Ground(1);
            break;
        }
        case 3:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 4: case 29:{
//            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/sphere.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_4K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)6,0))),true,true,density,.085);
            tests.Add_Ground();
            last_frame=2000;
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
        case 26:
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
            RIGID_BODY<TV>& box1=tests.Add_Analytic_Box(TV(3,20,20));
            box1.X()=TV(10,10,10);
            box1.is_static=true;
            tests.Add_Ground();
            last_frame=250;
            break;}        
        case 31: {
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/buddha.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV((T)0,(T)4.95,(T)0))),true,true,density,5);
            RIGID_BODY<TV>& sphere=tests.Add_Analytic_Sphere(1.2,density,5);
            sphere.is_static=false;
            tests.Add_Ground();            
            kinematic_id=sphere.particle_index;
            rigid_body_collection.rigid_body_particle.kinematic(sphere.particle_index)=true;
            curve.Add_Control_Point(0,FRAME<TV>(TV(0,4.8,-50)));
            curve.Add_Control_Point(40,FRAME<TV>(TV(0,4.8,350)));
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
        default:
            LOG::cerr<<"Unrecognized test number "<<test_number<<std::endl;exit(1);}

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
        case 33:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,1e5,0.4,0.05);
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
        case 39:{
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
            }
            break;}
        case 7:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            for(int i=1; i<=deformable_body_collection.particles.X.m; i++) deformable_body_collection.particles.X(i).y=3;
            break;}
        case 4: {
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01);
            break;}
        case 29:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)0,(T)0,(T).01);
            forces_are_removed=true;
            break;}
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
        case 24:
        case 25:
        case 26:
        case 27:
        case 28:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            break;} 
        case 31:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)5e5,(T).45,(T).01);
            solid_body_collection.template Find_Force<GRAVITY<TV>&>().gravity=0.5;
            break;} 
        case 30:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01);
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            for(int i=1; i<=deformable_body_collection.particles.X.m; i++){ deformable_body_collection.particles.V(i).x=(T)100;deformable_body_collection.particles.V(i).y=9.8;deformable_body_collection.particles.V(i).z=0;}
            break;} 
        case 9:{break;}
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

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
    if(test_number==28){
        int m=mattress_grid.counts.x;
        int n=mattress_grid.counts.y;
        int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        TV velocity_y = velocity_time<final_time?TV(0,attachment_velocity,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=-velocity_x;V(m+m*(j-1)+m*n*(ij-1))=velocity_x;}
    }
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
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(id==kinematic_id) twist=curve.Derivative(time);
    if(id==kinematic_id2) twist=curve2.Derivative(time);
    if(id==kinematic_id3) twist=curve3.Derivative(time);
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
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
void Update_Time_Varying_Material_Properties(const T time)
{   if(test_number==29 && time > .1){
        T critical=(T)3.0;
        T critical2=(T)3.5;
        T start_young=(T)4; T end_young=(T)5;
        if(time>critical && time<critical2) {
            
            DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
            FINITE_VOLUME<TV,3>& fv = deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
            CONSTITUTIVE_MODEL<T,3>& icm = fv.constitutive_model;
            T young = pow(10.0,start_young + (time-critical)/(critical2-critical)*(end_young-start_young));
            icm.Update_Lame_Constants(young,(T).45,(T).01);
            forces_are_removed=false;
        }    
    }
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################    
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if (dump_sv)
    {
        FINITE_VOLUME<TV,3>& force_field = solid_body_collection.deformable_body_collection.template Find_Force<FINITE_VOLUME<TV,3>&>();
        ARRAY<DIAGONAL_MATRIX<T,3> >& sv = force_field.Fe_hat;
        
        for (int i=1; i<=sv.m; i++)
        {
            svout << sv(i).x11 << " " << sv(i).x22 << " " << sv(i).x33 << std::endl;
        }
    }
} 
//#####################################################################
// Function Bind_Intersecting_Particles
//#####################################################################
void Bind_Intersecting_Particles()
{
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
    
    if (dump_sv)
    {
        std::string output_file = STRING_UTILITIES::string_sprintf("Standard_Tests/Test_%d/SV_%d",test_number,frame);
        svout.open(output_file.c_str());
    }

    if(test_number==32 && frame==1100) Bind_Intersecting_Particles();
}    
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,T stiffness,T poissons_ratio,T damping, T cutoff = 0.4, T efc = 20)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* icm=0;

    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc*stiffness*stiffness_multiplier);
    else if(use_extended_svk) icm=new SVK_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,efc*stiffness*stiffness_multiplier);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,cutoff,.6,efc*stiffness*stiffness_multiplier);
    else if(use_extended_mooney_rivlin) icm=new MOONEY_RIVLIN_3D_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20*stiffness*stiffness_multiplier);
    else if(use_mooney_rivlin) icm=new MOONEY_RIVLIN_3D2<T>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20*stiffness*stiffness_multiplier);
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
};
}
#endif
