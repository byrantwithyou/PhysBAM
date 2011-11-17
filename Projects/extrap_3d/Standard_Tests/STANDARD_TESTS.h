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
//    6. Precursor to a smash test
//    8. Falling mattress
//   17. Matress, no gravity, random start
//   18. Matress, no gravity, point start
//   24. Big 6 sides stretch
//   25. Big 8 corners stretch
//   26. Big stretch/bend
//   27. Force inversion
//#####################################################################
#ifndef __STANDARD_TESTS__
#define __STANDARD_TESTS__

#include <PhysBAM_Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <PhysBAM_Geometry/Solids_Geometry/DEFORMABLE_GEOMETRY_COLLECTION.h>
#include <PhysBAM_Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/COROTATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_COROTATED_BLEND.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_REFINED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/NEO_HOOKEAN_EXTRAPOLATED_HYPERBOLA.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Constitutive_Models/ROTATED_LINEAR.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/GRAVITY.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <PhysBAM_Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
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

    SOLIDS_STANDARD_TESTS<TV> tests;
    
    GRID<TV> mattress_grid;
    T attachment_velocity;
    bool semi_implicit;
    bool test_forces;
    bool use_extended_neohookean;
    bool use_extended_neohookean_refined;
    bool use_extended_neohookean_hyperbola;
    bool use_corotated;
    bool use_corot_blend;
    bool dump_sv;
    int kinematic_id;
    INTERPOLATION_CURVE<T,FRAME<TV> > curve;
    bool print_matrix;
    int parameter;
    T stiffness_multiplier;
    T damping_multiplier;
    
    // Came from Test 48
    int number_side_panels;
    T aspect_ratio,side_length;
    T bending_stiffness_multiplier,bending_damping_multiplier,planar_damping_multiplier,axial_bending_stiffness_multiplier,axial_bending_damping_multiplier;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(*this,solid_body_collection),semi_implicit(false),test_forces(false),use_extended_neohookean(false),
        use_extended_neohookean_refined(false),use_corotated(false),use_corot_blend(false),dump_sv(false),print_matrix(false)
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
    void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Set_Particle_Is_Simulated(ARRAY<bool>& particle_is_simulated) PHYSBAM_OVERRIDE {}
    void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
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
        case 17: case 18: case 24: case 25: case 27:
            mattress_grid=GRID<TV>(10,10,10,(T)-1,(T)1,(T)-1,(T)1,(T)-1,(T)1);
            break;
        case 26:
            mattress_grid=GRID<TV>(40,5,5,(T)-4,(T)4,(T)-.5,(T).5,(T)-.5,(T).5);
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
    use_extended_neohookean_refined=parse_args->Is_Value_Set("-use_ext_neo_hyper");
    use_corotated=parse_args->Is_Value_Set("-use_corotated");
    use_corot_blend=parse_args->Is_Value_Set("-use_corot_blend");
    dump_sv=parse_args->Is_Value_Set("-dump_sv");
    
    semi_implicit=parse_args->Is_Value_Set("-semi_implicit");
    if(parse_args->Is_Value_Set("-project_nullspace")) solids_parameters.implicit_solve_parameters.project_nullspace_frequency=1;
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=parse_args->Get_Integer_Value("-projection_iterations");

    switch(test_number){
        case 1:
        case 2:
        case 3:
        case 4:
        case 8:
        case 17:
        case 18:
            solids_parameters.cfl=(T)5;
            break;
        case 24:
        case 25:
        case 26:
        case 27:
            attachment_velocity = 0.2;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            solids_parameters.implicit_solve_parameters.cg_iterations=900;
            solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;
            last_frame=1000;
            break;
        case 5:
        case 6:
            frame_rate=60;
            last_frame=(int)(3*frame_rate);
            solids_parameters.cfl=(T)5.9;
            solids_parameters.implicit_solve_parameters.cg_iterations=300;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T)1e-3;
            break;
        case 48:
            frame_rate=24;
            solids_parameters.implicit_solve_parameters.cg_iterations=600;
            solids_parameters.implicit_solve_parameters.cg_tolerance=(T).01;
            last_frame=200;//(int)(200*frame_rate);
            aspect_ratio=1;
            side_length=(T)1;
            number_side_panels=2;
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
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/maggot_8K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)3,0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 4:{
            tests.Create_Tetrahedralized_Volume(data_directory+"/Tetrahedralized_Volumes/armadillo_20K.tet",RIGID_BODY_STATE<TV>(FRAME<TV>(TV(0,(T)80,0))),true,true,density);
            tests.Add_Ground();
            break;}
        case 8:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4,0)));
            tests.Create_Mattress(mattress_grid,true,&initial_state);
            tests.Add_Ground();
            break;}
        case 17:
        case 18:
        case 24:
        case 25:
        case 26:
        case 27:{
            tests.Create_Mattress(mattress_grid,true,0);
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
        case 5:
        case 6:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            break;}
        case 4:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            solid_body_collection.Add_Force(new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true));
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e6,(T).45,(T).01);
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
        case 27:{
            TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume=deformable_body_collection.deformable_geometry.template Find_Structure<TETRAHEDRALIZED_VOLUME<T>&>();
            Add_Constitutive_Model(tetrahedralized_volume,(T)1e5,(T).45,(T).01);
            break;}         
        default:
            LOG::cerr<<"Missing implementation for test number "<<test_number<<std::endl;exit(1);}

    if(solid_body_collection.deformable_body_collection.mpi_solids){
        VECTOR<int,3> processes_per_dimension(2,1,1);
        deformable_body_collection.mpi_solids->Simple_Partition(deformable_body_collection,solid_body_collection.rigid_body_collection.rigid_geometry_collection,particles.X,processes_per_dimension);}

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
        int m=mattress_grid.counts.x;
	int n=mattress_grid.counts.y;
	int mn=mattress_grid.counts.z;
        TV velocity_x = velocity_time<final_time?TV(attachment_velocity,0,0):TV();
        for(int ij=1;ij<=mn;ij++)for(int j=1;j<=n;j++){V(1+m*(j-1)+m*n*(ij-1))=velocity_x;V(m+m*(j-1)+m*n*(ij-1))=-velocity_x;}
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
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    if(id==kinematic_id) frame=curve.Value(time);
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    if(id==kinematic_id) twist=curve.Derivative(time);
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
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TETRAHEDRALIZED_VOLUME<T>& tetrahedralized_volume,T stiffness,T poissons_ratio,T damping)
{
    ISOTROPIC_CONSTITUTIVE_MODEL<T,3>* icm=0;
    if(use_extended_neohookean) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,20*stiffness*stiffness_multiplier);
    else if(use_extended_neohookean_refined) icm=new NEO_HOOKEAN_EXTRAPOLATED_REFINED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.4,.6,20*stiffness*stiffness_multiplier);
    if(use_extended_neohookean_hyperbola) icm=new NEO_HOOKEAN_EXTRAPOLATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier,.1);
    else if(use_corotated) icm=new COROTATED<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else if(use_corot_blend) icm=new NEO_HOOKEAN_COROTATED_BLEND<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    else icm=new NEO_HOOKEAN<T,3>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier);
    solid_body_collection.Add_Force(Create_Finite_Volume(tetrahedralized_volume,icm));
}
//#####################################################################
};
}
#endif
