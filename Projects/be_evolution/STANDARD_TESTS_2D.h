//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class STANDARD_TESTS_2D
//#####################################################################
//  100. What is it
//#####################################################################
#ifndef __STANDARD_TESTS_2D__
#define __STANDARD_TESTS_2D__

#include <Tools/Images/PNG_FILE.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Log/LOG.h>
#include <Tools/Random_Numbers/RANDOM_NUMBERS.h>
#include <Geometry/Basic_Geometry/CYLINDER.h>
#include <Geometry/Constitutive_Models/STRAIN_MEASURE.h>
#include <Geometry/Geometry_Particles/DEBUG_PARTICLES.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects_Uniform/LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Implicit_Objects_Uniform/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Collisions/COLLISION_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/RIGID_BODY_BINDING.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <Solids/Forces_And_Torques/ETHER_DRAG.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
#include <Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <Dynamics/Solids_And_Fluids/SOLIDS_FLUIDS_EXAMPLE_UNIFORM.h>
#include <fstream>
namespace PhysBAM{

extern bool siggraph_hack_newton_converged;

template<class TV> class STANDARD_TESTS;
template<class T_input>
class STANDARD_TESTS<VECTOR<T_input,2> >:public SOLIDS_FLUIDS_EXAMPLE_UNIFORM<VECTOR<T_input,2> >
{
    typedef T_input T;
    typedef VECTOR<T,2> TV;typedef VECTOR<int,2> TV_INT;typedef VECTOR<T,3> TV3;
public:
    typedef SOLIDS_FLUIDS_EXAMPLE_UNIFORM<TV> BASE;
    using BASE::fluids_parameters;using BASE::solids_parameters;using BASE::output_directory;using BASE::last_frame;using BASE::frame_rate;using BASE::solid_body_collection;
    using BASE::stream_type;using BASE::solids_evolution;using BASE::parse_args;using BASE::test_number;using BASE::data_directory;using BASE::m;using BASE::s;using BASE::kg;

    std::ofstream svout;
    SOLIDS_STANDARD_TESTS<TV> tests;

    GRID<TV> mattress_grid;
    T attachment_velocity;
    bool test_forces;
    ARRAY<int> kinematic_ids;
    ARRAY<INTERPOLATION_CURVE<T,FRAME<TV> > > curves;
    INTERPOLATION_CURVE<T,T> scalar_curve;
    bool print_matrix;
    int resolution;
    T stiffness_multiplier;
    T damping_multiplier;
    ARRAY<int> externally_forced;
    ARRAY<int> constrained_particles;
    ARRAY<TV> constrained_velocities;
    TV_INT image_size;
    int rand_seed;
    bool use_rand_seed;
    RANDOM_NUMBERS<T> rand;
    bool use_residuals;
    bool use_newmark,use_newmark_be;
    bool project_nullspace;
    BACKWARD_EULER_EVOLUTION<TV>* backward_euler_evolution;
    bool use_penalty_collisions;
    bool use_constraint_collisions;
    bool no_line_search;
    bool no_descent;
    T penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length;
    bool enforce_definiteness;
    T unit_rho,unit_p,unit_N,unit_J;
    T density;
    T save_dt;
    T final_x;
    T ether_drag;
    CELL_ITERATOR<TV>* cell_iterator;
    ARRAY<TV3,TV_INT> image;
    GRID<TV> image_grid;
    bool use_vanilla_newton;
    ARRAY<INTERPOLATION_CURVE<T,TV> > kinematic_particle_positions;
    ARRAY<int> kinematic_particle_ids;

    STANDARD_TESTS(const STREAM_TYPE stream_type)
        :BASE(stream_type,0,fluids_parameters.NONE),tests(stream_type,data_directory,solid_body_collection),test_forces(false),
        print_matrix(false),resolution(0),stiffness_multiplier(1),damping_multiplier(1),image_size(500,500),rand_seed(1234),
        use_rand_seed(false),use_residuals(false),use_newmark(false),use_newmark_be(false),project_nullspace(false),
        backward_euler_evolution(new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this)),
        use_penalty_collisions(false),use_constraint_collisions(true),penalty_collisions_stiffness((T)1e4),penalty_collisions_separation((T)1e-4),
        penalty_collisions_length(1),enforce_definiteness(false),unit_rho(1),unit_p(1),unit_N(1),unit_J(1),density(pow<TV::m>(10)),final_x((T)-16.175),
        ether_drag(0),cell_iterator(0),use_vanilla_newton(false)
    {
        this->fixed_dt=1./240;
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
    void Limit_Solids_Dt(T& dt,const T time) PHYSBAM_OVERRIDE {dt=this->fixed_dt;}
    void Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Align_Deformable_Bodies_With_Rigid_Bodies() PHYSBAM_OVERRIDE {}
    void Preprocess_Solids_Substep(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Update_Solids_Parameters(const T time) PHYSBAM_OVERRIDE {}
    void Self_Collisions_Begin_Callback(const T time,const int substep) PHYSBAM_OVERRIDE {}
    void Filter_Velocities(const T dt,const T time,const bool velocity_update) PHYSBAM_OVERRIDE {}
    void Update_Time_Varying_Material_Properties(const T time) PHYSBAM_OVERRIDE {}
    // void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    // void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE {}
    void Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE {}
  //  bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id) PHYSBAM_OVERRIDE {return true;}
   // void Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id) PHYSBAM_OVERRIDE {}

//#####################################################################
// Function Register_Options
//#####################################################################
void Register_Options() PHYSBAM_OVERRIDE
{
    BASE::Register_Options();
    solids_parameters.implicit_solve_parameters.cg_projection_iterations=5;
    solids_parameters.implicit_solve_parameters.cg_iterations=1000;
    solids_parameters.implicit_solve_parameters.cg_tolerance=1e-3;
    backward_euler_evolution->newtons_method.tolerance=1;
    backward_euler_evolution->newtons_method.max_newton_step_size=1000;
    backward_euler_evolution->newtons_method.max_krylov_iterations=100;
    solids_parameters.use_rigid_deformable_contact=false;
    solids_parameters.rigid_body_collision_parameters.use_push_out=true;
    solids_parameters.triangle_collision_parameters.use_gauss_jacobi=true;
    solids_parameters.triangle_collision_parameters.repulsions_limiter_fraction=1;
    solids_parameters.triangle_collision_parameters.collisions_final_repulsion_limiter_fraction=.1;
    solids_parameters.triangle_collision_parameters.perform_self_collision=false;
    solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
    solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
    solids_parameters.triangle_collision_parameters.collisions_output_number_checked=false;
    solids_parameters.deformable_object_collision_parameters.collide_with_interior=true;
    parse_args->Add("-test_forces",&test_forces,"use fully implicit forces");
    parse_args->Add("-resolution",&resolution,"resolution","resolution used by multiple tests to change the parameters of the test");
    parse_args->Add("-stiffen",&stiffness_multiplier,"multiplier","stiffness multiplier for various tests");
    parse_args->Add("-dampen",&damping_multiplier,"multiplier","damping multiplier for various tests");
    parse_args->Add("-residuals",&use_residuals,"print residuals during timestepping");
    parse_args->Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
    parse_args->Add("-cgsolids",&solids_parameters.implicit_solve_parameters.cg_tolerance,"tolerance","CG tolerance for backward Euler");
    parse_args->Add("-use_newmark",&use_newmark,"use newmark");
    parse_args->Add("-use_newmark_be",&use_newmark_be,"use backward euler variant of newmark");
    parse_args->Add("-print_matrix",&print_matrix,"print Krylov matrix");
    parse_args->Add("-project_nullspace",&project_nullspace,"project out nullspace");
    parse_args->Add("-projection_iterations",&solids_parameters.implicit_solve_parameters.cg_projection_iterations,"iterations","number of iterations used for projection in cg");
    parse_args->Add("-seed",&rand_seed,&use_rand_seed,"seed","random seed to use");
    parse_args->Add("-test_system",&solids_parameters.implicit_solve_parameters.test_system,"test system");
    parse_args->Add("-ether_drag",&ether_drag,"drag","Ether drag");
    parse_args->Add("-image_size",&image_size,"size","image size for plots");
    parse_args->Add("-kry_it",&backward_euler_evolution->newtons_method.max_krylov_iterations,"iter","maximum iterations for Krylov solver");
    parse_args->Add("-kry_tol",&backward_euler_evolution->newtons_method.krylov_tolerance,"tol","tolerance for Krylov solver");
    parse_args->Add("-newton_it",&backward_euler_evolution->newtons_method.max_iterations,"iter","maximum iterations for Newton");
    parse_args->Add("-newton_tol",&backward_euler_evolution->newtons_method.tolerance,"tol","tolerance for Newton");
    parse_args->Add("-newton_cd_tol",&backward_euler_evolution->newtons_method.countdown_tolerance,"tol","tolerance for Newton");
    parse_args->Add("-newton_max_step",&backward_euler_evolution->newtons_method.max_newton_step_size,"size","Limit newton step to this size");
    parse_args->Add("-no_descent",&no_descent,"Don't ensure descent direction");
    parse_args->Add("-debug_newton",&backward_euler_evolution->newtons_method.debug,"Enable diagnostics in Newton's method");
    parse_args->Add("-kry_fail",&backward_euler_evolution->newtons_method.fail_on_krylov_not_converged,"terminate if Krylov solver fails to converge");
    parse_args->Add("-angle_tol",&backward_euler_evolution->newtons_method.angle_tolerance,"tol","gradient descent tolerance");
    parse_args->Add("-final_x",&final_x,"position","final x position");
    parse_args->Add_Not("-mr",&backward_euler_evolution->newtons_method.use_cg,"use minres instead of cg");
    parse_args->Add("-no_line_search",&no_line_search,"disable line search");
    parse_args->Add("-gss",&backward_euler_evolution->newtons_method.use_golden_section_search,"use golden section search instead of wolfe conditions line search");
    parse_args->Add("-backtrack",&backward_euler_evolution->newtons_method.use_backtracking,"use backtracking line search instead of wolfe conditions line search");
    parse_args->Add("-use_penalty",&use_penalty_collisions,"use penalty collisions");
    parse_args->Add_Not("-no_constraints",&use_constraint_collisions,"disable constrained optimization for collisions");
    parse_args->Add("-penalty_stiffness",&penalty_collisions_stiffness,"tol","penalty collisions stiffness");
    parse_args->Add("-penalty_separation",&penalty_collisions_separation,"tol","penalty collisions separation");
    parse_args->Add("-penalty_length",&penalty_collisions_length,"tol","penalty collisions length scale");
    parse_args->Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");
    parse_args->Add("-use_tri_col",&solids_parameters.triangle_collision_parameters.perform_self_collision,"use triangle collisions");
    parse_args->Add("-use_vanilla_newton",&use_vanilla_newton,"use triangle collisions");
}
//#####################################################################
// Function Parse_Options
//#####################################################################
void Parse_Options() PHYSBAM_OVERRIDE
{
    BASE::Parse_Options();
    tests.data_directory=data_directory;
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    output_directory=STRING_UTILITIES::string_sprintf("Test_%d",test_number);
    if(use_rand_seed) rand.Set_Seed(rand_seed);
    solids_parameters.implicit_solve_parameters.project_nullspace_frequency=project_nullspace;
    if(use_newmark || use_newmark_be) backward_euler_evolution=0;
    else{delete solids_evolution;solids_evolution=backward_euler_evolution;}
    if(backward_euler_evolution && backward_euler_evolution->newtons_method.use_golden_section_search)
        backward_euler_evolution->newtons_method.use_wolfe_search=false;
    if(backward_euler_evolution && backward_euler_evolution->newtons_method.use_backtracking)
        backward_euler_evolution->newtons_method.use_wolfe_search=false;
    if(backward_euler_evolution && no_line_search)
        backward_euler_evolution->newtons_method.use_wolfe_search=false;
    if(backward_euler_evolution && no_descent)
        backward_euler_evolution->newtons_method.use_gradient_descent_failsafe=false;
    if(use_vanilla_newton) backward_euler_evolution->newtons_method.Make_Vanilla_Newton();

    unit_rho=kg/pow<TV::m>(m);
    unit_N=kg*m/(s*s);
    unit_p=unit_N/(m*m);
    unit_J=unit_N*m;
    penalty_collisions_length*=m;
    penalty_collisions_separation*=m;
    ether_drag/=s;
    penalty_collisions_stiffness*=unit_J;
    density*=unit_rho;
    if(backward_euler_evolution){
        backward_euler_evolution->newtons_method.tolerance*=unit_N*s;
        backward_euler_evolution->newtons_method.krylov_tolerance/=sqrt(unit_N*s);
        backward_euler_evolution->minimization_objective.collision_thickness*=m;}
    solids_parameters.use_trapezoidal_rule_for_velocities=!use_newmark_be;

    switch(test_number){
        case 1:
            if(!resolution) resolution=10;
            mattress_grid=GRID<TV>(TV_INT()+resolution+1,RANGE<TV>::Centered_Box());
            break;
        case 100:
            mattress_grid=GRID<TV>(TV_INT(3,2),RANGE<TV>(TV(),TV(2,1)));
            image_grid.Initialize(image_size,RANGE<TV>::Centered_Box()*5,true);
            cell_iterator=new CELL_ITERATOR<TV>(image_grid);
            image.Resize(image_grid.Cell_Indices());
            this->fixed_dt=1./24;
            frame_rate=1/(this->fixed_dt*image_size.Product());
            last_frame=1;
            for(int i=0;i<9;i++) if(i!=4) constrained_particles.Append(i);
            constrained_velocities.Resize(constrained_particles.m,true,true,TV());
            break;
        case 101:
            mattress_grid=GRID<TV>(TV_INT(3,2),RANGE<TV>(TV(),TV(2,1)));
            kinematic_particle_ids.Append(0);
            kinematic_particle_ids.Append(2);
            kinematic_particle_ids.Append(3);
            kinematic_particle_ids.Append(4);
            kinematic_particle_ids.Append(5);
            kinematic_particle_positions.Resize(5);
            kinematic_particle_positions(0).Add_Control_Point(0,TV(0,0));
            kinematic_particle_positions(1).Add_Control_Point(0,TV(2,0));
            kinematic_particle_positions(2).Add_Control_Point(0,TV(0,1));
            kinematic_particle_positions(3).Add_Control_Point(0,TV(1,1));
            kinematic_particle_positions(4).Add_Control_Point(0,TV(2,1));
            kinematic_particle_positions(0).Add_Control_Point(1,TV(0,0));
            kinematic_particle_positions(1).Add_Control_Point(1,TV(final_x,0));
            kinematic_particle_positions(2).Add_Control_Point(1,TV(0,1));
            kinematic_particle_positions(3).Add_Control_Point(1,TV(final_x/2,1));
            kinematic_particle_positions(4).Add_Control_Point(1,TV(final_x,1));
            break;
    }

    if(use_constraint_collisions) use_penalty_collisions=false;
    if(use_penalty_collisions || use_constraint_collisions){
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
        solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;}

    solid_body_collection.Print_Residuals(use_residuals);
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
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

    switch(test_number){
        case 1:{
            RIGID_BODY_STATE<TV> initial_state(FRAME<TV>(TV(0,4)*m));
            tests.Create_Mattress(mattress_grid,true,&initial_state,density);
            tests.Add_Ground();
            break;}
        case 100:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            tests.Add_Ground();
            break;}
        case 101:{
            tests.Create_Mattress(mattress_grid,true,0,density);
            tests.Add_Ground();
            break;}
        default:
            LOG::cerr<<"Initial Data: Unrecognized test number "<<test_number<<std::endl;exit(1);}

    // add structures and rigid bodies to collisions
    if(automatically_add_to_collision_structures) deformable_body_collection.collisions.collision_structures.Append_Elements(deformable_body_collection.structures);
    solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append_Elements(deformable_body_collection.structures);

    // correct number nodes
    for(int i=0;i<deformable_body_collection.structures.m;i++) deformable_body_collection.structures(i)->Update_Number_Nodes();

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
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    Get_Initial_Data();

    switch(test_number){
        case 1:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Gravity();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        case 100:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        case 101:{
            TRIANGULATED_AREA<T>& ta=deformable_body_collection.template Find_Structure<TRIANGULATED_AREA<T>&>();
            Add_Constitutive_Model(ta,(T)1e5*unit_p,(T).45,(T)0*s);
            break;}
        default:
            LOG::cerr<<"Missing bodies implementation for test number "<<test_number<<std::endl;exit(1);}

    if(ether_drag) solid_body_collection.Add_Force(new ELASTIC_ETHER_DRAG<TV>(deformable_body_collection.particles,true,ether_drag,1,save_dt));

    if(use_penalty_collisions)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object))
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);
            solid_body_collection.Add_Force(new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,
                    iot,penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length));}
    else if(use_constraint_collisions && backward_euler_evolution)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object)){
                lio->levelset.interpolation=new CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> >();
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);}
            backward_euler_evolution->minimization_objective.collision_objects.Append(iot);
            backward_euler_evolution->minimization_objective.coefficient_of_friction.Append(rigid_body_collection.Rigid_Body(b).coefficient_of_friction);}
    else
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision)
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    if(enforce_definiteness) solid_body_collection.Enforce_Definiteness(true);
    if(backward_euler_evolution) backward_euler_evolution->minimization_objective.Disable_Current_Colliding_Pairs(0);

    for(int i=0;i<solid_body_collection.solids_forces.m;i++) solid_body_collection.solids_forces(i)->use_implicit_velocity_independent_forces=true;
    for(int i=0;i<solid_body_collection.rigid_body_collection.rigids_forces.m;i++)
        solid_body_collection.rigid_body_collection.rigids_forces(i)->use_implicit_velocity_independent_forces=true;
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++) deformable_body_collection.deformables_forces(i)->use_implicit_velocity_independent_forces=true;
}

//#####################################################################
// Function Set_External_Velocities
//#####################################################################
void Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    V.Subset(constrained_particles)=constrained_velocities;
    for(int i=0;i<kinematic_particle_ids.m;i++)
        V(kinematic_particle_ids(i))=kinematic_particle_positions(i).Derivative(velocity_time);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
void Set_External_Positions(ARRAY_VIEW<TV> X,const T time) PHYSBAM_OVERRIDE
{
    for(int i=0;i<kinematic_particle_ids.m;i++)
        X(kinematic_particle_ids(i))=kinematic_particle_positions(i).Value(time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
void Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time) PHYSBAM_OVERRIDE
{
    V.Subset(constrained_particles).Fill(TV());
    V.Subset(externally_forced).Fill(TV());
    V.Subset(kinematic_particle_ids).Fill(TV());
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
    for(int i=0;i<kinematic_ids.m;i++)
        if(id==kinematic_ids(i)){
            frame=curves(i).Value(time);
            break;}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
bool Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    for(int i=0;i<kinematic_ids.m;i++)
        if(id==kinematic_ids(i)){
            twist=curves(i).Derivative(time);
            return true;}
    return false;
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
void Preprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    save_dt=dt;
    if(test_forces){
        solid_body_collection.deformable_body_collection.Test_Energy(time);
        solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}

    if(test_number==100) particles.X(4)=cell_iterator->Location();
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
void Postprocess_Substep(const T dt,const T time) PHYSBAM_OVERRIDE
{
    if(test_number==100){
        if(siggraph_hack_newton_converged) image(cell_iterator->index)=TV3(0,1,0);
        else image(cell_iterator->index)=TV3(1,0,0);
        cell_iterator->Next();
        if(!cell_iterator->Valid()) PNG_FILE<T>::Write("image.png",image);}
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
void Preprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
void Postprocess_Frame(const int frame) PHYSBAM_OVERRIDE
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),TV3(1,0,0));
    for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),TV3(0,1,0));
    for(int i=0;i<kinematic_particle_ids.m;i++) Add_Debug_Particle(particles.X(kinematic_particle_ids(i)),TV3(1,1,0));
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
void Add_Constitutive_Model(TRIANGULATED_AREA<T>& ta,T stiffness,T poissons_ratio,T damping)
{
    solid_body_collection.Add_Force(Create_Finite_Volume(ta,new COROTATED_FIXED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier)));

    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(damping*damping_multiplier){
        DEFORMABLES_FORCES<TV>* force=Create_Finite_Volume(ta,new COROTATED_FIXED<T,2>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier));
        force->use_implicit_velocity_independent_forces=true;
        force->Update_Position_Based_State(0,true);
        solid_body_collection.Add_Force(new RALEIGH_DAMPING_FORCE<TV>(particles,force,damping*damping_multiplier,1,save_dt));}
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
GRAVITY<TV>& Add_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true);
    solid_body_collection.Add_Force(g);
    return *g;
}
};
}
#endif
