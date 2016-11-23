//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Tools/Interpolation/INTERPOLATION_CURVE.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <Geometry/Implicit_Objects/SMOOTH_LEVELSET_IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Constitutive_Models/COROTATED_FIXED.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/ELASTIC_ETHER_DRAG.h>
#include <Deformables/Forces/FINITE_VOLUME.h>
#include <Deformables/Forces/RALEIGH_DAMPING_FORCE.h>
#include <Solids/Examples_And_Drivers/SOLIDS_EXAMPLE.h>
#include <Solids/Forces_And_Torques/GRAVITY.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include <Solids/Standard_Tests/SOLIDS_STANDARD_TESTS.h>
#include <fstream>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),test_forces(false),print_matrix(false),
    resolution(0),stiffness_multiplier(1),thickness_multiplier(1),curvature_stiffness_multiplier(1),damping_multiplier(1),input_poissons_ratio(-1),input_youngs_modulus(0),
    input_friction(.3),ether_drag(0),rand_seed(1234),
    use_rand_seed(false),use_newmark(false),use_newmark_be(false),project_nullspace(false),
    backward_euler_evolution(new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this)),
    use_penalty_collisions(false),use_constraint_collisions(true),penalty_collisions_stiffness((T)1e4),
    penalty_collisions_separation((T)1e-4),penalty_collisions_length(1),enforce_definiteness(false),
    unit_rho(1),unit_p(1),unit_N(1),unit_J(1),density(pow<TV::m>(10)),use_penalty_self_collisions(true),
    self_collide_surface_only(false),use_vanilla_newton(false),gauss_order(3),threads(1)
{
    this->fixed_dt=1./240;
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
    parse_args.Add("-test_forces",&test_forces,"use fully implicit forces");
    parse_args.Add("-print_matrix",&print_matrix,"print Krylov matrix");
    parse_args.Add("-resolution",&resolution,"resolution","resolution used by multiple tests to change the parameters of the test");
    parse_args.Add("-stiffen",&stiffness_multiplier,"multiplier","stiffness multiplier");
    parse_args.Add("-thicken",&thickness_multiplier,"multiplier","thickness multiplier");
    parse_args.Add("-kappa_stiffen",&curvature_stiffness_multiplier,"multiplier","stiffness multiplier for curvature");
    parse_args.Add("-dampen",&damping_multiplier,"multiplier","damping multiplier");
    parse_args.Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
    parse_args.Add("-cgsolids",&solids_parameters.implicit_solve_parameters.cg_tolerance,"tolerance","CG tolerance for backward Euler");
    parse_args.Add("-use_newmark",&use_newmark,"use newmark");
    parse_args.Add("-use_newmark_be",&use_newmark_be,"use backward euler variant of newmark");
    parse_args.Add("-project_nullspace",&project_nullspace,"project out nullspace");
    parse_args.Add("-projection_iterations",&solids_parameters.implicit_solve_parameters.cg_projection_iterations,"iterations","number of iterations used for projection in cg");
    parse_args.Add("-seed",&rand_seed,&use_rand_seed,"seed","random seed to use");
    parse_args.Add("-test_system",&solids_parameters.implicit_solve_parameters.test_system,"test system");
    parse_args.Add("-poissons_ratio",&input_poissons_ratio,"ratio","poissons_ratio");
    parse_args.Add("-youngs_modulus",&input_youngs_modulus,"stiffness","youngs modulus");
    parse_args.Add("-friction",&input_friction,"friction","amount of friction");
    parse_args.Add("-ether_drag",&ether_drag,"drag","Ether drag");
    parse_args.Add("-kry_it",&backward_euler_evolution->newtons_method.max_krylov_iterations,"iter","maximum iterations for Krylov solver");
    parse_args.Add("-kry_tol",&backward_euler_evolution->newtons_method.krylov_tolerance,"tol","tolerance for Krylov solver");
    parse_args.Add("-newton_it",&backward_euler_evolution->newtons_method.max_iterations,"iter","maximum iterations for Newton");
    parse_args.Add("-newton_tol",&backward_euler_evolution->newtons_method.tolerance,"tol","tolerance for Newton");
    parse_args.Add("-newton_cd_tol",&backward_euler_evolution->newtons_method.countdown_tolerance,"tol","tolerance for Newton");
    parse_args.Add("-newton_max_step",&backward_euler_evolution->newtons_method.max_newton_step_size,"size","Limit newton step to this size");
    parse_args.Add("-no_descent",&no_descent,"Don't ensure descent direction");
    parse_args.Add("-debug_newton",&backward_euler_evolution->newtons_method.debug,"Enable diagnostics in Newton's method");
    parse_args.Add("-kry_fail",&backward_euler_evolution->newtons_method.fail_on_krylov_not_converged,"terminate if Krylov solver fails to converge");
    parse_args.Add("-newton_fail",&backward_euler_evolution->fail_on_newton_not_converged,"terminate if Newton solver fails to converge");
    parse_args.Add("-angle_tol",&backward_euler_evolution->newtons_method.angle_tolerance,"tol","gradient descent tolerance");
    parse_args.Add_Not("-mr",&backward_euler_evolution->newtons_method.use_cg,"use minres instead of cg");
    parse_args.Add("-no_line_search",&no_line_search,"disable line search");
    parse_args.Add("-gss",&backward_euler_evolution->newtons_method.use_golden_section_search,"use golden section search instead of wolfe conditions line search");
    parse_args.Add("-backtrack",&backward_euler_evolution->newtons_method.use_backtracking,"use backtracking line search instead of wolfe conditions line search");
    parse_args.Add("-use_penalty",&use_penalty_collisions,"use penalty collisions");
    parse_args.Add_Not("-no_constraints",&use_constraint_collisions,"disable constrained optimization for collisions");
    parse_args.Add_Not("-no_collisions_in_solve",&backward_euler_evolution->minimization_objective.collisions_in_solve,"disable collisions in solve");
    parse_args.Add("-penalty_stiffness",&penalty_collisions_stiffness,"tol","penalty collisions stiffness");
    parse_args.Add("-penalty_separation",&penalty_collisions_separation,"tol","penalty collisions separation");
    parse_args.Add("-penalty_length",&penalty_collisions_length,"tol","penalty collisions length scale");
    parse_args.Add("-enf_def",&enforce_definiteness,"enforce definiteness in system");

    parse_args.Add("-use_tri_col",&solids_parameters.triangle_collision_parameters.perform_self_collision,"use triangle collisions");
    parse_args.Add("-no_self_interior",&self_collide_surface_only,"do not process penalty self collisions against interior particles");
    parse_args.Add("-use_vanilla_newton",&use_vanilla_newton,"use vanilla newton");
    parse_args.Add("-use_inf_norm",&backward_euler_evolution->minimization_system.use_l_inf_norm,"use l infinity norm to test for convergence");
    parse_args.Add("-density",&density,"density","density");
    parse_args.Add("-gauss",&gauss_order,"order","order of gaussian quadrature to use (max 7)");
    parse_args.Add("-threads",&threads,"threads","number of threads");
    parse_args.Parse();

    tests.data_directory=data_directory;
    LOG::cout<<"Running Standard Test Number "<<test_number<<std::endl;
    output_directory=LOG::sprintf("Test_%d",test_number);
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
    Set_Number_Of_Threads(threads);

    unit_rho=kg/pow<TV::m>(m);
    unit_N=kg*m/(s*s);
    unit_p=unit_N/(m*m);
    unit_J=unit_N*m;
    thickness_multiplier*=m;
    penalty_collisions_length*=m;
    penalty_collisions_separation*=m;
    input_youngs_modulus*=unit_p;
    ether_drag/=s;
    penalty_collisions_stiffness*=unit_J;
    density*=unit_rho;
    if(backward_euler_evolution){
        backward_euler_evolution->newtons_method.tolerance*=unit_N*s;
        backward_euler_evolution->newtons_method.krylov_tolerance/=sqrt(unit_N*s);
        backward_euler_evolution->minimization_objective.collision_thickness*=m;
        backward_euler_evolution->test_diff=test_forces;}

    solids_parameters.use_trapezoidal_rule_for_velocities=!use_newmark_be;

    if(use_penalty_collisions || use_constraint_collisions){
        solids_parameters.triangle_collision_parameters.perform_per_collision_step_repulsions=false;
        solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions=false;
        solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions=false;}
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{
}
//#####################################################################
// Function Set_Number_Of_Threads
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Number_Of_Threads(int threads)
{
#ifdef USE_OPENMP
    PHYSBAM_ASSERT(threads>0);
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp master
    {
        PHYSBAM_ASSERT(threads==omp_get_num_threads());
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#endif
}
//#####################################################################
// Function Post_Initialization
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Post_Initialization()
{
}
//#####################################################################
// Function Postprocess_Solids_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Postprocess_Solids_Substep(const T time,const int substep)
{
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Apply_Constraints(const T dt,const T time)
{
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_External_Forces(ARRAY_VIEW<TWIST<TV> > wrench,const T time)
{
}
//#####################################################################
// Function Add_External_Impulses_Before
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_External_Impulses_Before(ARRAY_VIEW<TV> V,const T time,const T dt)
{
}
//#####################################################################
// Function Add_External_Impulses
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_External_Impulses(ARRAY_VIEW<TV> V,const T time,const T dt)
{
}
//#####################################################################
// Function Add_External_Impulse
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_External_Impulse(ARRAY_VIEW<TV> V,const int node,const T time,const T dt)
{
}
//#####################################################################
// Function Limit_Solids_Dt
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Limit_Solids_Dt(T& dt,const T time)
{
    dt=this->fixed_dt;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_External_Velocities(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_External_Positions(ARRAY_VIEW<FRAME<TV> > frame,const T time)
{
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TWIST<TV> > twist,const T velocity_time,const T current_position_time)
{
}
//#####################################################################
// Function Align_Deformable_Bodies_With_Rigid_Bodies
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Align_Deformable_Bodies_With_Rigid_Bodies()
{
}
//#####################################################################
// Function Preprocess_Solids_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Preprocess_Solids_Substep(const T time,const int substep)
{
}
//#####################################################################
// Function Update_Solids_Parameters
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Update_Solids_Parameters(const T time)
{
}
//#####################################################################
// Function Self_Collisions_Begin_Callback
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Self_Collisions_Begin_Callback(const T time,const int substep)
{
}
//#####################################################################
// Function Filter_Velocities
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Filter_Velocities(const T dt,const T time,const bool velocity_update)
{
}
//#####################################################################
// Function Zero_Out_Enslaved_Position_Nodes
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Zero_Out_Enslaved_Position_Nodes(ARRAY_VIEW<TV> X,const T time)
{
}
//#####################################################################
// Function Preprocess_Frame
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Preprocess_Frame(const int frame)
{
}
//#####################################################################
// Function Add_External_Forces
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_External_Forces(ARRAY_VIEW<TV> F,const T time)
{
}
//#####################################################################
// Function Update_Time_Varying_Material_Properties
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Update_Time_Varying_Material_Properties(const T time)
{
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Postprocess_Substep(const T dt,const T time)
{
}
//#####################################################################
// Function Postprocess_Frame
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Postprocess_Frame(const int frame)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    for(int i=0;i<constrained_particles.m;i++) Add_Debug_Particle(particles.X(constrained_particles(i)),VECTOR<T,3>(1,0,0));
    for(int i=0;i<externally_forced.m;i++) Add_Debug_Particle(particles.X(externally_forced(i)),VECTOR<T,3>(0,1,0));
    for(int i=0;i<kinematic_points.m;i++) Add_Debug_Particle(particles.X(kinematic_points(i)),VECTOR<T,3>(0,1,1));
}
//#####################################################################
// Function Get_Initial_Data
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Get_Initial_Data_After(bool automatically_add_to_collision_structures)
{
    // deformable bodies
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    SOFT_BINDINGS<TV>& soft_bindings=solid_body_collection.deformable_body_collection.soft_bindings;

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
template<class TV> void STANDARD_TESTS_BASE<TV>::
Initialize_Bodies_After()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    if(ether_drag) solid_body_collection.Add_Force(new ELASTIC_ETHER_DRAG<TV>(deformable_body_collection.particles,true,ether_drag,1,save_dt));

    if(use_penalty_collisions)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object))
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);
            solid_body_collection.Add_Force(new IMPLICIT_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,
                    iot,penalty_collisions_stiffness,penalty_collisions_separation,penalty_collisions_length));}
    if(use_constraint_collisions && backward_euler_evolution)
        for(int b=0;b<rigid_body_collection.rigid_body_particles.number;b++){
            IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> > *iot=rigid_body_collection.Rigid_Body(b).implicit_object;
            if(LEVELSET_IMPLICIT_OBJECT<TV>* lio=dynamic_cast<LEVELSET_IMPLICIT_OBJECT<TV>*>(iot->object_space_implicit_object)){
                lio->levelset.interpolation=new CUBIC_SPLINE_INTERPOLATION_UNIFORM<TV,T,FACE_LOOKUP_UNIFORM<TV> >();
                iot=new IMPLICIT_OBJECT_TRANSFORMED<TV,FRAME<TV> >(new SMOOTH_LEVELSET_IMPLICIT_OBJECT<TV>(lio->levelset.grid,lio->levelset.phi),true,iot->transform);}
            backward_euler_evolution->minimization_objective.collision_objects.Append(iot);
            backward_euler_evolution->minimization_objective.coefficient_of_friction.Append(rigid_body_collection.Rigid_Body(b).coefficient_of_friction);}
    else if(!use_penalty_collisions)
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            deformable_body_collection.collisions.collision_structures.Append(deformable_body_collection.structures(i));
            if(solids_parameters.triangle_collision_parameters.perform_self_collision)
                solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.Append(deformable_body_collection.structures(i));}

    if(use_penalty_self_collisions){
        for(int b=0;b<deformable_body_collection.structures.m;b++){
            if(T_OBJECT* object=dynamic_cast<T_OBJECT*>(deformable_body_collection.structures(b))){
                DEFORMABLE_PARTICLES<TV>& undeformed_particles=*particles.Clone();
                T_SURFACE& surface=object->Get_Boundary_Object();
                T_SURFACE& undeformed_surface=*new T_SURFACE(surface.mesh,undeformed_particles);
                LEVELSET_IMPLICIT_OBJECT<TV>& undeformed_levelset=*tests.Initialize_Implicit_Surface(undeformed_surface,10);
                DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>* coll=
                    new DEFORMABLE_OBJECT_COLLISION_PENALTY_FORCES<TV>(particles,undeformed_particles,*object,
                        undeformed_surface,undeformed_levelset,penalty_collisions_stiffness,penalty_collisions_separation);
                if(self_collide_surface_only){
                    coll->colliding_particles=surface.mesh.elements.Flattened();
                    Prune_Duplicates(coll->colliding_particles);}
                int force_id=solid_body_collection.Add_Force(coll);
                if(backward_euler_evolution) backward_euler_evolution->minimization_objective.deformables_forces_lazy.Set(force_id);}}}

    if(enforce_definiteness) solid_body_collection.Enforce_Definiteness(true);
    if(backward_euler_evolution) backward_euler_evolution->minimization_objective.Disable_Current_Colliding_Pairs(0);

    solids_evolution->fully_implicit=true;
    for(int i=0;i<deformable_body_collection.deformables_forces.m;i++)
        if(COLLISION_FORCE<TV>* cf=dynamic_cast<COLLISION_FORCE<TV>*>(solid_body_collection.deformable_body_collection.deformables_forces(i)))
            cf->coefficient_of_friction=input_friction;
}
//#####################################################################
// Function Set_External_Velocities
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_External_Velocities(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    V.Subset(constrained_particles)=constrained_velocities;
    for(int i=0;i<kinematic_points.m;i++) V(kinematic_points(i))=point_curves(i).Derivative(velocity_time);
}
//#####################################################################
// Function Set_External_Positions
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_External_Positions(ARRAY_VIEW<TV> X,const T time)
{
    for(int i=0;i<kinematic_points.m;i++) X(kinematic_points(i))=point_curves(i).Value(time);
}
//#####################################################################
// Function Zero_Out_Enslaved_Velocity_Nodes
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Zero_Out_Enslaved_Velocity_Nodes(ARRAY_VIEW<TV> V,const T velocity_time,const T current_position_time)
{
    V.Subset(constrained_particles).Fill(TV());
    V.Subset(externally_forced).Fill(TV());
    V.Subset(kinematic_points).Fill(TV());
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Read_Output_Files_Solids(const int frame)
{
    BASE::Read_Output_Files_Solids(frame);
    solid_body_collection.Update_Simulated_Particles();
}
//#####################################################################
// Function Set_Kinematic_Positions
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Set_Kinematic_Positions(FRAME<TV>& frame,const T time,const int id)
{
    for(int i=0;i<kinematic_points.m;i++)
        if(id==kinematic_points(i)){
            frame=curves(i).Value(time);
            break;}
}
//#####################################################################
// Function Set_Kinematic_Velocities
//#####################################################################
template<class TV> bool STANDARD_TESTS_BASE<TV>::
Set_Kinematic_Velocities(TWIST<TV>& twist,const T time,const int id)
{
    for(int i=0;i<kinematic_points.m;i++)
        if(id==kinematic_points(i)){
            twist=curves(i).Derivative(time);
            return true;}
    return false;
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Preprocess_Substep(const T dt,const T time)
{
    save_dt=dt;
    if(test_forces){
        solid_body_collection.deformable_body_collection.Test_Energy(time);
        solid_body_collection.deformable_body_collection.Test_Force_Derivatives(time);}
}
//#####################################################################
// Function Add_Constitutive_Model
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Constitutive_Model(T_OBJECT& object,T stiffness,T poissons_ratio,T damping)
{
    if(input_poissons_ratio!=-1) poissons_ratio=input_poissons_ratio;
    if(input_youngs_modulus!=0) stiffness=input_youngs_modulus;
    solid_body_collection.Add_Force(Create_Finite_Volume(object,new COROTATED_FIXED<T,TV::m>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier)));

    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    if(damping && damping_multiplier){
        DEFORMABLES_FORCES<TV>* force=Create_Finite_Volume(object,new COROTATED_FIXED<T,TV::m>(stiffness*stiffness_multiplier,poissons_ratio,damping*damping_multiplier));
        force->Update_Position_Based_State(0,true,true);
        solid_body_collection.Add_Force(new RALEIGH_DAMPING_FORCE<TV>(particles,force,damping*damping_multiplier,1,save_dt));}
}
//#####################################################################
// Function Add_Gravity
//#####################################################################
template<class TV> GRAVITY<TV>& STANDARD_TESTS_BASE<TV>::
Add_Gravity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    GRAVITY<TV>* g=new GRAVITY<TV>(deformable_body_collection.particles,solid_body_collection.rigid_body_collection,true,true);
    solid_body_collection.Add_Force(g);
    return *g;
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}

