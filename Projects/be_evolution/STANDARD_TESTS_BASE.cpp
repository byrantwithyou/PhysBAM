//#####################################################################
// Copyright 2011.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Parsing/PARSE_ARGS.h>
#include <Geometry/Basic_Geometry/TRIANGLE_3D.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Geometry/Topology_Based_Geometry/TETRAHEDRALIZED_VOLUME.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_AREA.h>
#include <Geometry/Topology_Based_Geometry/TRIANGULATED_SURFACE.h>
#include <Rigids/Forces_And_Torques/RIGID_PENALTY_WITH_FRICTION.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Bindings/SOFT_BINDINGS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS_AND_COLLISIONS_GEOMETRY.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Forces/IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Forces/SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Collisions/PENALTY_FORCE_COLLECTION.h>
#include <Solids/Forces_And_Torques/RIGID_DEFORMABLE_PENALTY_WITH_FRICTION.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_OBJECTIVE.h>
#include "STANDARD_TESTS_BASE.h"
namespace PhysBAM{
template<class TV> STANDARD_TESTS_BASE<TV>::
STANDARD_TESTS_BASE(const STREAM_TYPE stream_type_input,PARSE_ARGS& parse_args)
    :BASE(stream_type_input,parse_args),tests(stream_type_input,data_directory,solid_body_collection),
    backward_euler_evolution(new BACKWARD_EULER_EVOLUTION<TV>(solids_parameters,solid_body_collection,*this))
{
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
    parse_args.Add("-test_forces",&test_forces,"test force derivatives");
    parse_args.Add("-print_energy",&solid_body_collection.print_energy,"print energy statistics");
    parse_args.Add("-cgsolids",&solids_parameters.implicit_solve_parameters.cg_tolerance,"tolerance","CG tolerance for backward Euler");
    parse_args.Add("-projection_iterations",&solids_parameters.implicit_solve_parameters.cg_projection_iterations,"iterations","number of iterations used for projection in cg");
    parse_args.Add("-test_system",&solids_parameters.implicit_solve_parameters.test_system,"test system");
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
    parse_args.Add_Not("-no_collisions_in_solve",&backward_euler_evolution->minimization_objective.collisions_in_solve,"disable collisions in solve");
    parse_args.Add("-use_tri_col",&solids_parameters.triangle_collision_parameters.perform_self_collision,"use triangle collisions");
    parse_args.Add("-use_vanilla_newton",&use_vanilla_newton,"use triangle collisions");
    parse_args.Add("-threads",&threads,"threads","Number of threads");
    parse_args.Add("-rd_stiffness",&rd_penalty_stiffness,"stiffness","rigid-deformable penalty force stiffness");
    parse_args.Add("-rd_friction",&rd_penalty_friction,"friction","rigid-deformable penalty force friction");
    parse_args.Add("-grad_ls",&backward_euler_evolution->newtons_method.use_gradient_magnitude_objective,"do line searches on norm of gradient");
    parse_args.Add("-rd",&use_rd,"enable rigid-deformable penalty force friction");
    parse_args.Add("-dd",&use_dd,"enable deformable-deformable penalty force friction");
    parse_args.Add("-rr",&use_rr,"enable rigid-rigid penalty force friction");
    parse_args.Add("-di",&use_di,"enable deformable-object penalty force friction");
    parse_args.Add("-rd_k",&rd_k,&use_rd_k,"stiffness","override stiffness for rigid-deformable penalty force friction");
    parse_args.Add("-dd_k",&dd_k,&use_dd_k,"stiffness","override stiffness for deformable-deformable penalty force friction");
    parse_args.Add("-rr_k",&rr_k,&use_rr_k,"stiffness","override stiffness for rigid-rigid penalty force friction");
    parse_args.Add("-di_k",&di_k,&use_di_k,"stiffness","override stiffness for deformable-object penalty force friction");
    parse_args.Add("-rd_mu",&rd_mu,&use_rd_mu,"friction","override friction for rigid-deformable penalty force friction");
    parse_args.Add("-dd_mu",&dd_mu,&use_dd_mu,"friction","override friction for deformable-deformable penalty force friction");
    parse_args.Add("-rr_mu",&rr_mu,&use_rr_mu,"friction","override friction for rigid-rigid penalty force friction");
    parse_args.Add("-di_mu",&di_mu,&use_di_mu,"friction","override friction for deformable-object penalty force friction");
    parse_args.Add("-bisection",&use_bisection,"use bisection relaxation");
    parse_args.Parse(true);

#ifdef USE_OPENMP
    omp_set_num_threads(threads);
#pragma omp parallel
#pragma omp single
    {
        if(omp_get_num_threads()!=threads) PHYSBAM_FATAL_ERROR();
        LOG::cout<<"Running on "<<threads<<" threads"<<std::endl;
    }
#else
    PHYSBAM_ASSERT(threads==1);
#endif

    tests.data_directory=data_directory;
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
    if(backward_euler_evolution){
        backward_euler_evolution->newtons_method.tolerance*=unit_N*s;
        backward_euler_evolution->newtons_method.krylov_tolerance/=sqrt(unit_N*s);
        backward_euler_evolution->minimization_objective.collision_thickness*=m;
        backward_euler_evolution->test_diff=test_forces;}
    solids_parameters.use_trapezoidal_rule_for_velocities=!use_newmark_be;
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> STANDARD_TESTS_BASE<TV>::
~STANDARD_TESTS_BASE()
{
}
//#####################################################################
// Function After_Get_Initial_Data
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
After_Get_Initial_Data(bool automatically_add_to_collision_structures)
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
// Function After_Initialize_Bodies
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
After_Initialize_Bodies()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    if(use_rd || use_di || use_dd || use_rr){
        Init_Penalty_Collection();
        backward_euler_evolution->asymmetric_system=true;}

    if((use_rd || use_rr) && rigid_body_collection.rigid_body_particles.number>0){
        move_rb_diff.Resize(rigid_body_collection.rigid_body_particles.number);
        backward_euler_evolution->minimization_objective.move_rb_diff=&move_rb_diff;}

    if(backward_euler_evolution) backward_euler_evolution->minimization_objective.Disable_Current_Colliding_Pairs(0);

    solids_evolution->fully_implicit=true;
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Preprocess_Substep(const T dt,const T time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    if(test_forces){
        deformable_body_collection.Test_Energy(time);
        deformable_body_collection.Test_Force_Derivatives(time);}
    if(pfd) pfd->Save_State();
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Postprocess_Substep(const T dt,const T time)
{
}
//#####################################################################
// Function Add_Collision_Object
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Add_Collision_Object(IMPLICIT_OBJECT<TV>* io)
{
    if(!use_di) return;
    Init_Penalty_Collection();
    pfd->di_penalty->ios.Append(io);
}
//#####################################################################
// Function Init_Penalty_Collection
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Init_Penalty_Collection()
{
    if(pfd) return;
    if(!use_rd && !use_di && !use_dd && !use_rr) return;
    detection_grid.Initialize(TV_INT()+100,RANGE<TV>::Centered_Box()*10,true);
    pfd=new PENALTY_FORCE_COLLECTION<TV>(solid_body_collection,
        solid_body_collection.deformable_body_collection.simulated_particles,this->move_rb_diff);
    pfd->Init(&solids_parameters.triangle_collision_parameters,
        use_di,use_dd,use_rd,use_rr);
    if(use_di) pfd->di_penalty->stiffness_coefficient=use_di_k?di_k:rd_penalty_stiffness;
    if(use_dd) pfd->dd_penalty->stiffness_coefficient=use_dd_k?dd_k:rd_penalty_stiffness;
    if(use_rd) pfd->rd_penalty->stiffness_coefficient=use_rd_k?rd_k:rd_penalty_stiffness;
    if(use_rr) pfd->rr_penalty->stiffness_coefficient=use_rr_k?rr_k:rd_penalty_stiffness;
    if(use_di) pfd->di_penalty->friction=use_di_mu?di_mu:rd_penalty_friction;
    if(use_dd) pfd->dd_penalty->friction=use_dd_mu?dd_mu:rd_penalty_friction;
    if(use_rd) pfd->rd_penalty->friction=use_rd_mu?rd_mu:rd_penalty_friction;
    if(use_rr) pfd->rr_penalty->friction=use_rr_mu?rr_mu:rd_penalty_friction;
    if(use_di) pfd->di_penalty->use_bisection=use_bisection;
    if(use_rd) pfd->rd_penalty->use_bisection=use_bisection;
    if(use_rr) pfd->rr_penalty->use_bisection=use_bisection;

    if(backward_euler_evolution) backward_euler_evolution->minimization_objective.pfd=pfd;
}
//#####################################################################
// Function Read_Output_Files_Solids
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Read_Output_Files_Solids(const int frame)
{
    BASE::Read_Output_Files_Solids(frame);
    if(pfd)
    {
        Read_From_File(LOG::sprintf("%s/%d/pfd_data",output_directory.c_str(),frame),pfd->grid);
        pfd->restarted=true;
    }
    if(pfd->di_penalty)
        Read_From_File(LOG::sprintf("%s/%d/di_data",output_directory.c_str(),frame),*pfd->di_penalty);
    if(pfd->rr_penalty)
        Read_From_File(LOG::sprintf("%s/%d/rr_data",output_directory.c_str(),frame),*pfd->rr_penalty);
    if(pfd->rd_penalty)
        Read_From_File(LOG::sprintf("%s/%d/rd_data",output_directory.c_str(),frame),*pfd->rd_penalty);
    if(pfd->dd_penalty)
        Read_From_File(LOG::sprintf("%s/%d/dd_data",output_directory.c_str(),frame),*pfd->dd_penalty);
}
//#####################################################################
// Function Write_Output_Files
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Write_Output_Files(const int frame) const
{
    BASE::Write_Output_Files(frame);
    if(pfd)
    {
#pragma omp task
        if(pfd->di_penalty)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/di_data",output_directory.c_str(),frame),*pfd->di_penalty);
#pragma omp task
        if(pfd->rr_penalty)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/rr_data",output_directory.c_str(),frame),*pfd->rr_penalty);
#pragma omp task
        if(pfd->rd_penalty)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/rd_data",output_directory.c_str(),frame),*pfd->rd_penalty);
#pragma omp task
        if(pfd->dd_penalty)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/dd_data",output_directory.c_str(),frame),*pfd->dd_penalty);
        if(pfd)
            Write_To_File(stream_type,LOG::sprintf("%s/%d/pfd_data",output_directory.c_str(),frame),pfd->grid);
    }
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
