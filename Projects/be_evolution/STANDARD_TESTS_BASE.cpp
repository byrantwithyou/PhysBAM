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
    parse_args.Add("-rd_stiffness",&rd_penalty_stiffness,&use_rd_penalty,"stiffness","rigid-deformable penalty force stiffness");
    parse_args.Add("-rd_friction",&rd_penalty_friction,"friction","rigid-deformable penalty force friction");
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
    const_repulsion_thickness*=m;
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
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m-1>::OBJECT T_SURFACE;
    typedef typename TOPOLOGY_BASED_SIMPLEX_POLICY<TV,TV::m>::OBJECT T_OBJECT;

    if(use_rd_penalty && rigid_body_collection.rigid_body_particles.number>0){
        backward_euler_evolution->asymmetric_system=true;
        move_rb_diff.Resize(rigid_body_collection.rigid_body_particles.number);
        backward_euler_evolution->minimization_objective.move_rb_diff=&move_rb_diff;

        rd_penalty=new RIGID_DEFORMABLE_PENALTY_WITH_FRICTION<TV>(
            particles,rigid_body_collection,move_rb_diff,
            rd_penalty_stiffness,rd_penalty_friction);
        rd_penalty->get_candidates=[this](){Get_RD_Collision_Candidates();};
        solid_body_collection.Add_Force(rd_penalty);

        rr_penalty=new RIGID_PENALTY_WITH_FRICTION<TV>(
            rigid_body_collection,move_rb_diff,
            rd_penalty_stiffness,rd_penalty_friction);
        rr_penalty->get_candidates=[this](){Get_RR_Collision_Candidates();};
        solid_body_collection.Add_Force(rr_penalty);

        dd_penalty=new SELF_COLLISION_PENALTY_FORCE_WITH_FRICTION<TV>(
            particles,rd_penalty_stiffness,rd_penalty_friction);
        dd_penalty->get_candidates=[this](){Get_DD_Collision_Candidates();};
        solid_body_collection.Add_Force(dd_penalty);
        deformable_body_collection.triangle_collisions.compute_edge_edge_collisions=false;
        for(int i=0;i<deformable_body_collection.structures.m;i++){
            auto* s=deformable_body_collection.structures(i);
            if(auto* p=dynamic_cast<T_SURFACE*>(s))
                dd_penalty->Add_Surface(*p);
            else if(auto* p=dynamic_cast<T_OBJECT*>(s))
                dd_penalty->Add_Surface(p->Get_Boundary_Object());}
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.Build_Collision_Geometry();
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free=particles.X;
        repulsion_thickness.Resize(particles.number,true,true,const_repulsion_thickness);
        recently_modified.Resize(particles.number,true,true,true);

        deformable_body_collection.triangle_repulsions_and_collisions_geometry.Initialize(solids_parameters.triangle_collision_parameters);
        deformable_body_collection.triangle_repulsions.Initialize(solids_parameters.triangle_collision_parameters);
        deformable_body_collection.triangle_collisions.Initialize(solids_parameters.triangle_collision_parameters);

        di_penalty=new IMPLICIT_OBJECT_PENALTY_FORCE_WITH_FRICTION<TV>(
            particles,rd_penalty_stiffness,rd_penalty_friction);
        di_penalty->get_candidates=[this](){Get_DI_Collision_Candidates();};
        solid_body_collection.Add_Force(di_penalty);}

    if(backward_euler_evolution) backward_euler_evolution->minimization_objective.Disable_Current_Colliding_Pairs(0);

    solids_evolution->fully_implicit=true;
}
//#####################################################################
// Function Get_RD_Collision_Candidates
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Get_RD_Collision_Candidates()
{
    // TODO use BOX_HIERARCHY
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    for(int d=0;d<particles.number;d++){
        TV X=particles.X(d);
        for(int r=0;r<rigid_body_collection.rigid_body_particles.number;r++){
            if(rigid_body_collection.Rigid_Body(r).Implicit_Geometry_Extended_Value(X)<0)
                rd_penalty->Add_Pair(d,r);}}
}
//#####################################################################
// Function Get_DD_Collision_Candidates
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Get_DD_Collision_Candidates()
{
    typedef typename BASIC_SIMPLEX_POLICY<TV,TV::m>::SIMPLEX_FACE T_FACE;
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    ARRAY<TV>& Xn=deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free;
    deformable_body_collection.triangle_collisions.Update_Swept_Hierachies_And_Compute_Pairs(
        particles.X,Xn,recently_modified,const_repulsion_thickness);

    LOG::printf("num candidate pairs: %i\n",deformable_body_collection.triangle_collisions.point_face_pairs_internal.m);
    for(int i=0;i<deformable_body_collection.triangle_collisions.point_face_pairs_internal.m;i++){ // p f
        VECTOR<int,TV::m+1> pf=deformable_body_collection.triangle_collisions.point_face_pairs_internal(i);
        int p=pf(0);
        TV_INT f=pf.Remove_Index(0);

        // Detection is expensive; make sure we are not already known.
        auto s_e=dd_penalty->object_from_element.Get(f.Sorted());
        if(dd_penalty->hash.Contains({p,s_e.x})) continue;
        const auto& ts=*dd_penalty->surfaces(s_e.x);
        PHYSBAM_ASSERT(f==ts.mesh.elements(s_e.y));

        // Particle must have exited.
        if(ts.Get_Element(s_e.y).Signed_Distance(particles.X(p))>0) continue;

        // Do the expensive check.
        T_FACE face(Xn.Subset(f));
        T collision_time=0;
        TV normal;
        VECTOR<T,TV::m+1> weights;
        VECTOR<TV,TV::m> V_f(particles.X.Subset(f)-Xn.Subset(f));
        bool in=face.Point_Face_Collision(Xn(p),particles.X(p)-Xn(p),V_f,1,
            const_repulsion_thickness,collision_time,normal,weights,false);
        if(!in) continue;
        
        dd_penalty->Add_Pair(p,s_e.x,weights.Remove_Index(0),s_e.y);}
}
//#####################################################################
// Function Get_DI_Collision_Candidates
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Get_DI_Collision_Candidates()
{
    // TODO use BOX_HIERARCHY
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    for(int b=0;b<di_penalty->ios.m;b++){
        const IMPLICIT_OBJECT<TV>& io=*di_penalty->ios(b);
        for(int p=0;p<particles.number;p++){
            if(io.Extended_Phi(particles.X(p))<0)
                di_penalty->Add_Pair(p,b);}}
}
//#####################################################################
// Function Get_RR_Collision_Candidates
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Get_RR_Collision_Candidates()
{
}
//#####################################################################
// Function Preprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Preprocess_Substep(const T dt,const T time)
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;
    if(test_forces){
        deformable_body_collection.Test_Energy(time);
        deformable_body_collection.Test_Force_Derivatives(time);}
    if(dd_penalty){
        deformable_body_collection.triangle_repulsions_and_collisions_geometry.X_self_collision_free=particles.X;
        repulsion_thickness.Resize(particles.number,true,true,const_repulsion_thickness);
        recently_modified.Resize(particles.number);
        recently_modified.Fill(true);}
}
//#####################################################################
// Function Postprocess_Substep
//#####################################################################
template<class TV> void STANDARD_TESTS_BASE<TV>::
Postprocess_Substep(const T dt,const T time)
{
    if(rd_penalty)
        rd_penalty->Update_Attachments_And_Prune_Pairs();
    if(di_penalty)
        di_penalty->Update_Attachments_And_Prune_Pairs();
    if(dd_penalty)
        dd_penalty->Update_Attachments_And_Prune_Pairs();
    if(rr_penalty)
        rr_penalty->Update_Attachments_And_Prune_Pairs();
}
template class STANDARD_TESTS_BASE<VECTOR<float,2> >;
template class STANDARD_TESTS_BASE<VECTOR<float,3> >;
template class STANDARD_TESTS_BASE<VECTOR<double,2> >;
template class STANDARD_TESTS_BASE<VECTOR<double,3> >;
}
