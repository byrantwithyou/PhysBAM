//#####################################################################
// Copyright 2006-2009, Zhaosheng Bao, Ron Fedkiw, Geoffrey Irving, Sergey Koltakov, Nipun Kwatra, Michael Lentine, Craig Schroeder, Andrew Selle, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class NEWMARK_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_RESIDUAL.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/LANCZOS_ITERATION.h>
#include <PhysBAM_Tools/Krylov_Solvers/SYMMQMR.h>
#include <PhysBAM_Tools/Log/DEBUG_PRINT.h>
#include <PhysBAM_Tools/Log/DEBUG_SUBSTEPS.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Geometry/Implicit_Objects/IMPLICIT_OBJECT_TRANSFORMED.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/SOFT_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/DEFORMABLE_OBJECT_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Collisions_And_Interactions/TRIANGLE_REPULSIONS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/BINDING_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Forces/INCOMPRESSIBLE_FINITE_VOLUME.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_1D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_2D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Articulated_Rigid_Bodies/ARTICULATED_RIGID_BODY_3D.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Collisions/RIGID_BODY_INTERSECTIONS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Forces_And_Torques/RIGID_LINEAR_SPRINGS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_COLLISION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Bodies/RIGID_BODY_EVOLUTION_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Rigids/Rigid_Body_Clusters/RIGID_BODY_CLUSTER_BINDINGS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGID_DEFORMABLE_COLLISIONS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Collisions/RIGIDS_NEWMARK_COLLISION_CALLBACKS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/NEWMARK_EVOLUTION.h>
#include <stdexcept>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> NEWMARK_EVOLUTION<TV>::
NEWMARK_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input)
    :SOLIDS_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input),rigids_evolution_callbacks(*new RIGIDS_NEWMARK_COLLISION_CALLBACKS<TV>(*this)),repulsions(0),
    use_existing_contact(false),print_matrix(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> NEWMARK_EVOLUTION<TV>::
~NEWMARK_EVOLUTION()
{
    delete &rigids_evolution_callbacks;
}
//#####################################################################
// Function Prepare_Backward_Euler_System
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Prepare_Backward_Euler_System(BACKWARD_EULER_SYSTEM<TV>& system,const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    rigid_body_collection.Update_Angular_Velocity(); // make sure omega = I^{-1} L

    B_full.Resize(particles.Size(),false,false);rigid_B_full.Resize(rigid_body_particles.Size(),false,false);

    GENERALIZED_VELOCITY<TV> B_all(B_full,rigid_B_full,solid_body_collection);
    GENERALIZED_VELOCITY<TV> V_all(particles.V,rigid_body_particles.twist,solid_body_collection);
    KRYLOV_SOLVER<T>::Ensure_Size(krylov_vectors,V_all,1); // Ensure Finish_Backward_Euler_Step can run successfully

    B_full.Subset(solid_body_collection.deformable_body_collection.simulated_particles).Fill(TV());rigid_B_full.Fill(TWIST<TV>());
    solid_body_collection.example_forces_and_velocities->Add_External_Forces(B_full,current_velocity_time+dt);
    solid_body_collection.example_forces_and_velocities->Add_External_Forces(rigid_B_full,current_velocity_time+dt);
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(particles.V);
    solid_body_collection.Add_Velocity_Independent_Forces(B_full,rigid_B_full,current_velocity_time+dt); // this is a nop for binding forces
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Distribute_Force_To_Parents(rigid_B_full);

    if(solid_body_collection.deformable_body_collection.soft_bindings.Need_Bindings_Mapped()){
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        solid_body_collection.deformable_body_collection.soft_bindings.Map_Forces_From_Parents(B_full,rigid_B_full);
        solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(B_full);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        for(int k=0;k<solid_body_collection.solids_forces.m;k++) if(dynamic_cast<BINDING_SPRINGS<TV>*>(&*solid_body_collection.solids_forces(k)))
            solid_body_collection.solids_forces(k)->Add_Force_Differential(particles.X,B_full,current_velocity_time+dt);
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data_Global(B_full);
        solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full,rigid_B_full);}

    Initialize_World_Space_Masses();
    if(articulated_rigid_body.constrain_pd_directions){
        for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass_inverse(p)*rigid_B_full(p);}
        articulated_rigid_body.Poststabilization_Projection(rigid_B_full,true);
        for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);rigid_B_full(p)=world_space_rigid_mass(p)*rigid_B_full(p);}}

    for(int i=0;i<solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
        B_full(p)=particles.V(p)+dt*particles.one_over_mass(p)*B_full(p);}
    for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_B_full(p)=rigid_body_particles.twist(p)+world_space_rigid_mass_inverse(p)*rigid_B_full(p)*dt;}

    V_all=B_all;
    rigid_body_collection.Update_Angular_Momentum();
    Diagnostics(dt,current_position_time,0,0,604,"Before boundary conditions");
    system.Set_Global_Boundary_Conditions(V_all,X_save,rigid_frame_save,rigid_velocity_save,rigid_angular_momentum_save,V_save,
        solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
    if(velocity_update && !solids_parameters.use_post_cg_constraints){
        Apply_Constraints(dt,current_velocity_time);}

    // TODO: This is completely broken for trapezoid and likely for BE as well.
    if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators() && solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
        saved_pd=rigid_body_particles.twist;
        solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(saved_pd,true);
        saved_pd=rigid_body_particles.twist-saved_pd;}
    if(!velocity_update) solid_body_collection.example_forces_and_velocities->Add_External_Impulses_Before(B_full,current_position_time,(T)2*dt); // 2*dt is position dt TODO: what time?
}
//#####################################################################
// Function Finish_Backward_Euler_Step
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Finish_Backward_Euler_Step(KRYLOV_SYSTEM_BASE<T>& system,const T dt,const T current_position_time,const bool velocity_update)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
    GENERALIZED_VELOCITY<TV>& F=debug_cast<GENERALIZED_VELOCITY<TV>&>(*krylov_vectors(0));
    rigid_deformable_collisions->rigid_body_collisions.Remove_Contact_Joints(); // TODO: Fix me.

    if(velocity_update && solids_parameters.use_post_cg_constraints){ // return rhs + dt Fd V^n+1 for friction processing
        if(solids_parameters.no_contact_friction){
            /*ARRAY<TV> V_no_projection;V_no_projection.Resize(particles.Size(),false,false);
            ARRAY<TV> V_projected;V_projected.Resize(particles.Size(),false,false);
            ARRAY<TWIST<TV> > rigid_V_no_projection;rigid_V_no_projection.Resize(rigid_body_collection.rigid_body_particles.Size(),false,false);
            ARRAY<TWIST<TV> > rigid_V_projected;rigid_V_projected.Resize(rigid_body_collection.rigid_body_particles.Size(),false,false);
            for(int i=0;i<solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){
                int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                V_no_projection(p)=B_full(p)+dt*particles.one_over_mass(p)*F_full(p);}
            for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){
                int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
                rigid_V_no_projection(p)=rigid_B_full(p)+world_space_rigid_mass_inverse(p)*rigid_F_full(p)*dt;}
            V_projected=V_no_projection;
            GENERALIZED_VELOCITY<TV> V_projected_all(V_projected,rigid_V_projected,solid_body_collection);
            system.Project(V_projected_all);*/}
        else{
            for(int i=0;i<solid_body_collection.deformable_body_collection.dynamic_particles.m;i++){int p=solid_body_collection.deformable_body_collection.dynamic_particles(i);
                particles.V(p)=B_full(p)+dt*particles.one_over_mass(p)*F.V.array(p);}
            for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
                rigid_body_particles.twist(p)=rigid_B_full(p)+world_space_rigid_mass_inverse(p)*F.rigid_V.array(p)*dt;}

            // No friction for these, so reproject them.
            solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(rigid_body_particles.twist,true);
            solid_body_collection.example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(particles.V,current_position_time,current_position_time);
            Zero_Out_Enslaved_Velocity_Nodes(rigid_body_particles.twist,current_position_time,current_position_time);
            
            if(solid_body_collection.rigid_body_collection.articulated_rigid_body.Has_Actuators() && solid_body_collection.rigid_body_collection.articulated_rigid_body.constrain_pd_directions){
                solid_body_collection.rigid_body_collection.articulated_rigid_body.Poststabilization_Projection(rigid_body_particles.twist,true);
                rigid_body_particles.twist+=saved_pd;}
            Diagnostics(dt,current_position_time,0,0,610,"After undo projections");}}
    rigid_body_collection.Update_Angular_Momentum();

    if(!velocity_update) solid_body_collection.example_forces_and_velocities->Add_External_Impulses(particles.V,current_position_time,(T)2*dt); // 2*dt is for

    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(); // TODO: MPI safe?

    rigid_body_collection.Update_Angular_Momentum(); // make sure L = I omega
}
//#####################################################################
// Function Backward_Euler_Step_Velocity_Helper
//#####################################################################
// assumes all solids_forces are linear in velocity, with a symmetric positive definite Jacobian.
template<class TV> void NEWMARK_EVOLUTION<TV>::
Backward_Euler_Step_Velocity_Helper(const T dt,const T current_velocity_time,const T current_position_time,const bool velocity_update)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=rigid_body_collection.rigid_body_particles;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    BACKWARD_EULER_SYSTEM<TV> system(*this,solid_body_collection,dt,current_velocity_time,current_position_time,&solid_body_collection.rigid_body_collection.articulated_rigid_body,
        (velocity_update && solids_parameters.enforce_repulsions_in_cg)?repulsions:0,mpi_solids,velocity_update);

    Prepare_Backward_Euler_System(system,dt,current_velocity_time,current_position_time,velocity_update);

    static CONJUGATE_GRADIENT<T> cg;
    static CONJUGATE_RESIDUAL<T> cr;
    static SYMMQMR<T> symmqmr;
    KRYLOV_SOLVER<T>* solver=0;
    const char* solver_name=0;
    if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cg){solver=&cg;solver_name="CG";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_cr){solver=&cr;solver_name="CONJUGATE_RESIDUAL";}
    else if(solids_parameters.implicit_solve_parameters.evolution_solver_type==krylov_solver_symmqmr){solver=&symmqmr;solver_name="SYMMQMR";}
    solver->print_diagnostics=solid_body_collection.print_diagnostics;
    solver->print_residuals=solid_body_collection.print_residuals;
    solver->iterations_used=&solid_body_collection.iterations_used_diagnostic;
    solver->restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
    system.project_nullspace_frequency=solids_parameters.implicit_solve_parameters.project_nullspace_frequency;

    GENERALIZED_VELOCITY<TV> V(particles.V,rigid_body_particles.twist,solid_body_collection),B(B_full,rigid_B_full,solid_body_collection);
    GENERALIZED_MASS<TV> mass(solid_body_collection); // TODO: Doing duplicate computation of mass.

    if(solids_parameters.implicit_solve_parameters.spectral_analysis){
        V=B;
        KRYLOV_SOLVER<T>::Ensure_Size(krylov_vectors,V,2);
        LANCZOS_ITERATION<T>::Print_Spectral_Information(system,V,*krylov_vectors(0),*krylov_vectors(1),
            (T)1e-2*solids_parameters.implicit_solve_parameters.cg_tolerance,
            solids_parameters.implicit_solve_parameters.lanczos_iterations);}

    LOG::Time(solver_name);
    Diagnostics(dt,current_position_time,0,0,606,"Before solve");
    static int solve_id=-1;solve_id++;
    if(print_matrix){
        LOG::cout<<"solve id "<<solve_id<<std::endl;
        KRYLOV_SOLVER<T>::Ensure_Size(krylov_vectors,V,2);
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("M-%i.txt",solve_id).c_str()).Write("M",system,*krylov_vectors(0),*krylov_vectors(1));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("P-%i.txt",solve_id).c_str()).Write_Projection("P",system,*krylov_vectors(0));
        OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("b-%i.txt",solve_id).c_str()).Write("b",B);}

    if(solids_parameters.implicit_solve_parameters.test_system){
        KRYLOV_SOLVER<T>::Ensure_Size(krylov_vectors,V,3);
        system.Test_System(*krylov_vectors(0),*krylov_vectors(1),*krylov_vectors(2));}
    if(!solver->Solve(system,V,B,krylov_vectors,solids_parameters.implicit_solve_parameters.cg_tolerance,1,solids_parameters.implicit_solve_parameters.cg_iterations) && solids_parameters.implicit_solve_parameters.throw_exception_on_backward_euler_failure)
        throw std::runtime_error("Backward Euler Failed");
    if(print_matrix) OCTAVE_OUTPUT<T>(STRING_UTILITIES::string_sprintf("x-%i.txt",solve_id).c_str()).Write("x",V);
    Diagnostics(dt,current_position_time,0,0,607,"After solve");
    LOG::Stop_Time();

    if(velocity_update && solids_parameters.use_post_cg_constraints)
        system.Force(V,debug_cast<GENERALIZED_VELOCITY<TV>&>(*krylov_vectors(0)));

    Finish_Backward_Euler_Step(system,dt,current_position_time,velocity_update);
}
//#####################################################################
// Function Average_And_Exchange_Position
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Average_And_Exchange_Position()
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    assert(X_save.m==particles.Size());
    assert(rigid_frame_save.m==rigid_body_collection.rigid_body_particles.Size());
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    const ARRAY<int>& simulated_rigid_body_particles=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles;
    INDIRECT_ARRAY<ARRAY<TV> > X_save_fragment(X_save,simulated_particles);
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > X_fragment(particles.X,simulated_particles);
    for(int i=0;i<X_fragment.Size();i++){TV X_average=(T).5*(X_fragment(i)+X_save_fragment(i));X_save_fragment(i)=X_fragment(i);X_fragment(i)=X_average;}
    ARRAY<int> rigid_body_indices(simulated_rigid_body_particles);rigid_body_indices.Append_Elements(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies);
    for(int i=0;i<rigid_body_indices.Size();i++){int p=rigid_body_indices(i);
        FRAME<TV> tmp=FRAME<TV>::Interpolation(rigid_body_collection.rigid_body_particles.frame(p),rigid_frame_save(p),(T).5);
        rigid_frame_save(p)=rigid_body_collection.rigid_body_particles.frame(p);
        rigid_body_collection.rigid_body_particles.frame(p)=tmp;}
    for(int i=0;i<rigid_body_indices.m;i++) rigid_body_collection.Rigid_Body(rigid_body_indices(i)).Update_Angular_Velocity();
}
//#####################################################################
// Function Trapezoidal_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt/2
template<class TV> void NEWMARK_EVOLUTION<TV>::
Trapezoidal_Step_Velocity(const T dt,const T time)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    // save V at time
    Save_Velocity();
    // update V implicitly to time+dt/2
    Backward_Euler_Step_Velocity_Helper(dt/2,time,time+dt/2,true);
    // set up rigid_V_save for extrapolation step
    rigid_V_save.Resize(rigid_body_collection.rigid_body_particles.Size(),false,false);
    for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        rigid_V_save(p).linear=rigid_velocity_save(p).linear;
        rigid_V_save(p).angular=rigid_body_collection.Rigid_Body(p).World_Space_Inertia_Tensor_Inverse_Times(rigid_angular_momentum_save(p));}
    // Use V_n instead of V_save rather than copying it around another time.  Also simplifies state dependencies.
    // extrapolate V to time+dt based on V at time and time+dt/2
    GENERALIZED_VELOCITY<TV> V_n(V_save,rigid_V_save,solid_body_collection),V(particles.V,rigid_body_collection.rigid_body_particles.twist,solid_body_collection);
    V*=(T)2;V-=V_n;

    // enforce boundary conditions again
    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt/2);
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particles.twist,time+dt,time+dt/2);
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
    rigid_body_collection.Update_Angular_Momentum(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies);
}
//#####################################################################
// Function Backward_Euler_Step_Velocity
//#####################################################################
// assumes X is fixed at time+dt
template<class TV> void NEWMARK_EVOLUTION<TV>::
Backward_Euler_Step_Velocity(const T dt,const T time)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;

    Backward_Euler_Step_Velocity_Helper(dt,time,time+dt,true);
    // enforce boundary conditions again
    if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);
    Set_External_Velocities(particles.V,time+dt,time+dt);
    rigid_body_collection.Update_Angular_Velocity(rigid_body_collection.simulated_rigid_body_particles);
    kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particles.twist,time+dt,time+dt);
    
    // enforce boundary conditions again
    solid_body_collection.rigid_body_collection.rigid_body_cluster_bindings.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    rigid_body_collection.Update_Angular_Momentum(rigid_body_collection.simulated_rigid_body_particles);
}
//#####################################################################
// Function Make_Incompressible
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Make_Incompressible(const T dt,const bool correct_volume)
{
    for(int f=0;f<solid_body_collection.deformable_body_collection.deformables_forces.m;f++)
        if(INCOMPRESSIBLE_FINITE_VOLUME_BASE<TV>* fvm=dynamic_cast<INCOMPRESSIBLE_FINITE_VOLUME_BASE<TV>*>(&*solid_body_collection.deformable_body_collection.deformables_forces(f))){
            fvm->Set_Neumann_Boundary_Conditions(&solid_body_collection.deformable_body_collection.collisions.particle_states,repulsions);
            fvm->Make_Incompressible(dt,correct_volume);
            Diagnostics(dt,time,1,0,700,"make incompressible");}
}
//#####################################################################
// Function Advance_One_Time_Step_Position
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time, const bool solids)
{
    PHYSBAM_ASSERT(solids_parameters.use_post_cg_constraints || !solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position Start dt=%f time=%f",dt,time),2,2);
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body; // Needn't be a pointer
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.

    if((solids_parameters.triangle_collision_parameters.perform_self_collision && (solid_body_collection.deformable_body_collection.triangle_repulsions_and_collisions_geometry.structures.m || solids_parameters.triangle_collision_parameters.initialize_collisions_without_objects))
       && (!repulsions && solids_parameters.triangle_collision_parameters.perform_per_time_step_repulsions))
        repulsions=&solid_body_collection.deformable_body_collection.triangle_repulsions;

    Diagnostics(dt,time,0,0,1,"begin integration");
    example_forces_and_velocities.Advance_One_Time_Step_Begin_Callback(dt,time);

    solids_evolution_callbacks->Update_Solids_Parameters(time);
    if(solids){
        rigid_body_collisions->Initialize_Data_Structures();

        // save position and velocity for later trapezoidal rule
        Save_Velocity();
        Save_Position(X_save,rigid_frame_save);

        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(false);}

    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
    solid_body_collection.Update_Position_Based_State(time+dt,true);

    // get momentum difference for v^n -> v^{n+1/2} udpate
    if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);

    Backward_Euler_Step_Velocity_Helper(dt/2,time,time,false); // update V implicitly to time+dt/2

    Diagnostics(dt,time,1,0,5,"backward Euler");
        
    if(solids_parameters.use_projections_in_position_update){
        Restore_Velocity();
        Apply_Projections_In_Position_Update(dt,time);
        Diagnostics(dt,time,1,0,6,"apply projections in position update");}
        
    if(solids_parameters.verbose) Print_Maximum_Velocities(time);
    if(!solids) return; // early exit for fluids only in parallel
        
    Make_Incompressible(dt,true); // adjust velocity to fix volume
    solids_evolution_callbacks->Filter_Velocities(dt,time+dt,false); // use time+dt since these velocities are used to step to time+dt
    if(repulsions){
        repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,true,false);
        Diagnostics(dt,time,1,2,8,"Repulsions");}
    Compute_Momentum_Differences();

    // add collision impulses to time n velocities and save
    Euler_Step_Position(dt,time);
    Diagnostics(dt,time,1,2,10,"Euler step position");
    Exchange_Velocity();
    Diagnostics(dt,time,0,2,11,"restore velocity");

    Process_Collisions(dt,time,advance_rigid_bodies);
    Diagnostics(dt,time,0,2,12,"add elastic collisions");

    Restore_Position(X_save,rigid_frame_save);
    Diagnostics(dt,time,0,0,13,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,14,"poststabilization");}
    Save_Velocity();

    // update positions, apply contact, arb prestabilization and push out
    if(solids_parameters.rigid_body_collision_parameters.perform_contact && advance_rigid_bodies) rigid_body_collisions->Compute_Contact_Graph(dt,time,&articulated_rigid_body);
    Update_Velocity_Using_Stored_Differences(dt/2,time);
    Diagnostics(dt,time,1,0,18,"update velocity using stored differences");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,185,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt/2,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,1,0,19,"solve velocities for pd");}

    Update_Positions_And_Apply_Contact_Forces(dt,time,false);
    Diagnostics(dt,time,1,2,20,"contact, prestabilization");
    if(solids_parameters.rigid_body_collision_parameters.use_push_out){
        if(solids_parameters.use_rigid_deformable_contact) rigid_deformable_collisions->Process_Push_Out();
        else{
            solid_body_collection.deformable_body_collection.collisions.Process_Push_Out();
            rigid_body_collisions->Process_Push_Out_Legacy();}
        Diagnostics(dt,time,1,2,21,"push out");}

    // This is necessary since the velocity update below needs to use soft bound positions, but modifying positions within the velocity update is prohibited.
    if(repulsions) solid_body_collection.deformable_body_collection.soft_bindings.Clamp_Particles_To_Embedded_Positions(true);

    Restore_Velocity();
    Diagnostics(dt,time,0,2,22,"restore velocity");

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Position End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Apply_Projections_In_Position_Update
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Apply_Projections_In_Position_Update(const T dt,const T time)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    // iterate over rigid/rigid
    if(rigid_body_collection.dynamic_rigid_body_particles.m && solids_parameters.rigid_body_collision_parameters.enforce_rigid_rigid_contact_in_cg){
        for(typename HASHTABLE<TRIPLE<int,int,TV> >::ITERATOR iterator(rigid_deformable_collisions->rigid_body_collisions.rigid_body_particles_intersections);iterator.Valid();iterator.Next()){
            const TRIPLE<int,int,TV>& intersection=iterator.Key();
            const RIGID_BODY<TV> &particle_body=rigid_body_collection.Rigid_Body(intersection.x),
                &levelset_body=rigid_body_collection.Rigid_Body(intersection.y);
            TV relative_velocity=RIGID_BODY<TV>::Relative_Velocity(particle_body,levelset_body,intersection.z);
            TV normal=levelset_body.implicit_object->Extended_Normal(intersection.z);
            if(TV::Dot_Product(relative_velocity,normal)>0) rigid_deformable_collisions->rigid_body_collisions.rigid_body_particles_intersections.Delete(iterator.Key());}}

    // iterate over rigid/deformable
    for(int i=0;i<rigid_deformable_collisions->precompute_contact_projections.m;i++){
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION& precompute=*rigid_deformable_collisions->precompute_contact_projections(i);
        const RIGID_BODY<TV>& body=precompute.rigid_body;
        for(int j=0;j<precompute.particles.m;j++){const int p=precompute.particles(j);
            TV V_rel=body.Pointwise_Object_Velocity(particles.X(p))-particles.V(p);
            precompute.V_rel_target(j)=TV();
            if(TV::Dot_Product(V_rel,precompute.N(j))<0){
                precompute.particles.Remove_Index_Lazy(j);
                j--;}}}

    Backward_Euler_Step_Velocity_Helper(dt,time,time,true);
}
//#####################################################################
// Function Write_Position_Update_Projection_Data
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Write_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix)
{
    if(rigid_deformable_collisions->rigid_body_collisions.rigid_body_particles_intersections.Size())
        FILE_UTILITIES::Write_To_File(stream_type,prefix+"projection_data_rigid_rigid",rigid_deformable_collisions->rigid_body_collisions.rigid_body_particles_intersections);
    if(!rigid_deformable_collisions->precompute_contact_projections.m) return;
    std::ostream* output=FILE_UTILITIES::Safe_Open_Output(prefix+"projection_data_rigid_deformable");
    TYPED_OSTREAM typed_output(*output,stream_type);
    Write_Binary(typed_output,rigid_deformable_collisions->precompute_contact_projections.m);
    for(int i=0;i<rigid_deformable_collisions->precompute_contact_projections.m;i++){
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION& p=*rigid_deformable_collisions->precompute_contact_projections(i);
        Write_Binary(typed_output,p.rigid_body.particle_index,p.particles,p.V_rel_target,p.N_over_NT_K_N,p.r,p.N,p.rN,p.A,p.A_inverted);}
    delete output;
}
//#####################################################################
// Function Read_Position_Update_Projection_Data
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Read_Position_Update_Projection_Data(const STREAM_TYPE stream_type,const std::string& prefix)
{
    if(FILE_UTILITIES::File_Exists(prefix+"projection_data_rigid_rigid"))
        FILE_UTILITIES::Read_From_File(stream_type,prefix+"projection_data_rigid_rigid",rigid_deformable_collisions->rigid_body_collisions.rigid_body_particles_intersections);
    if(!FILE_UTILITIES::File_Exists(prefix+"projection_data_rigid_deformable")) return;
    std::istream* input=FILE_UTILITIES::Safe_Open_Input(prefix+"projection_data_rigid_deformable");
    TYPED_ISTREAM typed_input(*input,stream_type);
    int precompute_contact_projections_size;
    Read_Binary(typed_input,precompute_contact_projections_size);
    for(int i=0;i<precompute_contact_projections_size;i++){
        int index;
        Read_Binary(typed_input,index);
        typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION* p=
            new typename RIGID_DEFORMABLE_COLLISIONS<TV>::PRECOMPUTE_CONTACT_PROJECTION(solid_body_collection.rigid_body_collection.Rigid_Body(index),false);
        Read_Binary(typed_input,p->particles,p->V_rel_target,p->N_over_NT_K_N,p->r,p->N,p->rN,p->A,p->A_inverted);}
    delete input;
}
//#####################################################################
// Function Process_Collisions
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Process_Collisions(const T dt,const T time,const bool advance_rigid_bodies)
{
    if(solids_parameters.use_rigid_deformable_contact)
        rigid_deformable_collisions->Add_Elastic_Collisions(dt,time,rigid_frame_save,rigid_angular_momentum_difference,rigid_velocity_difference,rigid_velocity_save,
            rigid_angular_momentum_save,X_save,V_save);
    else if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
        if(solids_parameters.rigid_body_collision_parameters.perform_collisions && advance_rigid_bodies) rigid_body_collisions->Add_Elastic_Collisions(dt,time);
        int interactions=solid_body_collection.deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,
            solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
        if(interactions) LOG::Stat("collision body collisions",interactions);}
}
//#####################################################################
// Function Advance_One_Time_Step_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids)
{
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity Start dt=%f time=%f",dt,time),2,2);

    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=true; //solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;  TODO: Fix this.

    if(solids){
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,0,2,24,"poststabilization");}
        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions) solid_body_collection.deformable_body_collection.collisions.Activate_Collisions(true);

        // initialize data needed for rigid/deformable contact projection in CG
        if(solids_parameters.use_rigid_deformable_contact) rigid_deformable_collisions->Initialize_All_Contact_Projections();
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after creating joints.  before trapezoidal velocities dt=%f time=%f",dt,time),2,2);}

    if(solids_parameters.use_trapezoidal_rule_for_velocities){
        Average_And_Exchange_Position(); // move to positions at time+dt/2 for trapezoidal step
        Diagnostics(dt,time,0,1,25,"average and exchange positions");
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt/2);

        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
        solid_body_collection.Update_Position_Based_State(time+dt/2,(solids_parameters.allow_altitude_spring_change_between_updates?true:false));
        Make_Incompressible(dt,false); // make velocity divergence free
        if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        Trapezoidal_Step_Velocity(dt,time);
        solids_evolution_callbacks->Filter_Velocities(dt,time+dt,true);
        Diagnostics(dt,time,2,1,29,"trazepoid rule");

        Restore_Position(X_save,rigid_frame_save); // move to final positions at time time+dt
        Diagnostics(dt,time,2,2,30,"restore position");
        if(advance_rigid_bodies){
            articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
            Diagnostics(dt,time,2,2,31,"poststabilization");}
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);}
    else{
        example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);

        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data_Global(solid_body_collection.deformable_body_collection.particles.X);
        solid_body_collection.Update_Position_Based_State(time+dt,(solids_parameters.allow_altitude_spring_change_between_updates?true:false));
        Make_Incompressible(dt,false); // make velocity divergence free
        if(articulated_rigid_body.Has_Actuators()) example_forces_and_velocities.Set_PD_Targets(dt,time);
        Backward_Euler_Step_Velocity(dt,time); // TODO: Tamar & Craig, do you need post stab?
        Diagnostics(dt,time,2,2,29,"backward Euler");
        solids_evolution_callbacks->Filter_Velocities(dt,time+dt,true);}

    if(solids){
        PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("after removing joints.  after trapezoidal velocities dt=%f time=%f",dt,time),2,2);

        if(!solids_parameters.no_contact_friction){
            if(solids_parameters.use_post_cg_constraints) Apply_Constraints(dt,time);
            else if(repulsions){
                repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
                Diagnostics(dt,time,2,2,40,"self repulsions");}}
        solid_body_collection.deformable_body_collection.collisions.Reset_Object_Collisions();
        if(advance_rigid_bodies && solids_parameters.rigid_body_evolution_parameters.clamp_rigid_body_velocities){
            Clamp_Velocities(); // TODO: Examples should do this during the Advance_One_Time_Step_End_Callback example callback
            Diagnostics(dt,time,2,2,41,"clamp velocities");}
        if(solids_parameters.verbose) Print_Maximum_Velocities(time);}

    example_forces_and_velocities.Advance_One_Time_Step_End_Callback(dt,time);

    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Advance_One_Time_Step_Velocity End dt=%f time=%f",dt,time),2,2);
}
//#####################################################################
// Function Apply_Constraints
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Apply_Constraints(const T dt,const T time)
{
    ARTICULATED_RIGID_BODY<TV>& articulated_rigid_body=solid_body_collection.rigid_body_collection.articulated_rigid_body;
    const bool advance_rigid_bodies=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m!=0;
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,0,0,345,"poststabilization");
        if(articulated_rigid_body.Has_Actuators() && !articulated_rigid_body.constrain_pd_directions){
            articulated_rigid_body.Compute_Position_Based_State(dt,time);
            articulated_rigid_body.Solve_Velocities_for_PD(time,dt,solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);}
        Diagnostics(dt,time,2,2,35,"solve velocities for pd");}
    Save_Position(X_save_for_constraints,rigid_frame_save_for_constraints);
    if(solids_parameters.rigid_body_collision_parameters.use_persistant_contact){rigid_deformable_collisions->Process_Precomputed_Contact_With_Rigid_Bodies();use_existing_contact=true;}
    Diagnostics(dt,time,2,2,137,"consistent contact");
    Update_Positions_And_Apply_Contact_Forces(dt,time,true);use_existing_contact=false;
    Diagnostics(dt,time,4,2,37,"contact, prestabilization");
    if(rigid_body_collisions->prune_stacks_from_contact) rigid_body_collisions->Construct_Stacks();
    if(rigid_body_collisions->prune_contact_using_velocity) rigid_body_collisions->Compute_Contact_Frequency();
    Restore_Position(X_save_for_constraints,rigid_frame_save_for_constraints);
    Diagnostics(dt,time,2,2,38,"restore position");
    if(advance_rigid_bodies){
        articulated_rigid_body.Apply_Poststabilization(solids_parameters.implicit_solve_parameters.test_system,solids_parameters.implicit_solve_parameters.print_matrix);
        Diagnostics(dt,time,2,2,39,"poststabilization");}
    // modify velocity with inelastic and friction repulsions
    if(repulsions) repulsions->Adjust_Velocity_For_Self_Repulsion_Using_History(dt,false,true);
}
//#####################################################################
// Function Print_Maximum_Velocities
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Print_Maximum_Velocities(const T time) const
{
    LOG::cout<<"time = "<<time<<std::endl;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    int max_index=-1;T max_magnitude_squared=-FLT_MAX;const INDIRECT_ARRAY<ARRAY_VIEW<TV>,ARRAY<int>&>& V=particles.V.Subset(solid_body_collection.deformable_body_collection.dynamic_particles);
    for(int i=0;i<V.Size();i++){T magnitude_squared=V(i).Magnitude_Squared();if(magnitude_squared>max_magnitude_squared){max_magnitude_squared=magnitude_squared;max_index=i;}}
    if(max_index>-1){
        int p=solid_body_collection.deformable_body_collection.dynamic_particles(max_index);T max_magnitude=sqrt(max_magnitude_squared),max_magnitude_global=max_magnitude;
        if(solid_body_collection.deformable_body_collection.mpi_solids) max_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_magnitude_global);
        LOG::cout<<"maximum velocity = "<<max_magnitude_global;
        if(solid_body_collection.deformable_body_collection.mpi_solids) LOG::cout<<", local = "<<max_magnitude;
        LOG::cout<<" ("<<p<<")"<<std::endl;}
    int max_linear_index=-1,max_angular_index=-1;T max_linear_magnitude_squared=-FLT_MAX,max_angular_magnitude_squared=-FLT_MAX;
    for(int i=0;i<solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles.m;i++){const int p=solid_body_collection.rigid_body_collection.dynamic_rigid_body_particles(i);
        T linear_magnitude_squared=rigid_body_particles.twist(p).linear.Magnitude_Squared(),angular_magnitude_squared=rigid_body_particles.twist(p).angular.Magnitude_Squared();
        if(linear_magnitude_squared>max_linear_magnitude_squared){max_linear_magnitude_squared=linear_magnitude_squared;max_linear_index=p;}
        if(angular_magnitude_squared>max_angular_magnitude_squared){max_angular_magnitude_squared=angular_magnitude_squared;max_angular_index=p;}}
    if(max_linear_index>=0){
        T max_linear_magnitude=sqrt(max_linear_magnitude_squared),max_angular_magnitude=sqrt(max_angular_magnitude_squared);
        T max_linear_magnitude_global=max_linear_magnitude,max_angular_magnitude_global=max_angular_magnitude;
        if(solid_body_collection.deformable_body_collection.mpi_solids){
            max_linear_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_linear_magnitude_global);
            max_angular_magnitude_global=solid_body_collection.deformable_body_collection.mpi_solids->Reduce_Max_Global(max_angular_magnitude_global);}
        LOG::cout<<"maximum rigid linear velocity = "<<max_linear_magnitude_global;
        if(solid_body_collection.deformable_body_collection.mpi_solids) LOG::cout<<", local = "<<max_linear_magnitude;
        LOG::cout<<" ("<<max_linear_index<<")\n";
        LOG::cout<<"maximum rigid angular velocity = "<<max_angular_magnitude_global;
        if(solid_body_collection.deformable_body_collection.mpi_solids) LOG::cout<<", local = "<<max_angular_magnitude;
        LOG::cout<<" ("<<max_angular_index<<")"<<std::endl;}
}
//#####################################################################
// Function Diagnostics
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Diagnostics(const T dt,const T time,const int velocity_time,const int position_time,int step,const char* description)
{
    static const char* time_table[]={"n","(n+1/2)","(n+1)","(n+3/2)","(n+2)"};
    solid_body_collection.Print_Energy(time+position_time*(T).5*time,step);
    PHYSBAM_DEBUG_WRITE_SUBSTEP(STRING_UTILITIES::string_sprintf("Finished step %i (%s).  State: x^%s  v^%s.   dt=%f time=%f",step,description,
        time_table[position_time],time_table[velocity_time],dt,time),2,3);
}
//#####################################################################
// Function Update_Positions_And_Apply_Contact_Forces
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Positions_And_Apply_Contact_Forces(const T dt,const T time,const bool use_saved_pairs)
{
    ARTICULATED_RIGID_BODY<TV>* articulated_rigid_body=&solid_body_collection.rigid_body_collection.articulated_rigid_body;

    Euler_Step_Position(dt,time);

    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i));
        if(!rigid_body.is_static) rigid_body.Update_Bounding_Box();}
    for(int i=0;i<solid_body_collection.rigid_body_collection.kinematic_rigid_bodies.m;i++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(solid_body_collection.rigid_body_collection.kinematic_rigid_bodies(i));
        rigid_body.Update_Bounding_Box();}

    if(solids_parameters.use_rigid_deformable_contact){
        rigid_deformable_collisions->Process_Contact(dt,time,articulated_rigid_body,use_saved_pairs,use_existing_contact,rigid_frame_save,rigid_velocity_difference,
            rigid_angular_momentum_difference,X_save,solids_parameters.rigid_body_collision_parameters.collision_body_thickness);
        // TODO: rigid/deformable shock propagation step
        if(solids_parameters.rigid_body_collision_parameters.perform_contact && solids_parameters.rigid_body_collision_parameters.use_shock_propagation)
            rigid_body_collisions->Shock_Propagation_Using_Graph(dt,time,articulated_rigid_body,use_saved_pairs);}
    else{
        // TODO: move into rigid/deformable collisions and get interleaved with rigid/deformable contact
        if(solids_parameters.rigid_body_collision_parameters.perform_contact)
            rigid_body_collisions->Process_Contact_Using_Graph(dt,time,articulated_rigid_body,solids_parameters.rigid_body_evolution_parameters.correct_contact_energy,use_saved_pairs);

        // rigid/rigid shock propagation
        if(solids_parameters.rigid_body_collision_parameters.perform_contact && solids_parameters.rigid_body_collision_parameters.use_shock_propagation)
            rigid_body_collisions->Shock_Propagation_Using_Graph(dt,time,articulated_rigid_body,use_saved_pairs);

        if(solids_parameters.deformable_object_collision_parameters.perform_collision_body_collisions){
            int interactions=0;
            if(use_existing_contact)
                interactions+=solid_body_collection.deformable_body_collection.collisions.Adjust_Existing_Nodes_For_Collision_Body_Collisions(0);
            else 
                interactions+=solid_body_collection.deformable_body_collection.collisions.Adjust_Nodes_For_Collision_Body_Collisions(solid_body_collection.deformable_body_collection.binding_list,solid_body_collection.deformable_body_collection.soft_bindings,X_save,dt,0);
            if(interactions) LOG::Stat("collision body collisions",interactions);}}

}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time,const int p)
{
    RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
    if(rigid_body.is_static) return;
    if(solid_body_collection.rigid_body_collection.rigid_body_particles.kinematic(p)) kinematic_evolution.Set_External_Velocities(rigid_body.Twist(),time+dt,p);
    rigid_body.Twist().linear+=rigid_velocity_difference(p);
    rigid_body.Angular_Momentum()+=rigid_angular_momentum_difference(p);
    rigid_body.Update_Angular_Velocity();
}
//#####################################################################
// Function Update_Velocity_Using_Stored_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Update_Velocity_Using_Stored_Differences(const T dt,const T time)
{
    for(int i=0;i<solid_body_collection.deformable_body_collection.simulated_particles.m;i++){int p=solid_body_collection.deformable_body_collection.simulated_particles(i);
        solid_body_collection.deformable_body_collection.particles.V(p)+=V_difference(p);}
    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        Update_Velocity_Using_Stored_Differences(dt,time,p);}
    for(int i=0;i<solid_body_collection.rigid_body_collection.kinematic_rigid_bodies.m;i++){int p=solid_body_collection.rigid_body_collection.kinematic_rigid_bodies(i);
        kinematic_evolution.Set_External_Velocities(solid_body_collection.rigid_body_collection.rigid_body_particles.twist(p),time+dt,p);}
}
//#####################################################################
// Function Compute_Momentum_Differences
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Compute_Momentum_Differences()
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    rigid_velocity_difference.Resize(rigid_body_collection.rigid_body_particles.Size());rigid_angular_momentum_difference.Resize(rigid_body_collection.rigid_body_particles.Size());
    V_difference.Resize(particles.Size());
    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        RIGID_BODY<TV>& rigid_body=rigid_body_collection.Rigid_Body(p);
        rigid_velocity_difference(p)=rigid_body.Twist().linear-rigid_velocity_save(p).linear;
        rigid_angular_momentum_difference(p)=rigid_body.Angular_Momentum()-rigid_angular_momentum_save(p);}
    for(int i=0;i<solid_body_collection.deformable_body_collection.simulated_particles.m;i++){
        int p=solid_body_collection.deformable_body_collection.simulated_particles(i);V_difference(p)=particles.V(p)-V_save(p);}
}
//#####################################################################
// Function Save_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Save_Velocity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    V_save.Resize(particles.Size(),false,false);
    rigid_velocity_save.Resize(rigid_body_collection.rigid_body_particles.Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_collection.rigid_body_particles.Size(),false,false);
    V_save.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_velocity_save(p)=rigid_body_collection.rigid_body_particles.twist(p);rigid_angular_momentum_save(p)=rigid_body_collection.rigid_body_particles.angular_momentum(p);}
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){rigid_velocity_save(i)=rigid_body_collection.rigid_body_particles.twist(i);rigid_angular_momentum_save(i)=rigid_body_collection.rigid_body_particles.angular_momentum(i);}}
}
//#####################################################################
// Function Restore_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Restore_Velocity() const
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    particles.V.Subset(solid_body_collection.deformable_body_collection.simulated_particles)=V_save.Subset(solid_body_collection.deformable_body_collection.simulated_particles);
    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        rigid_body_collection.rigid_body_particles.twist(p).linear=rigid_velocity_save(p).linear;
        rigid_body_collection.rigid_body_particles.angular_momentum(p)=rigid_angular_momentum_save(p);rigid_body_collection.Rigid_Body(p).Update_Angular_Velocity();}
    for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            rigid_body_collection.rigid_body_particles.twist(i).linear=rigid_velocity_save(i).linear;rigid_body_collection.rigid_body_particles.angular_momentum(i)=rigid_angular_momentum_save(i);body.Update_Angular_Velocity();}}
}
//#####################################################################
// Function Exchange_Velocity
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Exchange_Velocity()
{
    DEFORMABLE_BODY_COLLECTION<TV>& deformable_body_collection=solid_body_collection.deformable_body_collection;
    DEFORMABLE_PARTICLES<TV>& particles=deformable_body_collection.particles;RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;
    V_save.Resize(particles.Size(),false,false);
    rigid_velocity_save.Resize(rigid_body_particles.Size(),false,false);
    rigid_angular_momentum_save.Resize(rigid_body_particles.Size(),false,false);
    for(int i=0;i<solid_body_collection.deformable_body_collection.simulated_particles.m;i++){int p=solid_body_collection.deformable_body_collection.simulated_particles(i);
        exchange(V_save(p),deformable_body_collection.particles.V(p));}
    for(int i=0;i<solid_body_collection.rigid_body_collection.simulated_rigid_body_particles.m;i++){int p=solid_body_collection.rigid_body_collection.simulated_rigid_body_particles(i);
        exchange(rigid_velocity_save(p),rigid_body_particles.twist(p));
        exchange(rigid_angular_momentum_save(p),rigid_body_particles.angular_momentum(p));}
    for(int i=0;i<rigid_body_particles.Size();i++) if(solid_body_collection.rigid_body_collection.Is_Active(i)){RIGID_BODY<TV>& body=solid_body_collection.rigid_body_collection.Rigid_Body(i);
        if(!body.Is_Simulated()){
            exchange(rigid_velocity_save(i),rigid_body_particles.twist(i));
            exchange(rigid_angular_momentum_save(i),rigid_body_particles.angular_momentum(i));}}
}
//#####################################################################
// Function Initialize_Rigid_Bodies
//#####################################################################
template<class TV> void NEWMARK_EVOLUTION<TV>::
Initialize_Rigid_Bodies(const T frame_rate, const bool restart)
{
    RIGID_BODY_COLLECTION<TV>& rigid_body_collection=solid_body_collection.rigid_body_collection;
    // initialize kinematic object positions and velocities
    if(!restart){
        kinematic_evolution.Get_Current_Kinematic_Keyframes(1/frame_rate,time);
        kinematic_evolution.Set_External_Positions(rigid_body_collection.rigid_body_particles.frame,time);
        kinematic_evolution.Set_External_Velocities(rigid_body_collection.rigid_body_particles.twist,time,time);
        rigid_body_collection.Update_Angular_Momentum();
        for(int i=0;i<rigid_body_collection.rigid_body_particles.Size();i++) if(rigid_body_collection.Is_Active(i)){rigid_body_collection.rigid_body_particles.frame(i).r.Normalize();}}

    RIGID_BODY_COLLISIONS<TV>::Adjust_Bounding_Boxes(rigid_body_collection);
    // rigid body collisions
    if(!rigid_body_collisions)
        rigid_body_collisions=new RIGID_BODY_COLLISIONS<TV>(rigid_body_collection,solids_parameters.rigid_body_collision_parameters,rigids_evolution_callbacks,
            *solid_body_collection.example_forces_and_velocities);
    rigid_body_collisions->spatial_partition->Compute_Voxel_Size(solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_heuristic,
        solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_number_of_cells,solids_parameters.rigid_body_collision_parameters.rigid_collisions_spatial_partition_voxel_size_scale_factor);
    rigid_body_collisions->verbose=solids_parameters.verbose;

    // partitions and hierarchies
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition){
        rigid_body_collisions->intersections.Use_Particle_Partition(true,solids_parameters.rigid_body_collision_parameters.rigid_collisions_particle_partition_size*VECTOR<int,d>::All_Ones_Vector());
        rigid_body_collisions->intersections.Use_Particle_Partition_Center_Phi_Test(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_particle_partition_center_phi_test);}
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy){
        rigid_body_collisions->intersections.Use_Triangle_Hierarchy();
        if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_triangle_hierarchy_center_phi_test) rigid_body_collisions->intersections.Use_Triangle_Hierarchy_Center_Phi_Test();
        if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_use_edge_intersection) rigid_body_collisions->intersections.Use_Edge_Intersection();}

    // rigid deformable collisions
    if(!rigid_deformable_collisions)
        rigid_deformable_collisions=new RIGID_DEFORMABLE_COLLISIONS<TV>(solid_body_collection,*rigid_body_collisions,solids_parameters);

    // dynamics
    solids_parameters.rigid_body_evolution_parameters.rigid_cfl=solids_parameters.cfl;
    solids_parameters.rigid_body_evolution_parameters.rigid_minimum_dt=solids_parameters.rigid_body_evolution_parameters.minimum_rigid_body_time_step_fraction/frame_rate;
    solids_parameters.rigid_body_evolution_parameters.rigid_maximum_dt=solids_parameters.rigid_body_evolution_parameters.maximum_rigid_body_time_step_fraction/frame_rate;
    solids_evolution_callbacks->Update_Solids_Parameters(time);
    rigid_body_collisions->Initialize_Data_Structures(); // Must be done before we check interpenetration statistics

    // check for bad initial data (need to give it a chance to set up the collision manager first)
    if(solids_parameters.rigid_body_collision_parameters.rigid_collisions_print_interpenetration_statistics) rigid_body_collisions->Print_Interpenetration_Statistics();
    else if(!rigid_body_collisions->Check_For_Any_Interpenetration()) LOG::cout<<"No initial interpenetration"<<std::endl;
}
//#####################################################################
// Function Use_CFL
//#####################################################################
template<class TV> bool NEWMARK_EVOLUTION<TV>::
Use_CFL() const
{
    return true;
}
//#####################################################################
namespace PhysBAM{
template class NEWMARK_EVOLUTION<VECTOR<float,1> >;
template class NEWMARK_EVOLUTION<VECTOR<float,2> >;
template class NEWMARK_EVOLUTION<VECTOR<float,3> >;
template class NEWMARK_EVOLUTION<VECTOR<double,1> >;
template class NEWMARK_EVOLUTION<VECTOR<double,2> >;
template class NEWMARK_EVOLUTION<VECTOR<double,3> >;
}
