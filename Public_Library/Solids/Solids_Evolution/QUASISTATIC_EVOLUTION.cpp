//#####################################################################
// Copyright 2006-2008, Ron Fedkiw, Geoffrey Irving, Craig Schroeder, Tamar Shinar, Eftychios Sifakis, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class QUASISTATIC_EVOLUTION
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Arrays/INDIRECT_ARRAY.h>
#include <Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Log/LOG.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Parallel_Computation/MPI_SOLIDS.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids/SOLIDS_PARAMETERS.h>
#include <Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
#include <Solids/Solids_Evolution/QUASISTATIC_EVOLUTION.h>
using namespace PhysBAM;
//#####################################################################
// Class QUASISTATIC_SYSTEM
//#####################################################################
namespace{
template<class TV>
class QUASISTATIC_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    typedef INDIRECT_ARRAY<ARRAY<TV> > VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<TV> > > KRYLOV_VECTOR_T;

    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities;
    T time;
    MPI_SOLIDS<TV>* mpi_solids;

    QUASISTATIC_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities,const T time_input,MPI_SOLIDS<TV>* mpi_solids_input)
        :KRYLOV_SYSTEM_BASE<T>(false,true),solid_body_collection(solid_body_collection_input),example_forces_and_velocities(example_forces_and_velocities),time(time_input),mpi_solids(mpi_solids_input)
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bdX,KRYLOV_VECTOR_BASE<T>& bdF) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dX=debug_cast<const KRYLOV_VECTOR_T&>(bdX);KRYLOV_VECTOR_T& dF=debug_cast<KRYLOV_VECTOR_T&>(bdF);
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(dX.v.array);
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Positions(dX.v.array); // TODO: assumes bindings are linear, consider switching this to clamping velocities
    if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(dX.v.array);
    solid_body_collection.Force_Differential(dX.v.array,dF.v.array,time);
    if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(dX.v.array);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(dF.v.array);}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& dX) const PHYSBAM_OVERRIDE
    {Project(dX);}

    void Project(KRYLOV_VECTOR_BASE<T>& bdX) const PHYSBAM_OVERRIDE
    {KRYLOV_VECTOR_T& dX=debug_cast<KRYLOV_VECTOR_T&>(bdX);example_forces_and_velocities.Zero_Out_Enslaved_Position_Nodes(dX.v.array,time);}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bdX0,const KRYLOV_VECTOR_BASE<T>& bdX2) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dX0=debug_cast<const KRYLOV_VECTOR_T&>(bdX0),&dX2=debug_cast<const KRYLOV_VECTOR_T&>(bdX2);
    T inner_product=dX0.v.Dot(dX2.v);
    if(mpi_solids) inner_product=mpi_solids->Reduce_Add(inner_product);
    return inner_product;}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bdX) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dX=debug_cast<const KRYLOV_VECTOR_T&>(bdX);
    T convergence_norm=dX.v.Maximum_Magnitude();
    if(mpi_solids) convergence_norm=mpi_solids->Reduce_Max(convergence_norm);
    return convergence_norm;}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& dX) const PHYSBAM_OVERRIDE {}
};
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> QUASISTATIC_EVOLUTION<TV>::
QUASISTATIC_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities_input)
    :SOLIDS_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input,example_forces_and_velocities_input),balance_external_forces_only(false)
{}
//#####################################################################
// Function One_Newton_Step_Toward_Steady_State
//#####################################################################
template<class TV> void QUASISTATIC_EVOLUTION<TV>::
One_Newton_Step_Toward_Steady_State(const T time,ARRAY<TV>& dX_full)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    const ARRAY<int>& dynamic_particles=solid_body_collection.deformable_body_collection.dynamic_particles;

    dX_full.Resize(particles.Size()); // an initial guess might be passed in for dX, otherwise it's zero
    KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY<TV> > > dX(dX_full,dynamic_particles),B(B_full,dynamic_particles);

    B_full.Subset(solid_body_collection.deformable_body_collection.dynamic_particles).Fill(TV());
    if(!balance_external_forces_only) solid_body_collection.Add_Velocity_Independent_Forces(B_full,rigid_B_full,time);
    example_forces_and_velocities.Add_External_Forces(B_full,time);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(B_full);
    B.v=-B.v;

    QUASISTATIC_SYSTEM<TV> system(solid_body_collection,example_forces_and_velocities,time,solid_body_collection.deformable_body_collection.mpi_solids);
    CONJUGATE_GRADIENT<T> cg;cg.restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
    cg.print_diagnostics=solid_body_collection.print_diagnostics;cg.print_residuals=solid_body_collection.print_residuals;
    cg.iterations_used=&solid_body_collection.iterations_used_diagnostic;
    if(!cg.Solve(system,dX,B,krylov_vectors,solids_parameters.implicit_solve_parameters.cg_tolerance,0,solids_parameters.implicit_solve_parameters.cg_iterations) && !solids_parameters.use_partially_converged_result)
        throw std::runtime_error("Quasistatic conjugate gradient failed");
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void QUASISTATIC_EVOLUTION<TV>::
Advance_One_Time_Step_Position(const T dt,const T time,const bool solids)
{
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;
    const ARRAY<int>& dynamic_particles=solid_body_collection.deformable_body_collection.dynamic_particles;
    const ARRAY<int>& simulated_particles=solid_body_collection.deformable_body_collection.simulated_particles;
    MPI_SOLIDS<TV>* mpi_solids=solid_body_collection.deformable_body_collection.mpi_solids;

    // for parallelism, assume ghost data is synchronized

    // prepare for force computation
    solid_body_collection.Enforce_Definiteness(true);
    example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);
    Set_External_Positions(particles.X,time+dt);binding_list.Clamp_Particles_To_Embedded_Positions();
    solid_body_collection.Update_Position_Based_State(time+dt,true);solid_body_collection.deformable_body_collection.Update_Collision_Penalty_Forces_And_Derivatives();

    // iterate to steady state
    T supnorm=0;int iteration;
    dX_full.Resize(particles.Size(),false,false);
    for(iteration=0;iteration<solids_parameters.newton_iterations;iteration++){
        INDIRECT_ARRAY<ARRAY<TV>,ARRAY<int>&> dX_subset=dX_full.Subset(simulated_particles);
        dX_subset.Fill(TV()); // initial guess is zero
        One_Newton_Step_Toward_Steady_State(time+dt,dX_full);
        particles.X.Subset(dynamic_particles)+=dX_subset;
        if(mpi_solids) mpi_solids->Exchange_Binding_Boundary_Data(particles.X);
        binding_list.Clamp_Particles_To_Embedded_Positions();
        if(mpi_solids) mpi_solids->Exchange_Force_Boundary_Data(particles.X);
        solid_body_collection.Update_Position_Based_State(time+dt,true);
        solid_body_collection.deformable_body_collection.Update_Collision_Penalty_Forces_And_Derivatives();
        GENERALIZED_VELOCITY<TV>& R=debug_cast<GENERALIZED_VELOCITY<TV>&>(*krylov_vectors(0));
        R.V.Fill(TV());
        solid_body_collection.Add_Velocity_Independent_Forces(R.V.array,R.rigid_V.array,time+dt);
        example_forces_and_velocities.Add_External_Forces(R.V.array,time+dt);
        binding_list.Distribute_Force_To_Parents(R.V.array);
        example_forces_and_velocities.Zero_Out_Enslaved_Position_Nodes(R.V.array,time+dt);
        supnorm=R.V.Maximum_Magnitude();
        if(mpi_solids) supnorm=mpi_solids->Reduce_Max(supnorm);
        if(solid_body_collection.print_residuals) LOG::cout<<"Newton iteration residual after "<<iteration+1<<" iterations = "<<supnorm<<std::endl;
        if(supnorm<=solids_parameters.newton_tolerance && solid_body_collection.print_diagnostics){
            LOG::cout<<"Newton converged in "<<iteration+1<<std::endl;break;}}
    if(iteration>=solids_parameters.newton_iterations && solid_body_collection.print_diagnostics)
        LOG::cout<<"Newton iteration did not converge in "<<solids_parameters.newton_iterations<<", error = "<<supnorm<<" > "
                 <<solids_parameters.newton_tolerance<<std::endl;
}
namespace PhysBAM{
template class QUASISTATIC_EVOLUTION<VECTOR<float,2> >;
template class QUASISTATIC_EVOLUTION<VECTOR<float,3> >;
template class QUASISTATIC_EVOLUTION<VECTOR<double,2> >;
template class QUASISTATIC_EVOLUTION<VECTOR<double,3> >;
}
