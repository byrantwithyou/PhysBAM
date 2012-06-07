//#####################################################################
// Copyright 2006-2008, Geoffrey Irving, Avi Robinson-Mosher, Craig Schroeder, Andrew Selle, Tamar Shinar, Joseph Teran.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_EVOLUTION
//#####################################################################
#include <PhysBAM_Tools/Arrays/ARRAY_VIEW.h>
#include <PhysBAM_Tools/Arrays/INDIRECT_ARRAY.h>
#include <PhysBAM_Tools/Arrays_Computations/MAGNITUDE.h>
#include <PhysBAM_Tools/Krylov_Solvers/CONJUGATE_GRADIENT.h>
#include <PhysBAM_Tools/Krylov_Solvers/IMPLICIT_SOLVE_PARAMETERS.h>
#include <PhysBAM_Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <PhysBAM_Tools/Log/LOG.h>
#include <PhysBAM_Tools/Vectors/TWIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Bindings/BINDING_LIST.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids/SOLIDS_PARAMETERS.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/BACKWARD_EULER_EVOLUTION.h>
#include <PhysBAM_Solids/PhysBAM_Solids/Solids_Evolution/GENERALIZED_VELOCITY.h>
namespace PhysBAM{
//#####################################################################
// Class BACKWARD_EULER_SYSTEM
//#####################################################################
namespace{
// assumes all solids_forces are linear in velocity, with a symmetric negative definite Jacobian.
template<class TV>
class BACKWARD_EULER_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    typedef INDIRECT_ARRAY<ARRAY_VIEW<TV> > VECTOR_T;
    typedef KRYLOV_VECTOR_WRAPPER<T,VECTOR_T> KRYLOV_VECTOR_T;

    const SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    T dt,time;
    INDIRECT_ARRAY<const ARRAY_VIEW<const T> > mass,one_over_mass;
    mutable ARRAY<TWIST<TV> > rigid_dV,rigid_dF; // TODO: support rigid body particles and get rid of these

    BACKWARD_EULER_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection_input,const T dt_input,const T time_input)
        :KRYLOV_SYSTEM_BASE<typename TV::SCALAR>(false,true),solid_body_collection(solid_body_collection_input),dt(dt_input),
        time(time_input),mass(solid_body_collection.deformable_body_collection.particles.mass,solid_body_collection.deformable_body_collection.dynamic_particles),
        one_over_mass(solid_body_collection.deformable_body_collection.particles.one_over_mass,solid_body_collection.deformable_body_collection.dynamic_particles)
    {}

    void Multiply(const KRYLOV_VECTOR_BASE<T>& bdV,KRYLOV_VECTOR_BASE<T>& bdF) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dV=debug_cast<const KRYLOV_VECTOR_T&>(bdV);KRYLOV_VECTOR_T& dF=debug_cast<KRYLOV_VECTOR_T&>(bdF);
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(dV.v.array);
    solid_body_collection.Force_Differential(dV.v.array,dF.v.array,time);
    dF.v*=dt; // include embedded particles
    solid_body_collection.Add_Velocity_Dependent_Forces(dV.v.array,rigid_dV,dF.v.array,rigid_dF,time);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(dF.v.array);
    for(int p=0;p<dV.v.Size();p++) dF.v(p)=dV.v(p)-dt*one_over_mass(p)*dF.v(p);}

    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& dV) const PHYSBAM_OVERRIDE
    {Project(dV);}

    void Project(KRYLOV_VECTOR_BASE<T>& bdV) const PHYSBAM_OVERRIDE
    {KRYLOV_VECTOR_T& dV=debug_cast<KRYLOV_VECTOR_T&>(bdV);
    solid_body_collection.example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(dV.v.array,time,time);}

    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& bdV1,const KRYLOV_VECTOR_BASE<T>& bdV2) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dV1=debug_cast<const KRYLOV_VECTOR_T&>(bdV1),&dV2=debug_cast<const KRYLOV_VECTOR_T&>(bdV2);
    return dV1.v.Inner_Product_Double_Precision(mass,dV2.v);}

    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& bdV) const PHYSBAM_OVERRIDE
    {const KRYLOV_VECTOR_T& dV=debug_cast<const KRYLOV_VECTOR_T&>(bdV);return dV.v.Maximum_Magnitude();}

    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& dV) const PHYSBAM_OVERRIDE {}
};
}
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BACKWARD_EULER_EVOLUTION<TV>::
BACKWARD_EULER_EVOLUTION(SOLIDS_PARAMETERS<TV>& solids_parameters_input,SOLID_BODY_COLLECTION<TV>& solid_body_collection_input)
    :SOLIDS_EVOLUTION<TV>(solids_parameters_input,solid_body_collection_input)
{}
//#####################################################################
// Function One_Newton_Step_Backward_Euler
//#####################################################################
template<class TV> void BACKWARD_EULER_EVOLUTION<TV>::
One_Newton_Step_Backward_Euler(const T dt,const T time,ARRAY_VIEW<const TV> V_save,ARRAY<TV>& dV_full)
{
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > KRYLOV_VECTOR_T;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;

    dV_full.Resize(particles.Size()); // an initial guess might be passed in for dV, otherwise it's zero
    B_full.Resize(particles.Size(),false,false);
    INDIRECT_ARRAY<ARRAY_VIEW<TV> > V(particles.V,solid_body_collection.deformable_body_collection.dynamic_particles);
    KRYLOV_VECTOR_T dV(dV_full,solid_body_collection.deformable_body_collection.dynamic_particles),
        B(B_full,solid_body_collection.deformable_body_collection.dynamic_particles);
    KRYLOV_SOLVER<T>::Ensure_Size(krylov_vectors,dV,1);
    KRYLOV_VECTOR_T& F=debug_cast<KRYLOV_VECTOR_T&>(*krylov_vectors(0));
    INDIRECT_ARRAY<ARRAY_VIEW<T> > one_over_mass(particles.one_over_mass,solid_body_collection.deformable_body_collection.dynamic_particles);
    ARRAY<TWIST<TV> > rigid_F_full;

    F.v.Fill(TV());
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities();
    solid_body_collection.Add_All_Forces(F.v.array,rigid_F_full,time);example_forces_and_velocities.Add_External_Forces(F.v.array,time);
    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(F.v.array);
    for(int p=0;p<B.v.Size();p++) B.v.array(p)=V_save(p)-V.array(p)+dt*one_over_mass(p)*F.v.array(p);

    BACKWARD_EULER_SYSTEM<TV> system(solid_body_collection,dt,time);
    CONJUGATE_GRADIENT<T> cg;
    cg.restart_iterations=solids_parameters.implicit_solve_parameters.cg_restart_iterations;
    cg.print_diagnostics=solid_body_collection.print_diagnostics;
    cg.print_residuals=solid_body_collection.print_residuals;
    cg.iterations_used=&solid_body_collection.iterations_used_diagnostic;
    if(!cg.Solve(system,dV,B,krylov_vectors,solids_parameters.implicit_solve_parameters.cg_tolerance,0,
            solids_parameters.implicit_solve_parameters.cg_iterations) && !solids_parameters.use_partially_converged_result)
        throw std::runtime_error("Backward euler conjugate gradient failed");
}
//#####################################################################
// Function Advance_One_Time_Step
//#####################################################################
template<class TV> void BACKWARD_EULER_EVOLUTION<TV>::
Advance_One_Time_Step_Velocity(const T dt,const T time,const bool solids) // TODO: Split this into position/velocity parts
{
    typedef KRYLOV_VECTOR_WRAPPER<T,INDIRECT_ARRAY<ARRAY_VIEW<TV> > > KRYLOV_VECTOR_T;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>& example_forces_and_velocities=*solid_body_collection.example_forces_and_velocities;
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    BINDING_LIST<TV>& binding_list=solid_body_collection.deformable_body_collection.binding_list;

    solid_body_collection.Enforce_Definiteness(true);example_forces_and_velocities.Update_Time_Varying_Material_Properties(time+dt);
    example_forces_and_velocities.Set_External_Velocities(particles.V,time+dt,time+dt);
    Euler_Step_Position(dt,time); // initial guess
    solid_body_collection.Update_Position_Based_State(time+dt,false);

    ARRAY<TV> V_save(particles.V.Subset(solid_body_collection.deformable_body_collection.dynamic_particles));
    ARRAY<TWIST<TV> > rigid_R_full;

    // iterate to steady state
    T supnorm=0;int iteration;
    dV_full.Resize(particles.Size(),false,false);
    for(iteration=0;iteration<solids_parameters.newton_iterations;iteration++){
        dV_full.Subset(solid_body_collection.deformable_body_collection.dynamic_particles).Fill(TV());
        One_Newton_Step_Backward_Euler(dt,time+dt,V_save,dV_full);
        INDIRECT_ARRAY<ARRAY_VIEW<TV> > X(particles.X,solid_body_collection.deformable_body_collection.dynamic_particles),
            V(particles.V,solid_body_collection.deformable_body_collection.dynamic_particles),
            dV(dV_full,solid_body_collection.deformable_body_collection.dynamic_particles),
            R=debug_cast<KRYLOV_VECTOR_T&>(*krylov_vectors(0)).v;
        V+=dV;
        for(int p=0;p<X.Size();p++) X(p)+=dt*dV(p);
        binding_list.Clamp_Particles_To_Embedded_Positions();
        solid_body_collection.Update_Position_Based_State(time+dt,false);
        // compute residual
        R.array.Subset(solid_body_collection.deformable_body_collection.dynamic_particles).Fill(TV());

        solid_body_collection.Add_All_Forces(R.array,rigid_R_full,time+dt);
        example_forces_and_velocities.Add_External_Forces(R.array,time+dt);
        binding_list.Distribute_Force_To_Parents(R.array);
        INDIRECT_ARRAY<ARRAY_VIEW<T> > one_over_mass(particles.one_over_mass,solid_body_collection.deformable_body_collection.dynamic_particles);
        for(int p=0;p<R.Size();p++) R(p)=V_save(p)-V(p)+dt*one_over_mass(p)*R(p);
        example_forces_and_velocities.Zero_Out_Enslaved_Velocity_Nodes(R.array,time+dt,time+dt);
        supnorm=R.Maximum_Magnitude();
        if(solid_body_collection.print_residuals) LOG::cout<<"Newton iteration residual after "<<iteration+1<<" iterations = "<<supnorm<<std::endl;
        if(supnorm<=solids_parameters.newton_tolerance && solid_body_collection.print_diagnostics){
            LOG::cout<<"Newton converged in "<<iteration+1<<std::endl;break;}}
    if(iteration>=solids_parameters.newton_iterations && solid_body_collection.print_diagnostics)
        LOG::cout<<"Newton iteration did not converge in "<<solids_parameters.newton_iterations<<", error = "<<supnorm<<" > "<<solids_parameters.newton_tolerance<<std::endl;
    LOG::cout<<"maximum velocity = "<<particles.V.Subset(solid_body_collection.deformable_body_collection.dynamic_particles).Maximum_Magnitude()<<std::endl;
}
//#####################################################################
template class BACKWARD_EULER_EVOLUTION<VECTOR<float,2> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<float,3> >;
#ifndef COMPILE_WITHOUT_DOUBLE_SUPPORT
template class BACKWARD_EULER_EVOLUTION<VECTOR<double,2> >;
template class BACKWARD_EULER_EVOLUTION<VECTOR<double,3> >;
#endif
}
