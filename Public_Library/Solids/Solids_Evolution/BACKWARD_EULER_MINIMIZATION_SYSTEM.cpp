//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY.h>
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Bindings/BINDING_LIST.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Forces_And_Torques/EXAMPLE_FORCES_AND_VELOCITIES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
BACKWARD_EULER_MINIMIZATION_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities)
    :KRYLOV_SYSTEM_BASE<T>(false,false),solid_body_collection(solid_body_collection),dt(0),time(0),tmp(0),forced_collisions(),example_forces_and_velocities(example_forces_and_velocities)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
~BACKWARD_EULER_MINIMIZATION_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const GENERALIZED_VELOCITY<TV>& V=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV);
    GENERALIZED_VELOCITY<TV>& F=debug_cast<GENERALIZED_VELOCITY<TV>&>(BF);
    GENERALIZED_VELOCITY<TV>& t=debug_cast<GENERALIZED_VELOCITY<TV>&>(*tmp);
    DEFORMABLE_PARTICLES<TV>& particles=solid_body_collection.deformable_body_collection.particles;
    RIGID_BODY_PARTICLES<TV>& rigid_body_particles=solid_body_collection.rigid_body_collection.rigid_body_particles;

    t=V;
    solid_body_collection.deformable_body_collection.binding_list.Clamp_Particles_To_Embedded_Velocities(t.V.array);
    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        t.V.array(c.p).Project_Orthogonal_To_Unit_Direction(c.n);}

    F.V.array.Fill(TV());
    F.rigid_V.array.Fill(TWIST<TV>());
    solid_body_collection.Add_Implicit_Velocity_Independent_Forces(t.V.array,t.rigid_V.array,F.V.array,F.rigid_V.array,-dt*dt,time);

    for(int p=0;p<particles.number;p++) F.V.array(p)+=particles.mass(p)*t.V.array(p);

    for(int p=0;p<rigid_body_particles.number;p++){
        RIGID_BODY<TV>& rigid_body=solid_body_collection.rigid_body_collection.Rigid_Body(p);
        if(rigid_body.Has_Infinite_Inertia()) continue;
        F.rigid_V.array(p)+=rigid_body.Inertia_Times(t.rigid_V.array(p));}

    for(int i=0;i<collisions.m;i++){
        const COLLISION& c=collisions(i);
        TV& v=F.V.array(c.p),tt=t.V.array(c.p);
        v-=c.n*c.n.Dot(v)+dt*c.n_dE*c.H*tt+dt*c.n.Dot(tt)*c.H_dE+dt*c.H_dE.Dot(tt)*c.n;}

    solid_body_collection.deformable_body_collection.binding_list.Distribute_Force_To_Parents(F.V.array);
    solid_body_collection.deformable_body_collection.binding_list.Clear_Hard_Bound_Particles(F.V.array);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const
{
    const GENERALIZED_VELOCITY<TV>& V0=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV0),&V1=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV1);
    return V0.V.Dot_Double_Precision(V1.V)+V0.rigid_V.Dot_Double_Precision(V1.rigid_V);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
    GENERALIZED_VELOCITY<TV>& V=debug_cast<GENERALIZED_VELOCITY<TV>&>(BV);
    if(example_forces_and_velocities) example_forces_and_velocities->Zero_Out_Enslaved_Velocity_Nodes(V.V.array,time,time);
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<float,1> >;
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<float,2> >;
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<float,3> >;
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<double,1> >;
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<double,2> >;
template class BACKWARD_EULER_MINIMIZATION_SYSTEM<VECTOR<double,3> >;
}
