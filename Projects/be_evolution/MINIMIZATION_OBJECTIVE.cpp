//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Rigids/Rigid_Bodies/RIGID_BODY_COLLECTION.h>
#include <Deformables/Deformable_Objects/DEFORMABLE_BODY_COLLECTION.h>
#include <Deformables/Particles/DEFORMABLE_PARTICLES.h>
#include <Solids/Solids/SOLID_BODY_COLLECTION.h>
#include "MINIMIZATION_OBJECTIVE.h"
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MINIMIZATION_OBJECTIVE<TV>::
MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,T dt,T time)
    :KRYLOV_SYSTEM_BASE<T>(false,false),solid_body_collection(solid_body_collection),dt(dt),time(time),
    fake_rigid_x(solid_body_collection.rigid_body_collection.rigid_body_particles.twist.m),
    x1(solid_body_collection.deformable_body_collection.particles.X,fake_rigid_x,solid_body_collection),
    v1(solid_body_collection.deformable_body_collection.particles.V,solid_body_collection.rigid_body_collection.rigid_body_particles.twist,solid_body_collection),
    x0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),v0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),
    d(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),tmp0(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),
    tmp1(static_cast<GENERALIZED_VELOCITY<TV>&>(*x1.Clone_Default())),system(0,solid_body_collection,dt,time,time,0,0,0,true,true)
{
    x0=x1;
    v0=v1;
    a=1/(dt*dt);
    b=-a;
    c=-1/dt;
    d.Copy(b,x0);
    d.Copy(c,v0,d);
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MINIMIZATION_OBJECTIVE<TV>::
~MINIMIZATION_OBJECTIVE()
{
    delete &x0;
    delete &v0;
    delete &d;
    delete &tmp0;
    delete &tmp1;
}
//#####################################################################
// Function Compute
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const
{
    PHYSBAM_ASSERT(!h || h==this);

    tmp0.Copy(a,x,d);

    const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->v1.Copy(-1,x0,x);
    const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->v1*=1/dt;
    const_cast<MINIMIZATION_OBJECTIVE<TV>*>(this)->x1.Copy(1,x);
    solid_body_collection.Update_Position_Based_State(time,true);

    if(e){
        T ke=0,pe=0;
        solid_body_collection.Compute_Energy(time,ke,pe);
        *e=system.Inner_Product(tmp0,tmp0)/(2*a)+pe;}

    if(g){
        tmp1*=0;
        solid_body_collection.Add_Velocity_Independent_Forces(tmp1.V.array,tmp1.rigid_V.array,time);
        system.projection_data.mass.Multiply(tmp0,debug_cast<GENERALIZED_VELOCITY<TV>&>(*g),false);
        *g-=tmp1;}
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    const GENERALIZED_VELOCITY<TV>& V=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV);
    GENERALIZED_VELOCITY<TV>& F=debug_cast<GENERALIZED_VELOCITY<TV>&>(BF);
    solid_body_collection.Implicit_Velocity_Independent_Forces(V.V.array,V.rigid_V.array,F.V.array,F.rigid_V.array,-1,time);
    system.projection_data.mass.Multiply(V,tmp0,false);
    F.Copy(a,tmp0,F);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MINIMIZATION_OBJECTIVE<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const
{
    const GENERALIZED_VELOCITY<TV>& V1=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV1),&V2=debug_cast<const GENERALIZED_VELOCITY<TV>&>(BV2);
    return V1.V.Dot_Double_Precision(V2.V)+V1.rigid_V.Dot_Double_Precision(V2.rigid_V);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MINIMIZATION_OBJECTIVE<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Project(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MINIMIZATION_OBJECTIVE<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class MINIMIZATION_OBJECTIVE<VECTOR<float,2> >;
template class MINIMIZATION_OBJECTIVE<VECTOR<float,3> >;
template class MINIMIZATION_OBJECTIVE<VECTOR<double,2> >;
template class MINIMIZATION_OBJECTIVE<VECTOR<double,3> >;
