//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar, Chenfanfu Jiang, Chuyuan Fu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Read_Write/OCTAVE_OUTPUT.h>
#include <Geometry/Implicit_Objects/IMPLICIT_OBJECT.h>
#include <Hybrid_Methods/Examples_And_Drivers/MPM_EXAMPLE.h>
#include <Hybrid_Methods/Iterators/PARTICLE_GRID_WEIGHTS.h>
#include <Hybrid_Methods/System/KKT_KRYLOV_SYSTEM.h>
#include <Hybrid_Methods/System/KKT_KRYLOV_VECTOR.h>
namespace PhysBAM{
//#####################################################################
// Constructor
//#####################################################################
template<class TV> KKT_KRYLOV_SYSTEM<TV>::
KKT_KRYLOV_SYSTEM(const MPM_EXAMPLE<TV>& example)
    :KRYLOV_SYSTEM_BASE<T>(false,false),example(example)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> KKT_KRYLOV_SYSTEM<TV>::
~KKT_KRYLOV_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void KKT_KRYLOV_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const
{
    PHYSBAM_ASSERT(&BV!=&BF);
    const KKT_KRYLOV_VECTOR<TV>& V=debug_cast<const KKT_KRYLOV_VECTOR<TV>&>(BV);
    KKT_KRYLOV_VECTOR<TV>& F=debug_cast<KKT_KRYLOV_VECTOR<TV>&>(BF);
    F.u.Fill(TV());F.p.Fill((T)0);
    // -DT^T*p
    T dt=example.dt;
    const ARRAY<T,TV_INT>& r=example.density;
    int index=example.DT.offsets(0);
    for(int i=0;i<example.DT.m;i++){
        int end=example.DT.offsets(i+1);T y=V.p.array(example.valid_pressure_indices(i))*sqrt(dt);
        for(;index<end;index++){
            int cell=example.DT.A(index).j/TV::m;
            if(abs(sqrt(r.array(cell))>1e-16)) F.u.array.Flattened()(example.DT.A(index).j)-=example.DT.A(index).a*y/sqrt(r.array(cell));}}
    // (1/dt)*R*u
    T one_over_dt=(T)1/example.dt;
    for(int t=0;t<example.valid_pressure_indices.m;t++){
        int cell=example.valid_pressure_indices(t);
        F.u.array(cell)+=V.u.array(cell);}
    // -DT*u-(1/lambda*dt)*J.Inverse()*p
    index=example.DT.offsets(0);
    for(int i=0;i<example.DT.m;i++){
        int end=example.DT.offsets(i+1);T sum=(T)0;
        for(;index<end;index++){
            int cell=example.DT.A(index).j/TV::m;
            if(abs(r.array(cell))>1e-16) sum-=example.DT.A(index).a*V.u.array.Flattened()(example.DT.A(index).j)/sqrt(r.array(cell));}
        F.p.array(example.valid_pressure_indices(i))=sum*sqrt(dt);}
    for(int t=0;t<example.valid_pressure_indices.m;t++){
        TV_INT valid_cell=example.valid_pressure_cell_indices(t); 
        F.p(valid_cell)-=one_over_dt*example.one_over_lambda(valid_cell)/example.J(valid_cell)*V.p(valid_cell);}
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double KKT_KRYLOV_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const KKT_KRYLOV_VECTOR<TV>& X=debug_cast<const KKT_KRYLOV_VECTOR<TV>&>(x);
    const KKT_KRYLOV_VECTOR<TV>& Y=debug_cast<const KKT_KRYLOV_VECTOR<TV>&>(y);
    T r=0;
#pragma omp parallel for reduction(+:r)
    for(int k=0;k<example.valid_pressure_indices.m;k++){
        int i=example.valid_pressure_indices(k);
        r+=X.p.array(i)*Y.p.array(i);
        r+=X.u.array(i).Dot(Y.u.array(i));}
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR KKT_KRYLOV_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const
{
    return sqrt(Inner_Product(BR,BR));
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void KKT_KRYLOV_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& BV) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void KKT_KRYLOV_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void KKT_KRYLOV_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BR) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void KKT_KRYLOV_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const
{
}
template class KKT_KRYLOV_SYSTEM<VECTOR<float,2> >;
template class KKT_KRYLOV_SYSTEM<VECTOR<float,3> >;
template class KKT_KRYLOV_SYSTEM<VECTOR<double,2> >;
template class KKT_KRYLOV_SYSTEM<VECTOR<double,3> >;
}
