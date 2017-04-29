//#####################################################################
// Copyright 2016, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Log/FINE_TIMER.h>
#include <Core/Log/LOG.h>
#include <Core/Log/SCOPE.h>
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Core/Random_Numbers/RANDOM_NUMBERS.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_SYSTEM.h>
#include <Hybrid_Methods/Projection/MPM_PROJECTION_VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class TV> MPM_PROJECTION_SYSTEM<TV>::
MPM_PROJECTION_SYSTEM()
    :KRYLOV_SYSTEM_BASE<T>(true,false),dc_present(false)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class TV> MPM_PROJECTION_SYSTEM<TV>::
~MPM_PROJECTION_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    TIMER_SCOPE_FUNC;
    const MPM_PROJECTION_VECTOR<TV>& vx=debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(x);
    MPM_PROJECTION_VECTOR<TV>& vresult=debug_cast<MPM_PROJECTION_VECTOR<TV>&>(result);
    A.Times_Threaded(vx.v,vresult.v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_PROJECTION_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    TIMER_SCOPE_FUNC;
    const MPM_PROJECTION_VECTOR<TV>& vx=debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(x),vy=debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(y);
    double r=0;
#pragma omp parallel for reduction(+:r)
    for(int i=0;i<vx.v.m;i++)
        r+=vx.v(i)*vy.v(i);
    return r;
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_PROJECTION_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    TIMER_SCOPE_FUNC;
    const MPM_PROJECTION_VECTOR<TV>& vx=debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(x);
    T r=0;
#pragma omp parallel for reduction(max:r)
    for(int i=0;i<vx.v.m;i++)
        r=std::max(r,std::abs(vx.v(i)));
    return r;
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
    Project_Nullspace(x);
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Compute_Ones_Nullspace
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Compute_Ones_Nullspace()
{
    null_u.v.Resize(A.m);
    null_u.v.Fill(1/sqrt(A.m));
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
    if(!dc_present){
        MPM_PROJECTION_VECTOR<TV>& v=debug_cast<MPM_PROJECTION_VECTOR<TV>&>(x);
        v.Copy(-Inner_Product(v,null_u),null_u,v);}
} 
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    const MPM_PROJECTION_VECTOR<TV>& vr=debug_cast<const MPM_PROJECTION_VECTOR<TV>&>(r);
    MPM_PROJECTION_VECTOR<TV>& vz=debug_cast<MPM_PROJECTION_VECTOR<TV>&>(z);
    temp_vector.Resize(A.m);
    A.C->Solve_Forward_Substitution(vr.v,temp_vector,true);
    A.C->Solve_Backward_Substitution(temp_vector,vz.v,false,true);
}
//#####################################################################
// Function Test_Helper
//#####################################################################
template<class T> static void
Test_Helper(const ARRAY<T>& x,const ARRAY<T>& y,T tolerance,const char* msg)
{
    T a=x.Magnitude();
    T b=y.Magnitude();
    T c=(x-y).Magnitude();
    T r=c/std::max(std::max(a,b),(T)1e-30);
    if(r>=tolerance)
        LOG::printf("FAIL: %s: %g %g %g  rel %g\n",msg,a,b,c,r);
}
//#####################################################################
// Function Test
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Test() const
{
    T tolerance=(T)1e-5;
    PHYSBAM_ASSERT(A.m==A.n && A.m==gradient.n && mass.m==gradient.m);
    ARRAY<T> x(A.m),y(A.m),z(A.m),w(gradient.m);
    RANDOM_NUMBERS<T> random;
    random.Fill_Uniform(x,-1,1);
    gradient.Times(x,w);
    w/=mass;
    gradient.Transpose_Times(w,y);
    A.Times(x,z);
    Test_Helper(y,z,tolerance,"A = G^T M^(-1) G");
    MPM_PROJECTION_VECTOR<TV> v;
    v.v=z;
    Project(v);
    Test_Helper(v.v,z,tolerance,"P A = A");
}
namespace PhysBAM{
template class MPM_PROJECTION_SYSTEM<VECTOR<float,3> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<float,2> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<float,1> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,3> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,2> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,1> >;
}
