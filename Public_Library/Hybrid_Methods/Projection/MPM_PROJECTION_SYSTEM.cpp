//#####################################################################
// Copyright 2016, Lin Huang, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Core/Matrices/MATRIX.h>
#include <Core/Matrices/MATRIX_MXN.h>
#include <Core/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
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
    :KRYLOV_SYSTEM_BASE<T>(true,false)
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
    const MPM_PROJECTION_VECTOR<TV>& vx=dynamic_cast<const MPM_PROJECTION_VECTOR<TV>&>(x);
    MPM_PROJECTION_VECTOR<TV>& vresult=dynamic_cast<MPM_PROJECTION_VECTOR<TV>&>(result);
    A.Times(vx.v,vresult.v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class TV> double MPM_PROJECTION_SYSTEM<TV>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const MPM_PROJECTION_VECTOR<TV>& vx=dynamic_cast<const MPM_PROJECTION_VECTOR<TV>&>(x),vy=dynamic_cast<const MPM_PROJECTION_VECTOR<TV>&>(y);
    return vx.v.Dot_Product_Double_Precision(vx.v,vy.v);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class TV> typename TV::SCALAR MPM_PROJECTION_SYSTEM<TV>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const MPM_PROJECTION_VECTOR<TV>& vx=dynamic_cast<const MPM_PROJECTION_VECTOR<TV>&>(x);
    return vx.v.Maximum_Magnitude();
}
//#####################################################################
// Function Project
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
} 
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class TV> void MPM_PROJECTION_SYSTEM<TV>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    const MPM_PROJECTION_VECTOR<TV>& vr=dynamic_cast<const MPM_PROJECTION_VECTOR<TV>&>(r);
    MPM_PROJECTION_VECTOR<TV>& vz=dynamic_cast<MPM_PROJECTION_VECTOR<TV>&>(z);
    temp_vector.Resize(A.m);
    A.C->Solve_Forward_Substitution(vr.v,temp_vector,true);
    A.C->Solve_Backward_Substitution(temp_vector,vz.v,false,true);
}
namespace PhysBAM{
template class MPM_PROJECTION_SYSTEM<VECTOR<float,3> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<float,2> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<float,1> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,3> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,2> >;
template class MPM_PROJECTION_SYSTEM<VECTOR<double,1> >;
}
