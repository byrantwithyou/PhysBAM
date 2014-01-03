//#####################################################################
// Copyright 2009, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Krylov_Solvers/MATRIX_SYSTEM.h>
#include <Tools/Matrices/MATRIX.h>
#include <Tools/Matrices/MATRIX_MXN.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Vectors/VECTOR.h>
using namespace PhysBAM;
//#####################################################################
// Constructor
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
MATRIX_SYSTEM(const T_MATRIX& A_input)
    :KRYLOV_SYSTEM_BASE<T>(false,true),A(A_input)
{
}
//#####################################################################
// Destructor
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
~MATRIX_SYSTEM()
{
}
//#####################################################################
// Function Multiply
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    VECTOR_T& vresult=dynamic_cast<VECTOR_T&>(result);
    A.Times(vx.v,vresult.v);
}
//#####################################################################
// Function Inner_Product
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> double MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x),vy=dynamic_cast<const VECTOR_T&>(y);
    return vx.v.Dot_Product_Double_Precision(vx.v,vy.v);
}
//#####################################################################
// Function Convergence_Norm
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> T MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const
{
    const VECTOR_T& vx=dynamic_cast<const VECTOR_T&>(x);
    return vx.v.Maximum_Magnitude();
}
//#####################################################################
// Function Project
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Project(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Set_Boundary_Conditions
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const
{
}
//#####################################################################
// Function Project_Nullspace
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const
{
} 
//#####################################################################
// Function Initialize_Preconditioner
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Set_Preconditioner(const T_MATRIX_PRECON& preconditioner,VECTOR_T& vector)
{
    this->use_preconditioner=true;
    P=&preconditioner;
    temp_vector=&vector;
}
//#####################################################################
// Function Apply_Preconditioner_Helper
//#####################################################################
template<class T,class VECTOR_T> static void
Apply_Preconditioner_Helper(SPARSE_MATRIX_FLAT_MXN<T>* P,const VECTOR_T& vr,VECTOR_T& vz,VECTOR_T* temp_vector)
{
    P->Solve_Forward_Substitution(vr.v,temp_vector->v,true);
    P->Solve_Backward_Substitution(temp_vector->v,vz.v,false,true);
}
//#####################################################################
// Function Apply_Preconditioner_Helper
//#####################################################################
template<class MATRIX_T,class VECTOR_T> static void
Apply_Preconditioner_Helper(MATRIX_T* P,const VECTOR_T& vr,VECTOR_T& vz,VECTOR_T* temp_vector)
{
}
//#####################################################################
// Function Apply_Preconditioner
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON> void MATRIX_SYSTEM<T_MATRIX,T,VECTOR_T,T_MATRIX_PRECON>::
Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const
{
    const VECTOR_T& vr=dynamic_cast<const VECTOR_T&>(r);
    VECTOR_T& vz=dynamic_cast<VECTOR_T&>(z);
    Apply_Preconditioner_Helper(P,vr,vz,temp_vector);
}
namespace PhysBAM{
template struct MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<float>,float,KRYLOV_VECTOR_WRAPPER<float,ARRAY<float> > >;
template struct MATRIX_SYSTEM<SPARSE_MATRIX_FLAT_MXN<double>,double,KRYLOV_VECTOR_WRAPPER<double,ARRAY<double> > >;
template struct MATRIX_SYSTEM<MATRIX<float,1,1>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,1> >,MATRIX<float,1,1> >;
template struct MATRIX_SYSTEM<MATRIX<float,2,2>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,2> >,MATRIX<float,2,2> >;
template struct MATRIX_SYSTEM<MATRIX<float,3,3>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,3> >,MATRIX<float,3,3> >;
template struct MATRIX_SYSTEM<MATRIX<float,4,4>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,4> >,MATRIX<float,4,4> >;
template struct MATRIX_SYSTEM<MATRIX<float,5,5>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,5> >,MATRIX<float,5,5> >;
template struct MATRIX_SYSTEM<MATRIX<float,6,6>,float,KRYLOV_VECTOR_WRAPPER<float,VECTOR<float,6> >,MATRIX<float,6,6> >;
template struct MATRIX_SYSTEM<MATRIX<double,1,1>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,1> >,MATRIX<double,1,1> >;
template struct MATRIX_SYSTEM<MATRIX<double,2,2>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,2> >,MATRIX<double,2,2> >;
template struct MATRIX_SYSTEM<MATRIX<double,3,3>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,3> >,MATRIX<double,3,3> >;
template struct MATRIX_SYSTEM<MATRIX<double,4,4>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,4> >,MATRIX<double,4,4> >;
template struct MATRIX_SYSTEM<MATRIX<double,5,5>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,5> >,MATRIX<double,5,5> >;
template struct MATRIX_SYSTEM<MATRIX<double,6,6>,double,KRYLOV_VECTOR_WRAPPER<double,VECTOR<double,6> >,MATRIX<double,6,6> >;
template struct MATRIX_SYSTEM<MATRIX_MXN<float>,float,KRYLOV_VECTOR_WRAPPER<float,ARRAY<float,int> >,MATRIX_MXN<float> >;
template struct MATRIX_SYSTEM<MATRIX_MXN<double>,double,KRYLOV_VECTOR_WRAPPER<double,ARRAY<double,int> >,MATRIX_MXN<double> >;
}
