//#####################################################################
// Copyright 2015, Andre Pradhana, Greg Klar, Chenfanfu Jiang, Chuyuan Fu
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_KKT_KRYLOV_SYSTEM
//#####################################################################
#ifndef __MPM_KKT_KRYLOV_SYSTEM__
#define __MPM_KKT_KRYLOV_SYSTEM__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Vectors/VECTOR.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
namespace PhysBAM{
template<class TV> class MPM_KKT_EXAMPLE;
template<class TV> class MPM_KKT_KRYLOV_VECTOR;

template<class TV>
class MPM_KKT_KRYLOV_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    const MPM_KKT_EXAMPLE<TV>& example;
    SPARSE_MATRIX_FLAT_MXN<T> D;
    SPARSE_MATRIX_FLAT_MXN<T> B;
    SPARSE_MATRIX_FLAT_MXN<T> M;
    MPM_KKT_KRYLOV_SYSTEM(const MPM_KKT_EXAMPLE<TV>& example);
    virtual ~MPM_KKT_KRYLOV_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const override;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const override;
    void Build_Div_Matrix();
    void Build_B_Matrix(bool& have_non_zero);
    void Build_FEM_Mass_Matrix();
};
}
#endif
