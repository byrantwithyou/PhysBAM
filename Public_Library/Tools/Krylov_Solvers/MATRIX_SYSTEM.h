//#####################################################################
// Copyright 2007-2009, Geoffrey Irving, Tamar Shinar.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Namespace MATRIX_SYSTEM
//#####################################################################
#ifndef __MATRIX_SYSTEM__
#define __MATRIX_SYSTEM__

#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Log/DEBUG_UTILITIES.h>
#include <Tools/Matrices/MATRIX_FORWARD.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
#include <Tools/Utilities/TYPE_UTILITIES.h>
#include <Tools/Vectors/VECTOR_FORWARD.h>
namespace PhysBAM{
//#####################################################################
// Struct IS_MATRIX
//#####################################################################
template<class T_MATRIX> struct IS_MATRIX {static const bool value=false;};
template<class T> struct IS_MATRIX<MATRIX_MXN<T> > {static const bool value=true;};
template<class T> struct IS_MATRIX<SPARSE_MATRIX_NXN<T> > {static const bool value=true;};
template<class T> struct IS_MATRIX<SPARSE_MATRIX_FLAT_MXN<T> > {static const bool value=true;};

//#####################################################################
// Class MATRIX_SYSTEM
//#####################################################################
template<class T_MATRIX,class T,class VECTOR_T,class T_MATRIX_PRECON=T_MATRIX>
struct MATRIX_SYSTEM:public KRYLOV_SYSTEM_BASE<T>
{
    const T_MATRIX& A;
    const T_MATRIX_PRECON* P;
    mutable VECTOR_T* temp_vector;

    MATRIX_SYSTEM(const T_MATRIX& A_input);
    virtual ~MATRIX_SYSTEM();

    void Set_Preconditioner(const T_MATRIX_PRECON& preconditioner,VECTOR_T& vector);
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE;
};
//#####################################################################
}
#endif
