//#####################################################################
// Copyright 2008, Geoffrey Irving.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class KRYLOV_SYSTEM_BASE
//#####################################################################
#ifndef __KRYLOV_SYSTEM_BASE__
#define __KRYLOV_SYSTEM_BASE__

#include <Core/Arrays/ARRAY.h>
namespace PhysBAM{

template<class T> class KRYLOV_VECTOR_BASE;
template<class T>
class KRYLOV_SYSTEM_BASE
{
public:
    bool use_preconditioner;
    bool preconditioner_commutes_with_projection;

    KRYLOV_SYSTEM_BASE(const bool use_preconditioner,const bool preconditioner_commutes_with_projection);
    virtual ~KRYLOV_SYSTEM_BASE();

    virtual void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result,bool transpose=false) const=0;

    virtual double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const=0;

    virtual T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const=0;

    // modifies system Ax=b to PAPx=Pb
    virtual void Project(KRYLOV_VECTOR_BASE<T>& x) const=0;

    // implements Sx=Px+x0, where x0 is the desired component of x in the nullspace of P
    virtual void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const=0;

    // removes component of x in nullspace of A (used to project residual for stopping conditions)
    virtual void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const;
    const KRYLOV_VECTOR_BASE<T>& Precondition(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const;
    void Test_System(const KRYLOV_VECTOR_BASE<T>& t,bool assume_symmetry=true) const;
    void Compute_Nullspace(const KRYLOV_VECTOR_BASE<T>& tmp,ARRAY<KRYLOV_VECTOR_BASE<T>*>& null,int max_null) const;
    void Compute_Small_Eigenvectors(const KRYLOV_VECTOR_BASE<T>& tmp,ARRAY<KRYLOV_VECTOR_BASE<T>*>& null,
        ARRAY<KRYLOV_VECTOR_BASE<T>*>& eigenvectors,ARRAY<T>& eigenvalues,int max_eigen,T tol,int power_iter) const;
    virtual void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const;
//#####################################################################
};
}
#endif
