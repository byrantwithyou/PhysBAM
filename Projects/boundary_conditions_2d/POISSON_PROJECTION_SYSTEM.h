#ifndef __POISSON_PROJECTION_SYSTEM__
#define __POISSON_PROJECTION_SYSTEM__
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_WRAPPER.h>
#include <Tools/Matrices/SPARSE_MATRIX_FLAT_MXN.h>
#include <Tools/Vectors/VECTOR.h>
using namespace PhysBAM;

// TODO: Detect Neumann pockets.
// TODO: Preconditioner

template<class TV>
struct POISSON_PROJECTION_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef KRYLOV_VECTOR_WRAPPER<T,ARRAY<T> > VECTOR_T;

    SPARSE_MATRIX_FLAT_MXN<T> gradient,neg_divergence;
    SPARSE_MATRIX_FLAT_MXN<T> poisson;
    ARRAY<T> beta_inverse;
    ARRAY<ARRAY<T> > projections;
    mutable ARRAY<T> temp;

    POISSON_PROJECTION_SYSTEM();
    virtual ~POISSON_PROJECTION_SYSTEM();

    void Initialize();
    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const override;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const override;
};

#endif
