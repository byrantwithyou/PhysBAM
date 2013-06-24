//#####################################################################
// Copyright 2013, Chenfanfu Jiang
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_POISSON_SYSTEM
//#####################################################################
#ifndef __MPM_POISSON_SYSTEM__
#define __MPM_POISSON_SYSTEM__

#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Utilities/PHYSBAM_OVERRIDE.h>
namespace PhysBAM{

template<class T> class KRYLOV_VECTOR_BASE;
template<class TV> class MPM_PROJECTION;
template<class TV>
class MPM_POISSON_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;typedef VECTOR<int,TV::m> TV_INT;
    typedef KRYLOV_SYSTEM_BASE<T> BASE;
public:
    MPM_PROJECTION<TV>& proj;
    ARRAY<T,TV_INT> jacobi_scales;

    MPM_POISSON_SYSTEM(MPM_PROJECTION<TV>& proj);
    virtual ~MPM_POISSON_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_VECTOR_BASE<T>& result) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& x,const KRYLOV_VECTOR_BASE<T>& y) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& x) const PHYSBAM_OVERRIDE;
protected:
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& r,KRYLOV_VECTOR_BASE<T>& z) const PHYSBAM_OVERRIDE;
//#####################################################################
};
}
#endif
