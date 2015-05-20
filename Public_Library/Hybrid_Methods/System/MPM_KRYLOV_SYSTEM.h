//#####################################################################
// Copyright 2015, Craig Schroeder.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_KRYLOV_SYSTEM
//#####################################################################
#ifndef __MPM_KRYLOV_SYSTEM__
#define __MPM_KRYLOV_SYSTEM__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class MPM_EXAMPLE;
template<class TV> class MPM_KRYLOV_VECTOR;

template<class TV>
class MPM_KRYLOV_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    MPM_EXAMPLE<TV>& example;
    struct COLLISION
    {
        int object;
        int p;
        T phi,n_dE;
        TV n,H_dE;
        SYMMETRIC_MATRIX<T,TV::m> H;
        bool allow_sep;
    };

    MPM_KRYLOV_VECTOR<TV>& tmp;
    ARRAY<COLLISION> collisions;
    ARRAY<int> stuck_nodes;
    ARRAY<TV> stuck_velocity;
    HASHTABLE<int,int> forced_collisions;
    MPM_KRYLOV_SYSTEM(MPM_EXAMPLE<TV>& example);
    virtual ~MPM_KRYLOV_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
};
}
#endif
