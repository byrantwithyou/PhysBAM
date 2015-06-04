//#####################################################################
// Copyright 2015, Greg Klar, Andre Pradhana
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class FLUID_KRYLOV_SYSTEM
//#####################################################################
#ifndef __FLUID_KRYLOV_SYSTEM__
#define __FLUID_KRYLOV_SYSTEM__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Data_Structures/HASHTABLE.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Vectors/VECTOR.h>
namespace PhysBAM{
template<class TV> class MPM_EXAMPLE;
template<class TV> class FLUID_KRYLOV_VECTOR;

template<class TV>
class FLUID_KRYLOV_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
    typedef VECTOR<int,TV::m> TV_INT;
public:
    MPM_EXAMPLE<TV>& example;
    FLUID_KRYLOV_SYSTEM(MPM_EXAMPLE<TV>& example);
    virtual ~FLUID_KRYLOV_SYSTEM();

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