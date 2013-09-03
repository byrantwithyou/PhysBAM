//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_MINIMIZATION_SYSTEM
//#####################################################################
#ifndef __BACKWARD_EULER_MINIMIZATION_SYSTEM__
#define __BACKWARD_EULER_MINIMIZATION_SYSTEM__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class EXAMPLE_FORCES_AND_VELOCITIES;

template<class TV>
class BACKWARD_EULER_MINIMIZATION_SYSTEM:public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    T dt,time;
    KRYLOV_VECTOR_BASE<T>* tmp;

    struct COLLISION
    {
        int p;
        T phi,n_dE;
        TV n,H_dE;
        SYMMETRIC_MATRIX<T,TV::m> H;
    };

    ARRAY<COLLISION> collisions;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities;

    BACKWARD_EULER_MINIMIZATION_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities);

    virtual ~BACKWARD_EULER_MINIMIZATION_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
};
}
#endif
