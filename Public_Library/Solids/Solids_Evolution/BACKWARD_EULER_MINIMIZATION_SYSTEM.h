//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_MINIMIZATION_SYSTEM
//#####################################################################
#ifndef __BACKWARD_EULER_MINIMIZATION_SYSTEM__
#define __BACKWARD_EULER_MINIMIZATION_SYSTEM__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
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
    bool use_l_inf_norm;
    KRYLOV_VECTOR_BASE<T>* tmp;

    struct COLLISION
    {
        int object;
        int p;
        T phi,n_dE;
        TV n,H_dE;
        SYMMETRIC_MATRIX<T,TV::m> H;
    };

    ARRAY<COLLISION> collisions;
    HASHTABLE<int,int> forced_collisions;
    EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities;

    BACKWARD_EULER_MINIMIZATION_SYSTEM(SOLID_BODY_COLLECTION<TV>& solid_body_collection,EXAMPLE_FORCES_AND_VELOCITIES<TV>* example_forces_and_velocities);

    virtual ~BACKWARD_EULER_MINIMIZATION_SYSTEM();

    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const override;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV0,const KRYLOV_VECTOR_BASE<T>& BV1) const override;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const override;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const override;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const override;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const override;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const override;
};
}
#endif
