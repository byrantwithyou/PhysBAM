//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_MINIMIZATION_OBJECTIVE
//#####################################################################
#ifndef __BACKWARD_EULER_MINIMIZATION_OBJECTIVE__
#define __BACKWARD_EULER_MINIMIZATION_OBJECTIVE__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class BACKWARD_EULER_MINIMIZATION_SYSTEM;
template<class TV> class BACKWARD_EULER_SYSTEM;

template<class TV>
class BACKWARD_EULER_MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    mutable GENERALIZED_VELOCITY<TV> v1;
    ARRAY<TV> X0;
    ARRAY<FRAME<TV> > frame0;
    GENERALIZED_VELOCITY<TV> &v0,&tmp0,&tmp1;
    BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system;
    T dt,time;

    BACKWARD_EULER_MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system);
    virtual ~BACKWARD_EULER_MINIMIZATION_OBJECTIVE();

    void Reset();
    void Compute(const KRYLOV_VECTOR_BASE<T>& dv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const PHYSBAM_OVERRIDE;
};
}
#endif