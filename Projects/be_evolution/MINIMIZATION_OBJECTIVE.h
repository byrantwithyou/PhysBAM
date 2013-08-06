//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MINIMIZATION_OBJECTIVE
//#####################################################################
#ifndef __MINIMIZATION_OBJECTIVE__
#define __MINIMIZATION_OBJECTIVE__
#include <Tools/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
using namespace PhysBAM;
namespace PhysBAM{template<class TV> class SOLID_BODY_COLLECTION;}
template<class TV>
class MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>, public KRYLOV_SYSTEM_BASE<typename TV::SCALAR>
{
    typedef typename TV::SCALAR T;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    T dt,time,a,b,c;
    ARRAY<TWIST<TV> > fake_rigid_x;
    GENERALIZED_VELOCITY<TV> x1,v1,&x0,&v0,&d;
    GENERALIZED_VELOCITY<TV> &tmp0,&tmp1;
    BACKWARD_EULER_SYSTEM<TV> system;

    MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,T dt,T time);
    virtual ~MINIMIZATION_OBJECTIVE();

    void Compute(const KRYLOV_VECTOR_BASE<T>& x,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const PHYSBAM_OVERRIDE;
    void Multiply(const KRYLOV_VECTOR_BASE<T>& BV,KRYLOV_VECTOR_BASE<T>& BF) const PHYSBAM_OVERRIDE;
    double Inner_Product(const KRYLOV_VECTOR_BASE<T>& BV1,const KRYLOV_VECTOR_BASE<T>& BV2) const PHYSBAM_OVERRIDE;
    T Convergence_Norm(const KRYLOV_VECTOR_BASE<T>& BR) const PHYSBAM_OVERRIDE;
    void Project(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Project_Nullspace(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
    void Apply_Preconditioner(const KRYLOV_VECTOR_BASE<T>& V,KRYLOV_VECTOR_BASE<T>& R) const PHYSBAM_OVERRIDE;
    void Set_Boundary_Conditions(KRYLOV_VECTOR_BASE<T>& V) const PHYSBAM_OVERRIDE;
};
#endif
