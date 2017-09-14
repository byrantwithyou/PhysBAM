//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class MPM_OBJECTIVE
//#####################################################################
#ifndef __MPM_OBJECTIVE__
#define __MPM_OBJECTIVE__
#include <Core/Arrays/ARRAY.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Hybrid_Methods/System/MPM_KRYLOV_SYSTEM.h>
namespace PhysBAM{
template<class TV> class MPM_KRYLOV_SYSTEM;
template<class TV> class MPM_KRYLOV_VECTOR;

template<class TV>
class MPM_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>
{
    typedef typename TV::SCALAR T;
    typedef typename MPM_KRYLOV_SYSTEM<TV>::COLLISION COLLISION;
public:
    typedef typename MPM_COLLISION_OBJECT<TV>::COLLISION_TYPE COLLISION_TYPE;
    MPM_KRYLOV_SYSTEM<TV>& system;
    MPM_KRYLOV_VECTOR<TV> &v0,&v1,&tmp0,&tmp1,&tmp2;
    ARRAY<MATRIX<T,TV::m> > F0;
    ARRAY<SYMMETRIC_MATRIX<T,TV::m> > S0;
    ARRAY<TV> X0;
    T collision_thickness;

    MPM_OBJECTIVE(MPM_EXAMPLE<TV>& example);
    virtual ~MPM_OBJECTIVE();

    void Reset();
    void Compute(const KRYLOV_VECTOR_BASE<T>& dv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const override;
    void Compute_Unconstrained(const KRYLOV_VECTOR_BASE<T>& dv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const;
    bool Initial_Guess(KRYLOV_VECTOR_BASE<T>& dv,T tolerance,bool no_test) const;
    void Adjust_For_Collision(KRYLOV_VECTOR_BASE<T>& Bdv) const;
    void Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const override;
    void Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& dv,bool allow_sep) const;
    void Test_Diff(const KRYLOV_VECTOR_BASE<T>& dv);
    void Update_F(const MPM_KRYLOV_VECTOR<TV>& v) const;
    void Restore_F() const;
    bool Test_Add_Collision(MPM_KRYLOV_VECTOR<TV>& dv,int object,int p,
        bool prune,bool stick,ARRAY<COLLISION>& collisions,
        ARRAY<int>& stuck_nodes,ARRAY<TV>& stuck_velocity) const;
};
}
#endif
