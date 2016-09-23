//#####################################################################
// Copyright 2013.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class BACKWARD_EULER_MINIMIZATION_OBJECTIVE
//#####################################################################
#ifndef __BACKWARD_EULER_MINIMIZATION_OBJECTIVE__
#define __BACKWARD_EULER_MINIMIZATION_OBJECTIVE__
#include <Core/Arrays/ARRAY.h>
#include <Core/Data_Structures/HASHTABLE.h>
#include <Core/Data_Structures/PAIR.h>
#include <Tools/Krylov_Solvers/KRYLOV_SYSTEM_BASE.h>
#include <Tools/Krylov_Solvers/KRYLOV_VECTOR_BASE.h>
#include <Tools/Nonlinear_Equations/NONLINEAR_FUNCTION.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_MINIMIZATION_SYSTEM.h>
#include <Solids/Solids_Evolution/BACKWARD_EULER_SYSTEM.h>
namespace PhysBAM{
template<class TV> class SOLID_BODY_COLLECTION;
template<class TV> class BACKWARD_EULER_MINIMIZATION_SYSTEM;
template<class TV> class BACKWARD_EULER_SYSTEM;
template<class TV> class IMPLICIT_OBJECT;

template<class TV>
class BACKWARD_EULER_MINIMIZATION_OBJECTIVE:public NONLINEAR_FUNCTION<typename TV::SCALAR(KRYLOV_VECTOR_BASE<typename TV::SCALAR>&)>
{
    typedef typename TV::SCALAR T;
    typedef typename BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>::COLLISION COLLISION;
public:
    SOLID_BODY_COLLECTION<TV>& solid_body_collection;
    mutable GENERALIZED_VELOCITY<TV> v1;
    ARRAY<TV> X0;
    ARRAY<FRAME<TV> > frame0;
    GENERALIZED_VELOCITY<TV> &v0,&tmp0,&tmp1;
    BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system;
    T dt,time;
    ARRAY<IMPLICIT_OBJECT<TV>*> collision_objects;
    ARRAY<T> coefficient_of_friction;
    T collision_thickness;
    HASHTABLE<PAIR<int,int> > disabled_collision;
    HASHTABLE<int> solids_forces_lazy,rigids_forces_lazy,deformables_forces_lazy;
    mutable T last_energy;
    bool collisions_in_solve;

    BACKWARD_EULER_MINIMIZATION_OBJECTIVE(SOLID_BODY_COLLECTION<TV>& solid_body_collection,BACKWARD_EULER_MINIMIZATION_SYSTEM<TV>& minimization_system);
    virtual ~BACKWARD_EULER_MINIMIZATION_OBJECTIVE();

    void Reset();
    void Compute(const KRYLOV_VECTOR_BASE<T>& dv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const override;
    void Compute_Unconstrained(const KRYLOV_VECTOR_BASE<T>& dv,KRYLOV_SYSTEM_BASE<T>* h,KRYLOV_VECTOR_BASE<T>* g,T* e) const;
    void Initial_Guess(KRYLOV_VECTOR_BASE<T>& dv) const;
    void Adjust_For_Collision(KRYLOV_VECTOR_BASE<T>& Bdv) const;
    void Make_Feasible(KRYLOV_VECTOR_BASE<T>& dv) const override;
    void Project_Gradient_And_Prune_Constraints(KRYLOV_VECTOR_BASE<T>& dv,bool allow_sep) const;
    void Test_Diff(const KRYLOV_VECTOR_BASE<T>& dv);
    void Disable_Current_Colliding_Pairs(T thickness=0);
    T Update_Position_Based_State_Early_Out(T time,bool is_position_update,T energy_early_out,bool update_hessian) const;
};
}
#endif
